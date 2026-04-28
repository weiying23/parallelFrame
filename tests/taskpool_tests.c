#include <pthread.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>

#include "mythread/mythread.h"

#define SLOT_TP 0
#define INVALID_SLOT_LOW (-1)
#define INVALID_SLOT_HIGH 101

typedef enum {
  TEST_NONE = 0,
  TEST_GROUPED_BASIC,
  TEST_GROUPED_BLOCK_FULL,
  TEST_GROUPED_ZERO_TASK,
  TEST_GROUPED_ISOLATION,
  TEST_GROUPED_INVALID_ORDER,
  TEST_GROUPED_TRY_FULL,
  TEST_UNGROUPED_BASIC,
  TEST_BOUNDARY
} test_kind;

typedef struct {
  test_kind kind;
  const char *name;
  int grouped;
  int groups;
  int threads_per_group;
  int rounds;
  int tasks_per_group;
  int capacity;
  int flags;
  int task_sleep_us;
  int release_delay_us;
} test_config;

typedef struct {
  int global_index;
  int expected_group;
  int round_id;
  int task_id;
  int sleep_us;
} task_ctx;

typedef struct {
  volatile int *gate;
  int delay_us;
} gate_release_args;

static test_config g_cfg;
static task_ctx *g_tasks = NULL;
static int *g_expected_counts = NULL;
static int *g_exec_counts = NULL;
static volatile int g_failures = 0;
static volatile int g_total_expected = 0;
static volatile int g_total_executed = 0;
static volatile int g_worker_gates[MCLUST];

static void test_fail(const char *fmt,...){
  va_list ap;
  va_start(ap,fmt);
  fprintf(stderr,"[%s] FAIL: ",g_cfg.name?g_cfg.name:"taskpool");
  vfprintf(stderr,fmt,ap);
  fprintf(stderr,"\n");
  va_end(ap);
  __sync_add_and_fetch(&g_failures,1);
}

static void expect_int(const char *label,int got,int want){
  if(got!=want){
    test_fail("%s: got %d, want %d",label,got,want);
  }
}

static void expect_true(const char *label,int cond){
  if(!cond){
    test_fail("%s",label);
  }
}

static int total_task_slots(void){
  int total=g_cfg.groups*g_cfg.rounds*g_cfg.tasks_per_group;
  if(total<1) total=1;
  return total;
}

static int task_index(int round_id,int group_id,int task_id){
  return ((round_id*g_cfg.groups)+group_id)*g_cfg.tasks_per_group+task_id;
}

static double wall_seconds(void){
  struct timeval tv;
  gettimeofday(&tv,NULL);
  return (double)tv.tv_sec+(double)tv.tv_usec*1e-6;
}

static void *release_gate_after_delay(void *arg){
  gate_release_args *cfg=(gate_release_args*)arg;
  usleep((useconds_t)cfg->delay_us);
  *(cfg->gate)=1;
  return NULL;
}

static void mark_expected(int round_id,int group_id,int task_id){
  int idx=task_index(round_id,group_id,task_id);
  g_expected_counts[idx]++;
  __sync_add_and_fetch(&g_total_expected,1);
}

static void task_count_fn(void *arg){
  task_ctx *ctx=(task_ctx*)arg;
  if(g_cfg.grouped){
    if(!ti){
      test_fail("task %d ran without TLS thread info",ctx->global_index);
    }else{
      if(ti->igrp!=ctx->expected_group){
        test_fail("task %d ran in group %d, expected %d",
                  ctx->global_index,ti->igrp,ctx->expected_group);
      }
      if(ti->ind==0){
        test_fail("task %d ran on group main thread",ctx->global_index);
      }
    }
  }else{
    if(ti==&md.tm){
      test_fail("task %d ran on the non-grouped main thread",ctx->global_index);
    }
  }
  if(ctx->sleep_us>0){
    usleep((useconds_t)ctx->sleep_us);
  }
  int seen=__sync_add_and_fetch(&g_exec_counts[ctx->global_index],1);
  if(seen!=1){
    test_fail("task %d executed %d times",ctx->global_index,seen);
  }
  __sync_add_and_fetch(&g_total_executed,1);
}

static void submit_expected_task(int round_id,int group_id,int task_id){
  int idx=task_index(round_id,group_id,task_id);
  int rc=mt_taskpool_submit(SLOT_TP,task_count_fn,&g_tasks[idx]);
  if(rc!=0){
    test_fail("submit task[%d,%d,%d] failed: %d",round_id,group_id,task_id,rc);
    return;
  }
  mark_expected(round_id,group_id,task_id);
}

static void run_standard_round(int round_id,int group_id){
  int rc=mt_taskpool_begin(SLOT_TP);
  expect_int("mt_taskpool_begin",rc,0);
  if(rc!=0) return;

  if(!(g_cfg.kind==TEST_GROUPED_ZERO_TASK && round_id==1 && group_id==0)){
    for(int i=0;i<g_cfg.tasks_per_group;i++){
      submit_expected_task(round_id,group_id,i);
    }
  }

  rc=mt_taskpool_close(SLOT_TP);
  expect_int("mt_taskpool_close",rc,0);
  rc=mt_taskpool_wait(SLOT_TP);
  expect_int("mt_taskpool_wait",rc,0);
}

static void run_block_full_round(int group_id){
  pthread_t releaser;
  gate_release_args release_cfg={&g_worker_gates[group_id],g_cfg.release_delay_us};
  int releaser_started=0;

  g_worker_gates[group_id]=0;
  expect_int("mt_taskpool_begin",mt_taskpool_begin(SLOT_TP),0);
  submit_expected_task(0,group_id,0);

  if(pthread_create(&releaser,NULL,release_gate_after_delay,&release_cfg)==0){
    releaser_started=1;
  }else{
    test_fail("failed to create gate releaser thread");
    g_worker_gates[group_id]=1;
  }

  double t0=wall_seconds();
  int rc=mt_taskpool_submit(SLOT_TP,task_count_fn,&g_tasks[task_index(0,group_id,1)]);
  double dt=wall_seconds()-t0;
  expect_int("second submit under BLOCK",rc,0);
  if(rc==0){
    mark_expected(0,group_id,1);
  }
  if(releaser_started){
    pthread_join(releaser,NULL);
  }
  if(dt<((double)g_cfg.release_delay_us*0.5e-6)){
    test_fail("BLOCK submit returned too quickly: %.6f sec",dt);
  }

  rc=mt_taskpool_close(SLOT_TP);
  expect_int("mt_taskpool_close",rc,0);
  rc=mt_taskpool_wait(SLOT_TP);
  expect_int("mt_taskpool_wait",rc,0);
}

static void run_try_full_round(int group_id){
  g_worker_gates[group_id]=0;
  expect_int("mt_taskpool_begin",mt_taskpool_begin(SLOT_TP),0);
  submit_expected_task(0,group_id,0);

  int rc=mt_taskpool_submit(SLOT_TP,task_count_fn,&g_tasks[task_index(0,group_id,1)]);
  expect_int("second submit under TRY",rc,-4);

  g_worker_gates[group_id]=1;
  rc=mt_taskpool_close(SLOT_TP);
  expect_int("mt_taskpool_close",rc,0);
  rc=mt_taskpool_wait(SLOT_TP);
  expect_int("mt_taskpool_wait",rc,0);
}

static void run_invalid_order_checks(void){
  int rc=mt_taskpool_attach(SLOT_TP,g_cfg.capacity,g_cfg.flags);
  expect_int("initial attach",rc,0);

  rc=mt_taskpool_attach(SLOT_TP,g_cfg.capacity,g_cfg.flags);
  expect_int("duplicate attach",rc,-2);

  rc=mt_taskpool_submit(SLOT_TP,task_count_fn,&g_tasks[0]);
  expect_int("submit before begin",rc,-3);

  rc=mt_taskpool_begin(SLOT_TP);
  expect_int("first begin",rc,0);

  rc=mt_taskpool_begin(SLOT_TP);
  expect_int("second begin",rc,-3);

  rc=mt_taskpool_close(SLOT_TP);
  expect_int("close after begin",rc,0);

  rc=mt_taskpool_submit(SLOT_TP,task_count_fn,&g_tasks[0]);
  expect_int("submit after close",rc,-3);

  rc=mt_taskpool_detach(SLOT_TP);
  expect_int("detach before shutdown",rc,-2);

  rc=mt_taskpool_wait(SLOT_TP);
  expect_int("wait after close",rc,0);

  rc=mt_taskpool_shutdown(SLOT_TP);
  expect_int("shutdown",rc,0);

  rc=mt_taskpool_detach(SLOT_TP);
  expect_int("detach after shutdown",rc,0);
}

static void run_boundary_checks(void){
  expect_true("attach should reject negative slot",
              mt_taskpool_attach(INVALID_SLOT_LOW,4,MT_TASKPOOL_BLOCK)!=0);
  expect_true("attach should reject slot > 100",
              mt_taskpool_attach(INVALID_SLOT_HIGH,4,MT_TASKPOOL_BLOCK)!=0);
  expect_true("attach should reject zero capacity",
              mt_taskpool_attach(SLOT_TP,0,MT_TASKPOOL_BLOCK)!=0);

  expect_int("valid attach",mt_taskpool_attach(SLOT_TP,4,MT_TASKPOOL_BLOCK),0);
  expect_int("begin on empty pool",mt_taskpool_begin(SLOT_TP),0);
  expect_int("close on empty pool",mt_taskpool_close(SLOT_TP),0);
  expect_int("wait on empty pool",mt_taskpool_wait(SLOT_TP),0);
  expect_int("shutdown on empty pool",mt_taskpool_shutdown(SLOT_TP),0);
  expect_int("detach on empty pool",mt_taskpool_detach(SLOT_TP),0);
}

static void grouped_main_main(void){
  for(int round_id=0;round_id<g_cfg.rounds;round_id++){
    int state=round_id+1;
    mSetGrps(state);
    mWaitGrps(state);
  }
}

static void grouped_group_main(void){
  if(g_cfg.kind==TEST_GROUPED_INVALID_ORDER){
    gWaitMain(1);
    run_invalid_order_checks();
    gSetMain(1);
    return;
  }

  expect_int("attach",mt_taskpool_attach(SLOT_TP,g_cfg.capacity,g_cfg.flags),0);

  for(int round_id=0;round_id<g_cfg.rounds;round_id++){
    int state=round_id+1;
    gWaitMain(state);

    switch(g_cfg.kind){
      case TEST_GROUPED_BLOCK_FULL:
        run_block_full_round(ti->igrp);
        break;
      case TEST_GROUPED_TRY_FULL:
        run_try_full_round(ti->igrp);
        break;
      case TEST_GROUPED_BASIC:
      case TEST_GROUPED_ZERO_TASK:
      case TEST_GROUPED_ISOLATION:
        run_standard_round(round_id,ti->igrp);
        break;
      default:
        test_fail("unexpected grouped test kind %d",g_cfg.kind);
        break;
    }

    gSetMain(state);
  }

  expect_int("shutdown",mt_taskpool_shutdown(SLOT_TP),0);
  expect_int("detach",mt_taskpool_detach(SLOT_TP),0);
}

static void grouped_worker(void){
  while(!GetLocV(1,SLOT_TP,NULL)){
    ntdelay(1);
  }
  if(g_cfg.kind==TEST_GROUPED_BLOCK_FULL || g_cfg.kind==TEST_GROUPED_TRY_FULL){
    while(!g_worker_gates[ti->igrp]){
      ntdelay(1);
    }
  }
  expect_int("worker loop",mt_taskpool_worker_loop(SLOT_TP),0);
}

static void ungrouped_main(void){
  if(g_cfg.kind==TEST_BOUNDARY){
    run_boundary_checks();
    return;
  }

  expect_int("attach",mt_taskpool_attach(SLOT_TP,g_cfg.capacity,g_cfg.flags),0);
  for(int round_id=0;round_id<g_cfg.rounds;round_id++){
    run_standard_round(round_id,0);
  }
  expect_int("shutdown",mt_taskpool_shutdown(SLOT_TP),0);
  expect_int("detach",mt_taskpool_detach(SLOT_TP),0);
}

static void ungrouped_worker(void){
  while(!GetLocV(2,SLOT_TP,NULL)){
    ntdelay(1);
  }
  expect_int("worker loop",mt_taskpool_worker_loop(SLOT_TP),0);
}

void thread_run(){
  if(g_cfg.grouped){
    if(ti->igrp==-1){
      grouped_main_main();
    }else if(ti->ind==0){
      grouped_group_main();
    }else{
      grouped_worker();
    }
    return;
  }

  if(ti==&md.tm){
    ungrouped_main();
  }else{
    ungrouped_worker();
  }
}

static void init_task_arrays(void){
  int total=total_task_slots();
  memset((void*)g_worker_gates,0,sizeof(g_worker_gates));
  g_failures=0;
  g_total_expected=0;
  g_total_executed=0;
  g_tasks=(task_ctx*)calloc((size_t)total,sizeof(*g_tasks));
  g_expected_counts=(int*)calloc((size_t)total,sizeof(*g_expected_counts));
  g_exec_counts=(int*)calloc((size_t)total,sizeof(*g_exec_counts));
  if(!g_tasks || !g_expected_counts || !g_exec_counts){
    fprintf(stderr,"[%s] failed to allocate test buffers\n",g_cfg.name);
    exit(2);
  }

  for(int round_id=0;round_id<g_cfg.rounds;round_id++){
    for(int group_id=0;group_id<g_cfg.groups;group_id++){
      for(int task_id=0;task_id<g_cfg.tasks_per_group;task_id++){
        int idx=task_index(round_id,group_id,task_id);
        g_tasks[idx].global_index=idx;
        g_tasks[idx].expected_group=group_id;
        g_tasks[idx].round_id=round_id;
        g_tasks[idx].task_id=task_id;
        g_tasks[idx].sleep_us=g_cfg.task_sleep_us;
      }
    }
  }
}

static void free_task_arrays(void){
  free(g_exec_counts);
  free(g_expected_counts);
  free(g_tasks);
  g_exec_counts=NULL;
  g_expected_counts=NULL;
  g_tasks=NULL;
}

static void verify_counts(void){
  int total=total_task_slots();
  for(int i=0;i<total;i++){
    if(g_exec_counts[i]!=g_expected_counts[i]){
      test_fail("task slot %d executed %d times, expected %d",
                i,g_exec_counts[i],g_expected_counts[i]);
    }
  }
  if(g_total_executed!=g_total_expected){
    test_fail("total executed %d, expected %d",
              (int)g_total_executed,(int)g_total_expected);
  }
}

static int configure_case(const char *name){
  memset(&g_cfg,0,sizeof(g_cfg));

  if(strcmp(name,"grouped-basic")==0){
    g_cfg.kind=TEST_GROUPED_BASIC;
    g_cfg.name=name;
    g_cfg.grouped=1;
    g_cfg.groups=2;
    g_cfg.threads_per_group=5;
    g_cfg.rounds=4;
    g_cfg.tasks_per_group=24;
    g_cfg.capacity=8;
    g_cfg.flags=MT_TASKPOOL_BLOCK;
    return 0;
  }
  if(strcmp(name,"grouped-block-full")==0){
    g_cfg.kind=TEST_GROUPED_BLOCK_FULL;
    g_cfg.name=name;
    g_cfg.grouped=1;
    g_cfg.groups=2;
    g_cfg.threads_per_group=2;
    g_cfg.rounds=1;
    g_cfg.tasks_per_group=2;
    g_cfg.capacity=1;
    g_cfg.flags=MT_TASKPOOL_BLOCK;
    g_cfg.task_sleep_us=1000;
    g_cfg.release_delay_us=50000;
    return 0;
  }
  if(strcmp(name,"grouped-zero-task")==0){
    g_cfg.kind=TEST_GROUPED_ZERO_TASK;
    g_cfg.name=name;
    g_cfg.grouped=1;
    g_cfg.groups=2;
    g_cfg.threads_per_group=4;
    g_cfg.rounds=3;
    g_cfg.tasks_per_group=12;
    g_cfg.capacity=4;
    g_cfg.flags=MT_TASKPOOL_BLOCK;
    return 0;
  }
  if(strcmp(name,"grouped-isolation")==0){
    g_cfg.kind=TEST_GROUPED_ISOLATION;
    g_cfg.name=name;
    g_cfg.grouped=1;
    g_cfg.groups=3;
    g_cfg.threads_per_group=4;
    g_cfg.rounds=2;
    g_cfg.tasks_per_group=18;
    g_cfg.capacity=6;
    g_cfg.flags=MT_TASKPOOL_BLOCK;
    return 0;
  }
  if(strcmp(name,"grouped-invalid-order")==0){
    g_cfg.kind=TEST_GROUPED_INVALID_ORDER;
    g_cfg.name=name;
    g_cfg.grouped=1;
    g_cfg.groups=2;
    g_cfg.threads_per_group=1;
    g_cfg.rounds=1;
    g_cfg.tasks_per_group=1;
    g_cfg.capacity=4;
    g_cfg.flags=MT_TASKPOOL_BLOCK;
    return 0;
  }
  if(strcmp(name,"grouped-try-full")==0){
    g_cfg.kind=TEST_GROUPED_TRY_FULL;
    g_cfg.name=name;
    g_cfg.grouped=1;
    g_cfg.groups=2;
    g_cfg.threads_per_group=2;
    g_cfg.rounds=1;
    g_cfg.tasks_per_group=2;
    g_cfg.capacity=1;
    g_cfg.flags=MT_TASKPOOL_TRY;
    g_cfg.task_sleep_us=1000;
    return 0;
  }
  if(strcmp(name,"ungrouped-basic")==0){
    g_cfg.kind=TEST_UNGROUPED_BASIC;
    g_cfg.name=name;
    g_cfg.grouped=0;
    g_cfg.groups=1;
    g_cfg.threads_per_group=5;
    g_cfg.rounds=3;
    g_cfg.tasks_per_group=20;
    g_cfg.capacity=8;
    g_cfg.flags=MT_TASKPOOL_BLOCK;
    return 0;
  }
  if(strcmp(name,"boundary")==0){
    g_cfg.kind=TEST_BOUNDARY;
    g_cfg.name=name;
    g_cfg.grouped=0;
    g_cfg.groups=1;
    g_cfg.threads_per_group=1;
    g_cfg.rounds=1;
    g_cfg.tasks_per_group=1;
    g_cfg.capacity=4;
    g_cfg.flags=MT_TASKPOOL_BLOCK;
    return 0;
  }
  return -1;
}

static void print_usage(const char *argv0){
  fprintf(stderr,"Usage: %s <case>\n",argv0);
  fprintf(stderr,"Cases:\n");
  fprintf(stderr,"  grouped-basic\n");
  fprintf(stderr,"  grouped-block-full\n");
  fprintf(stderr,"  grouped-zero-task\n");
  fprintf(stderr,"  grouped-isolation\n");
  fprintf(stderr,"  grouped-invalid-order\n");
  fprintf(stderr,"  grouped-try-full\n");
  fprintf(stderr,"  ungrouped-basic\n");
  fprintf(stderr,"  boundary\n");
}

static int run_case(void){
  int manage_core_id=-1;
  int err=InitThreads(0,16,1,16,g_cfg.threads_per_group,g_cfg.groups,1,&manage_core_id);
  if(err!=0){
    fprintf(stderr,"[%s] InitThreads failed: %d\n",g_cfg.name,err);
    return 2;
  }

  StartThreads(thread_run);
  thread_run();
  EndThreads();
  return 0;
}

int main(int argc,char **argv){
  if(argc!=2){
    print_usage(argv[0]);
    return 2;
  }
  if(configure_case(argv[1])!=0){
    print_usage(argv[0]);
    return 2;
  }

  init_task_arrays();
  int rc=run_case();
  verify_counts();

  if(rc==0 && g_failures==0){
    printf("[PASS] %s\n",g_cfg.name);
  }else if(rc==0){
    rc=1;
  }

  free_task_arrays();
  return rc;
}
