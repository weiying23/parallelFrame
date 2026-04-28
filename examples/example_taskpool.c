#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mythread/mythread.h"

#define N_GROUPS 2
#define N_THREADS_PER_GROUP 5
#define SLOT_TP 0
#define TP_CAP 64
#define N_PER_GROUP 2000
#define N_ROUNDS 1

typedef struct {
  double *a;
  double *b;
  int i;
} TaskCtx;

static double *g_a=NULL;
static double *g_b=NULL;
static TaskCtx *g_ctx=NULL;

static void task_fn(void *p){
  TaskCtx *c=(TaskCtx*)p;
  int i=c->i;
  c->b[i]=c->a[i]*c->a[i];
}

static void worker_thread(){
  while(!GetLocV(1,SLOT_TP,NULL)){
    ntdelay(1);
  }
  mt_taskpool_worker_loop(SLOT_TP);
}

static void group_main(){
  int gid=ti->igrp;
  int base=gid*N_PER_GROUP;
  mt_taskpool_attach(SLOT_TP,TP_CAP,MT_TASKPOOL_BLOCK);
  for(int r=0;r<N_ROUNDS;r++){
    int s=r+1;
    gWaitMain(s);
    mt_taskpool_begin(SLOT_TP);
    if(!(r==1 && gid==0)){
      for(int i=0;i<N_PER_GROUP;i++){
        mt_taskpool_submit(SLOT_TP,task_fn,&g_ctx[base+i]);
      }
    }
    mt_taskpool_close(SLOT_TP);
    mt_taskpool_wait(SLOT_TP);
    gSetMain(s);
  }
  mt_taskpool_shutdown(SLOT_TP);
  mt_taskpool_detach(SLOT_TP);
}

static void main_main(){
  for(int r=0;r<N_ROUNDS;r++){
    int s=r+1;
    mSetGrps(s);
    mWaitGrps(s);
  }
}

void thread_run(){
  if(ti->igrp==-1){
    main_main();
  }else if(ti->ind==0){
    group_main();
  }else{
    worker_thread();
  }
}

int main(int argc,char **argv){
  (void)argc;
  (void)argv;
  int n_total=N_GROUPS*N_PER_GROUP;
  g_a=(double*)malloc((size_t)n_total*sizeof(*g_a));
  g_b=(double*)malloc((size_t)n_total*sizeof(*g_b));
  g_ctx=(TaskCtx*)malloc((size_t)n_total*sizeof(*g_ctx));
  if(!g_a||!g_b||!g_ctx) return 2;
  for(int i=0;i<n_total;i++){
    g_a[i]=sin((double)i*0.001);
    g_b[i]=0.0;
    g_ctx[i].a=g_a;
    g_ctx[i].b=g_b;
    g_ctx[i].i=i;
  }

  int ManageCoreId=-1;
  int err=InitThreads(0,16,1,16,N_THREADS_PER_GROUP,N_GROUPS,1,&ManageCoreId);
  if(err){
    printf("InitThreads err=%d\n",err);
    return 1;
  }
  StartThreads(thread_run);
  thread_run();
  EndThreads();

  int bad=0;
  for(int i=0;i<n_total;i++){
    double v=g_a[i]*g_a[i];
    if(fabs(g_b[i]-v)>1e-12){
      bad++;
      if(bad<10){
        printf("bad i=%d b=%f v=%f\n",i,g_b[i],v);
      }
    }
  }
  printf("check bad=%d/%d\n",bad,n_total);
  free(g_ctx);
  free(g_b);
  free(g_a);
  return bad?3:0;
}
