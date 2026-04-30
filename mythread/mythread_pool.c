#include "mythread_pool.h"
#include "mythread_locv.h"
#include "mythread_util.h"
#include <string.h>

typedef struct _mt_taskpool{
  int capacity;
  int flags;
  int slot;
  pthread_mutex_t mu;
  pthread_cond_t cv_epoch;
  pthread_cond_t cv_nonempty;
  pthread_cond_t cv_nonfull;
  int shutdown;
  int open;
  unsigned int epoch;
  int head;
  int tail;
  int count;
  int nworkers;
  int done_count;
  int workers_exited;
  mt_task *buf;
} mt_taskpool;

static void mt_mu_unlock(void *p){
  pthread_mutex_unlock((pthread_mutex_t*)p);
}

static inline int mt_taskpool_typ(void){
  if(ThreadG && gi && ti && ti->igrp>=0) return 1;
  return 2;
}

static inline int mt_taskpool_nworkers(void){
  if(!ti) return 0;
  int n=ti->Nthreads-1;
  if(n<0) n=0;
  return n;
}

static inline int mt_taskpool_valid_slot(int slot){
  return slot>=0 && slot<=100;
}

static mt_taskpool *mt_taskpool_get(int slot){
  int typ=mt_taskpool_typ();
  if(!mt_taskpool_valid_slot(slot)) return NULL;
  return (mt_taskpool*)GetLocV(typ,slot,NULL);
}

int mt_taskpool_attach(int slot,int capacity,int flags){
  int typ=mt_taskpool_typ();
  if(!mt_taskpool_valid_slot(slot)) return -5;
  if(capacity<=0) return -1;
  if(GetLocV(typ,slot,NULL)) return -2;
  mt_taskpool *tp=(mt_taskpool*)hmalloc(sizeof(*tp));
  if(!tp) return -3;
  memset(tp,0,sizeof(*tp));
  tp->capacity=capacity;
  tp->flags=flags;
  tp->slot=slot;
  tp->buf=(mt_task*)hmalloc((size_t)capacity*sizeof(*tp->buf));
  if(!tp->buf){
    free(tp);
    return -4;
  }
  pthread_mutex_init(&tp->mu,NULL);
  pthread_cond_init(&tp->cv_epoch,NULL);
  pthread_cond_init(&tp->cv_nonempty,NULL);
  pthread_cond_init(&tp->cv_nonfull,NULL);
  tp->shutdown=0;
  tp->open=0;
  tp->epoch=0;
  tp->head=tp->tail=tp->count=0;
  tp->nworkers=mt_taskpool_nworkers();
  tp->done_count=0;
  SetLocV(typ,slot,tp);
  return 0;
}

int mt_taskpool_detach(int slot){
  int typ=mt_taskpool_typ();
  mt_taskpool *tp=(mt_taskpool*)GetLocV(typ,slot,NULL);
  if(!tp) return -1;
  pthread_mutex_lock(&tp->mu);
  if(!tp->shutdown || tp->open || tp->count>0 || tp->done_count<tp->nworkers){
    pthread_mutex_unlock(&tp->mu);
    return -2;
  }
  while(tp->workers_exited<tp->nworkers){
    pthread_cleanup_push(mt_mu_unlock,&tp->mu);
    pthread_cond_wait(&tp->cv_epoch,&tp->mu);
    pthread_cleanup_pop(0);
  }
  pthread_mutex_unlock(&tp->mu);
  SetLocV(typ,slot,NULL);
  pthread_cond_destroy(&tp->cv_nonfull);
  pthread_cond_destroy(&tp->cv_nonempty);
  pthread_cond_destroy(&tp->cv_epoch);
  pthread_mutex_destroy(&tp->mu);
  free(tp->buf);
  free(tp);
  return 0;
}

int mt_taskpool_begin(int slot){
  mt_taskpool *tp=mt_taskpool_get(slot);
  if(!tp) return -1;
  pthread_mutex_lock(&tp->mu);
  if(tp->shutdown){
    pthread_mutex_unlock(&tp->mu);
    return -2;
  }
  if(tp->open || tp->count!=0){
    pthread_mutex_unlock(&tp->mu);
    return -3;
  }
  tp->nworkers=mt_taskpool_nworkers();
  tp->done_count=0;
  tp->workers_exited=0;
  tp->head=tp->tail=tp->count=0;
  tp->open=1;
  tp->epoch++;
  pthread_cond_broadcast(&tp->cv_epoch);
  pthread_mutex_unlock(&tp->mu);
  return 0;
}

int mt_taskpool_submit(int slot,mt_task_fn fn,void *ctx){
  mt_taskpool *tp=mt_taskpool_get(slot);
  if(!tp) return -1;
  if(!fn) return -2;
  pthread_mutex_lock(&tp->mu);
  while(1){
    if(tp->shutdown || !tp->open){
      pthread_mutex_unlock(&tp->mu);
      return -3;
    }
    if(tp->count<tp->capacity) break;
    if(tp->flags & MT_TASKPOOL_TRY){
      pthread_mutex_unlock(&tp->mu);
      return -4;
    }
    if(tp->flags & MT_TASKPOOL_SPIN){
      pthread_mutex_unlock(&tp->mu);
      ntdelay(1);
      pthread_mutex_lock(&tp->mu);
      continue;
    }
    pthread_cleanup_push(mt_mu_unlock,&tp->mu);
    pthread_cond_wait(&tp->cv_nonfull,&tp->mu);
    pthread_cleanup_pop(0);
  }
  tp->buf[tp->tail].fn=fn;
  tp->buf[tp->tail].ctx=ctx;
  tp->tail++;
  if(tp->tail>=tp->capacity) tp->tail=0;
  tp->count++;
  pthread_cond_signal(&tp->cv_nonempty);
  pthread_mutex_unlock(&tp->mu);
  return 0;
}

int mt_taskpool_close(int slot){
  mt_taskpool *tp=mt_taskpool_get(slot);
  if(!tp) return -1;
  pthread_mutex_lock(&tp->mu);
  tp->open=0;
  pthread_cond_broadcast(&tp->cv_nonempty);
  pthread_cond_broadcast(&tp->cv_epoch);
  pthread_mutex_unlock(&tp->mu);
  return 0;
}

int mt_taskpool_wait(int slot){
  mt_taskpool *tp=mt_taskpool_get(slot);
  if(!tp) return -1;
  pthread_mutex_lock(&tp->mu);
  while(tp->done_count<tp->nworkers && !tp->shutdown){
    pthread_cleanup_push(mt_mu_unlock,&tp->mu);
    pthread_cond_wait(&tp->cv_epoch,&tp->mu);
    pthread_cleanup_pop(0);
  }
  pthread_mutex_unlock(&tp->mu);
  return 0;
}

int mt_taskpool_shutdown(int slot){
  mt_taskpool *tp=mt_taskpool_get(slot);
  if(!tp) return -1;
  pthread_mutex_lock(&tp->mu);
  tp->shutdown=1;
  tp->open=0;
  pthread_cond_broadcast(&tp->cv_nonfull);
  pthread_cond_broadcast(&tp->cv_nonempty);
  pthread_cond_broadcast(&tp->cv_epoch);
  pthread_mutex_unlock(&tp->mu);
  return 0;
}

int mt_taskpool_worker_loop(int slot){
  mt_taskpool *tp=mt_taskpool_get(slot);
  if(!tp) return -1;
  unsigned int local_epoch=0;
  for(;;){
    pthread_mutex_lock(&tp->mu);
    while(!tp->shutdown && tp->epoch==local_epoch){
      pthread_cleanup_push(mt_mu_unlock,&tp->mu);
      pthread_cond_wait(&tp->cv_epoch,&tp->mu);
      pthread_cleanup_pop(0);
    }
    if(tp->shutdown){
      tp->workers_exited++;
      pthread_cond_broadcast(&tp->cv_epoch);
      pthread_mutex_unlock(&tp->mu);
      break;
    }
    local_epoch=tp->epoch;
    for(;;){
      while(!tp->shutdown && tp->count==0 && tp->open){
        pthread_cleanup_push(mt_mu_unlock,&tp->mu);
        pthread_cond_wait(&tp->cv_nonempty,&tp->mu);
        pthread_cleanup_pop(0);
      }
      if(tp->shutdown) break;
      if(tp->count==0 && !tp->open) break;
      mt_task t=tp->buf[tp->head];
      tp->head++;
      if(tp->head>=tp->capacity) tp->head=0;
      tp->count--;
      pthread_cond_signal(&tp->cv_nonfull);
      pthread_mutex_unlock(&tp->mu);
      t.fn(t.ctx);
      pthread_mutex_lock(&tp->mu);
    }
    tp->done_count++;
    if(tp->done_count>=tp->nworkers){
      pthread_cond_broadcast(&tp->cv_epoch);
    } else {
      pthread_cond_signal(&tp->cv_epoch);
    }
    pthread_mutex_unlock(&tp->mu);
  }
  return 0;
}
