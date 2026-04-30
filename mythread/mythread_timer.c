#include "mythread_timer.h"
#include <string.h>

static inline uint64_t atsc(void) {
#ifdef __ARM_ARCH
  uint64_t tsc;
  asm volatile("mrs %0, cntvct_el0" : "=r" (tsc));
  return tsc&0xFFFFFFFFFFFFFFL;
#else
  uint64_t a,d;
  asm volatile("rdtsc " : "=a" (a),"=d"(d));
  return (d<<32)|a;
#endif
}

void tscinit(void){
  memset(ti->tscds,0,sizeof(ti->tscds));
}

void tscb(int id){
  ti->tscbs[id]=atsc();
}

void tsce(int id){
  uint64_t d=atsc();
  ti->tscds[id]+=((d-ti->tscbs[id]));
}

void tsceb(int id){
  uint64_t d=atsc();
  ti->tscds[id-1]+=((d-ti->tscbs[id-1]));
  ti->tscbs[id]=d;
}

void tscb_(int *id_){
  int id=*id_;
  ti->tscbs[id]=atsc();
}

void tsce_(int *id_){
  int id=*id_;
  uint64_t d=atsc();
  ti->tscds[id]+=((d-ti->tscbs[id]));
}

void tsceb_(int *id_){
  int id=*id_;
  uint64_t d=atsc();
  ti->tscds[id-1]+=((d-ti->tscbs[id-1]));
  ti->tscbs[id]=d;
}

void prtsc(const char*tag){
  char buf[1024]={0};
  for(int i=0;i<16;i++){
    sprintf(buf+i*10,"|%9.5f     Z",ti->tscds[i]*(1./(100*1000*1000)));
  }
  fprintf(ti->fo,"%2s(s):%2.2d %2.2d %2d:%s\n",tag,mpi_id,ti->igrp,ti->ind,buf);
  fflush(ti->fo);
}

void prtsc_(void){
  prtsc("MT");
}
