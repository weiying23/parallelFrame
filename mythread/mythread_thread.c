#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include "mythread_thread.h"
#include "mythread_util.h"
#include <pthread.h>
#include <string.h>
#include <unistd.h>

#ifdef __APPLE__
#include <sys/types.h>
#include <sys/sysctl.h>
#else
#include <sched.h>
#endif

/* ── TLS globals ── */
__thread THREADINFO *ti=NULL;
__thread threadGroup *gi=NULL;
__thread void*td=NULL;
__thread void*gd=NULL;

/* ── global config ── */
int ThreadG=-1;
int mpi_id,NCorePClu=38,NCluPNode=16,NCorePGrp=38,NGrpPProc=4,NThPGrp=7,NProcPNode=16,ManageCoreId=-1,NThreads=3;
int sync_flag_init = 6666;
threadProc md={0};

/* ── static forward decls ── */
static void zStartThreads(TFunc tfun,void*para,int detach,int clear);
static void InitThread(HTHREADINFO pti,TFunc tfun,void*para,int detach,int nlocv);
static void ClearThread(HTHREADINFO pti);

#define VFREE(p) if(p)free(p)

/* ── weak symbols ── */
__attribute__((weak))
void thread_run(void) {
  printf("Warning: using default thread_run(), please define your own.\n");
}

__attribute__((weak)) int _getgdsize_(void){ return sizeof(float); }
__attribute__((weak)) int _gettdsize_(void){ return sizeof(float); }

/* ── CPU ── */
int getcpuid(void){
#ifdef __APPLE__
  int cpu = 0;
  size_t size = sizeof(cpu);
  if (sysctlbyname("hw.logicalcpu", &cpu, &size, NULL, 0) == -1) {
    printf("warning: could not get CPU affinity , continuing...\n");
    return -1;
  }
  return cpu;
#else
  cpu_set_t mask;
  CPU_ZERO(&mask);
  if (sched_getaffinity(0, sizeof(mask), &mask) == -1){
    printf("warning: could not get CPU affinity , continuing...\n");
    return -1;
  }
  for(int i=0;i<CPU_SETSIZE;i++){
    if(CPU_ISSET(i,&mask)){
      return i;
    }
  }
  return 0;
#endif
}

int bindcpu(int id){
#ifdef __APPLE__
  return 0;
#else
  cpu_set_t mask;
  CPU_ZERO(&mask);
  CPU_SET(id,&mask);
  if (sched_setaffinity(0, sizeof(mask), &mask) == -1){
    printf("warning: could not set CPU affinity, continuing...\n");
    return -1;
  }
  return 0;
#endif
}

void bindthread(void){
  if(ti) bindcpu(ti->indg);
}

/* ── thread info ── */
int gettid(void){ return ti->ind; }
int getnt(void){ return ti->Nthreads; }

void _threadmain_(HTHREADINFO ti){
  bindcpu(ti->indg);
  ti->tfun();
}

/* ── config ── */
void setthread(int NCorePClu_,int NThPGrp_,int NGrpPProc_,int NProcPNode_,int ManageCoreId_){
  if(NCorePClu_>0){
    NCorePClu=NCorePClu_;
    NCorePGrp=NCorePClu_;
  }
  if(NThPGrp_>0){
    NThPGrp=NThPGrp_;
  }
  if(NGrpPProc_>0){
    NGrpPProc=NGrpPProc_;
  }
  if(NProcPNode_>0){
    NProcPNode=NProcPNode_;
  }
  ManageCoreId=ManageCoreId_;
  ThreadG=(NGrpPProc>1);
}

/* ── init ── */
int InitThreads(int mpi_id_,int NCorePClu_,int NCluPNode_,int NCorePGrp_,int NThPGrp_,int NGrpPProc_,int NProcPNode_,int *ManageCoreId_){
  mpi_id      =mpi_id_;
  md.mpi_id   =mpi_id_;
#define VD(d,v) if(v>0)d=v;
  VD(NCorePClu   ,NCorePClu_);
  VD(NCluPNode   ,NCluPNode_);
  VD(NCorePGrp   ,NCorePGrp_);
  VD(NThPGrp     ,NThPGrp_);
  VD(NGrpPProc   ,NGrpPProc_);
  VD(NProcPNode  ,NProcPNode_);
#undef VD
  ManageCoreId=*ManageCoreId_;

  if(NGrpPProc>1)ThreadG=1;
  else {ThreadG=0;NGrpPProc=1;}
  if(NThPGrp>MCOREPC){
    printf("NThPGrp %d biger than surported MCOREPC %d\n",NThPGrp,MCOREPC);
    return (-1);
  }
  if(NThPGrp>NCorePClu){
    printf("NThPGrp %d biger than NCorePClu %d\n",NThPGrp,NCorePClu);
    return (-2);
  }
  if(NGrpPProc>MCLUST){
    printf("NGrpPProc %d biger than surported MCLUST %d\n",NGrpPProc,MCLUST);
    return (-3);
  }
  if(NCluPNode>MCLUST){
    printf("NCluPNode %d biger than surported MCLUST %d\n",NCluPNode,MCLUST);
    return (-4);
  }
  if(NProcPNode > MCLUST){
    printf("NProcPNode %d biger than surported MCLUST %d\n", NProcPNode , MCLUST);
    return (-5);
  }
  initmd();
  NThreads=ThreadG?1:md.Nthreads;
  if(ThreadG){
    for(int i=0;i<md.ngrp;i++){
      NThreads+=md.grps[i]->Nthreads;
    }
  }
  *ManageCoreId_=ManageCoreId;
  return 0;
}

#define NGG 8

void initmd(void){
  int i,j,cbase;
  if(md.PThreadInited) return;
  memset(&md,0,sizeof(md));
  md.PThreadInited=1;
  md.state=0;
  ti=&md.tm;
  if(ThreadG){
    md.ngrp=NGrpPProc;
    int mg=0,mt;
    cbase=0;
    if(NProcPNode<1)NProcPNode=1;
    int ipr=mpi_id%NProcPNode;
    cbase=ipr*NCorePClu*NGrpPProc;
    if(ManageCoreId <-2)ManageCoreId =0;
    else if(ManageCoreId >0&&ManageCoreId<NThPGrp)ManageCoreId =0;
    else if(ManageCoreId>NCorePClu){
      if(ManageCoreId<NCorePClu*NCluPNode){
        ManageCoreId =0;
      }
    }
    if(ManageCoreId ==0){
      ManageCoreId = cbase;
      mg=0;mt=0;
    }else if(ManageCoreId ==-1){
      ManageCoreId = cbase;
      mg=0;mt=0;
    }else if(ManageCoreId ==-2){
      mg=-1;
      mt=-1;
      ManageCoreId = (ipr+NCorePClu*NCluPNode);
    }else{
      if(ManageCoreId<NCorePClu){
        mg=0;mt=ManageCoreId;
        ManageCoreId+=cbase;
      }else{
        mg=-1;mt=-1;
        ManageCoreId = (ipr+NCorePClu*NCluPNode);
      }
    }
    for(i=0;i<NGrpPProc;i++){
      bindcpu(cbase+i*NCorePClu);
      threadGroup *pgi=md.grps[i]=HNEW(threadGroup);
      memset(pgi,0,sizeof(threadGroup));
      pgi->pmd=&md;
      pgi->ngrp=NGrpPProc;
      pgi->Nthreads=NThPGrp;
      pgi->idMainThread=0;
      int ib=0;
      if(mg==i || mt==-1){
        if(!mt){ ib++; pgi->Nthreads--; }
        gi=ti->pg=pgi;
      }
      pgi->igrp=i;
      pgi->NSpecial=0;
      int nsize=_getgdsize_();
      if(nsize)pgi->gd=hmalloc(nsize);
      for(j=0;j<pgi->Nthreads;j++){
        struct _THREADINFO *pti=pgi->threads+j;
        int nsize=_gettdsize_();
        if(nsize)pti->td=hmalloc(nsize);
        pti->pg=pgi;
        pti->threads=pgi->threads;
        pti->ind=j;
        pti->sync_flag = sync_flag_init;
        pti->Nthreads=pgi->Nthreads;
        pti->MainThread=(j==0);
        pti->igrp=pgi->igrp;
        pti->indg=cbase+pti->igrp*NCorePClu +j+ib;
        pti->sib=pti->ind+1;
        pti->state = 0;
        if(pti->ind%NGG){
          pti->sib= pti->sie=0;
        } else {
          pti->sie=pti->sib+NGG-1;
          if(pti->sie>pti->Nthreads)pti->sie=pti->Nthreads;
        }
      }
    }
    bindcpu(ManageCoreId);
    int nsize=_gettdsize_();
    if(nsize)td=ti->td=hmalloc(nsize);
    gd=gi->gd;
  }else{
    md.idMainThread=1;
    int ipr=mpi_id%NProcPNode;
    if(NProcPNode<=NCluPNode){
      cbase=ipr*NCorePClu;
    }else{
      int n=NCorePClu/NThPGrp;
      if(n>1){
        int icl=ipr/n;
        int icp=ipr%n;
        cbase=icl*NCorePClu+icp*NThPGrp;
      }
    }
    ManageCoreId=cbase;
    md.Nthreads=NThPGrp;
    for(j=0;j<NThPGrp-1;j++){
      struct _THREADINFO *pti=md.threads+j;
      int nsize=_gettdsize_();
      if(nsize){
        pti->td=hmalloc(nsize);
      }
      pti->ind=j;
      pti->sync_flag = sync_flag_init;
      pti->indg=cbase+j+1;
      pti->state = 0;
      pti->Nthreads=NThPGrp;
      pti->MainThread=(j==0);
    }
    int nsize=_gettdsize_();
    if(nsize)ti->td=hmalloc(nsize);
  }
}

/* ── thread entry ── */
void threadMain(HTHREADINFO pti){
  ti=pti;td=pti->td;
  if (ThreadG){
    gi=ti->pg;gd=gi->gd;
  }
  bindthread();
  _threadmain_(pti);
}

/* ── thread lifecycle ── */
void StartThreads(TFunc tfun){
  zStartThreads(tfun,0,0,1);
}

void EndThreads(void){
  zStartThreads(NULL,0,0,1);
}

static void InitThread(HTHREADINFO pti,TFunc tfun,void*para,int detach,int nlocv){
  ClearThread(pti);
  pti->tfun=tfun;
  pti->para=para;
  pti->detach=detach;
  pti->nlocv=nlocv;
  if(tfun==NULL)return;
}

static void ClearThread(HTHREADINFO pti){
  if(pti->threadid){
    void*vres;
    pthread_cancel(pti->threadid);
    pthread_join(pti->threadid,&vres);
    if(ThreadG){
      pti->pg->sstate[pti->ind*MSBG]=0;
    }
  }
  pti->threadid=0;
  if(pti->nlocv){
    int i;
    if(pti->nlocv>=8){
      void**p=(void**)pti->locv[7];
      int nn=pti->nlocv-7;
      for(i=0;i<nn;i++){
        VFREE(p[i]);
      }
      VFREE(p);pti->locv[7]=NULL;
      pti->nlocv=7;
    }
  }
  pti->nlocv=0;
}

static void zStartThreads(TFunc tfun,void*para,int detach,int clear){
  int i,j;
  HTHREADINFO tis;
  ti=&md.tm;
  ti->threadid=0;
  ti->igrp=-1;
  ti->indg=ManageCoreId;
  ti->ind= 0;
  ti->sync_flag = sync_flag_init;
  md.state = 0;
  ti->Nthreads=md.Nthreads;
  if(tfun){
    InitThread(&md.tm,0,0,0,0);
    if(ThreadG){
      for(i=0;i<NGrpPProc;i++){
        threadGroup *pgi=md.grps[i];
        tis=pgi->threads;
        for(j=0;j<pgi->Nthreads;j++){
          InitThread(tis+j,tfun,para,detach,0);
          pthread_create(&tis[j].threadid,NULL,(TSFunc)threadMain,&tis[j]);
          if(detach)pthread_detach(tis[j].threadid);
        }
      }
    }else{
      tis=md.threads;
      md.NSpecial=0;
      for(int i=0;i<md.Nthreads-1;i++){
        md.threads[i].Nthreads=md.Nthreads;
      }
      for(int i=0;i<md.Nthreads-1;i++){
        InitThread(tis+i,tfun,para,detach,0);
        pthread_create(&tis[i].threadid,NULL,(TSFunc)threadMain,&tis[i]);
        if(detach)pthread_detach(tis[i].threadid);
      }
    }
  }else{
    if(ThreadG){
      for(i=0;i<NGrpPProc;i++){
        threadGroup *pgi=md.grps[i];
        tis=pgi->threads;
        for(j=0;j<pgi->Nthreads;j++){
          ClearThread(tis+j);
        }
      }
    }else{
      tis=md.threads;
      for(i=0;i<md.Nthreads-1;i++){
        ClearThread(tis+i);
      }
    }
  }
}
