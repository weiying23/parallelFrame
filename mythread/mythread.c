#define _GNU_SOURCE
#include <pthread.h>
#include <stdio.h>
#include "mythread.h" //
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#ifdef __APPLE__
#include <sys/types.h>
#include <sys/sysctl.h>
#else
#include <sched.h>
#endif
#ifdef DEBUG
//#define DBGSYNC
#endif
#define NGG 8

#ifndef HNEW
#define HNEW(T)     (T*)hmalloc(sizeof(T)) 
#endif
__thread THREADINFO *ti=NULL;
__thread threadGroup *gi=NULL;
__thread void*td=NULL;
__thread void*gd=NULL;
int ThreadG=-1; //分组标志
int mpi_id,NCorePClu=38,NCluPNode=16,NCorePGrp=38,NGrpPProc=4 ,NThPGrp=7 ,NProcPNode=16,ManageCoreId=-1,NThreads=3;
static void zStartThreads(TFunc tfun,void*para,int detach,int clear);
static int select_locv(int typ,short **pnlocv,void ***plocv);

// #define FUNCTION_NUM (sizeof(fun_names) / sizeof(fun_names[0]))
int sync_flag_init = 6666;  // 同步标志初始值
__attribute__((weak)) 
void thread_run() {
  // 默认空实现，Fortran 用户应提供自己的实现
  printf("Warning: using default thread_run(), please define your own.\n");
};
__attribute__((weak))
int GetVInt(int volatile *volatile p){
  return *p;
}

int getcpuid(){
#ifdef __APPLE__
  int cpu = 0;
  size_t size = sizeof(cpu);
  if (sysctlbyname("hw.logicalcpu", &cpu, &size, NULL, 0) == -1) {
    printf("warning: could not get CPU affinity , continuing...\n");
    return -1;
  }
  return cpu;
#else
  cpu_set_t mask;  //CPU核的集合
  CPU_ZERO(&mask);    //置空
  if (sched_getaffinity(0, sizeof(mask), &mask) == -1){//设置线程CPU亲和力
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

void initthreads_(int *mpi_id_,int *NCorePClu_ ,int *NCluPNode_,int *NCorePGrp_,int *NThPGrp_,int *NGrpPProc_,int *NProcPNode_,int *ManageCoreId_,int *err){
  *err=InitThreads(*mpi_id_,*NCorePClu_,* NCluPNode_,* NCorePGrp_,*NThPGrp_,*NGrpPProc_,*NProcPNode_,ManageCoreId_);
}
void startthreads_(){  
  StartThreads(thread_run);
}
void endthreads_(){  
  EndThreads();
}
void initdatasize(int tsize,int gsize){
}
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

// 弱符号定义：用户代码可以覆盖
__attribute__((weak)) int _getgdsize_(){
  return sizeof(float);
}
__attribute__((weak)) int _gettdsize_(){
  return sizeof(float);
}

int gettid(){ return ti->ind; }

int getnt(){ return ti->Nthreads; }
void _threadmain_(HTHREADINFO ti){
  bindcpu(ti->indg);
  ti->tfun();
}
#undef T
#undef F
static int select_locv(int typ,short **pnlocv,void ***plocv){
  if(typ==0){
    if(!ti) return -1;
    *pnlocv=&ti->nlocv;
    *plocv=ti->locv;
    return 0;
  }
  if(typ==1){
    if(!gi) return -1;
    *pnlocv=&gi->nlocv;
    *plocv=gi->locv;
    return 0;
  }
  if(typ==2){
    *pnlocv=&md.nlocv;
    *plocv=md.locv;
    return 0;
  }
  return -1;
}
void SetLocV(int typ,int ind,void*p){
  short *pnlocv=NULL;
  void **pl=NULL;
  int nlocv=0;
  if(select_locv(typ,&pnlocv,&pl)!=0){
    return;
  }
  if(ind<0||ind>100){
    return ;
  }
  nlocv=*pnlocv;
  if(ind<8){
    pl[ind]=p;
  }else{
    void **pp=(void**)pl[7];
    if(nlocv<=ind){
      int old_extra=(nlocv>7)?(nlocv-7):0;
      int new_nlocv=ind+4;
      int new_extra=new_nlocv-7;
      pp=hrealloc(pp,(size_t)new_extra*sizeof(*pp));
      if(!pp){
        return;
      }
      memset(pp+old_extra,0,(size_t)(new_extra-old_extra)*sizeof(*pp));
      pl[7]=pp;
      *pnlocv=ind+4;
    }
    pp[ind-7]=p;
  }
}
void *GetLocV(int typ,int ind,void*p){
  short *pnlocv=NULL;
  void **pl=NULL;
  (void)p;
  if(select_locv(typ,&pnlocv,&pl)!=0){
    return NULL;
  }
  if(ind<0||ind>100){
    return NULL;
  }
  if(ind<8) return pl[ind];
  if(*pnlocv<=ind){
    return NULL;
  }
  void **pp=(void**)pl[7];
  if(!pp){
    return NULL;
  }
  return pp[ind-7];
}
int InitThreads(int mpi_id_,int NCorePClu_ ,int NCluPNode_,int NCorePGrp_,int NThPGrp_,int NGrpPProc_,int NProcPNode_,int *ManageCoreId_){  
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
void threadMain(HTHREADINFO pti){
  ti=pti;td=pti->td;
  if (ThreadG){
    gi=ti->pg;gd=gi->gd;
  }
  bindthread();
  _threadmain_(pti);
}
void StartThreads(TFunc tfun){  
  
  zStartThreads(tfun,0,0,1);
  //printf("startthreads over\n");
}
void EndThreads(){  
 // printf("End multithread comput\n");

  zStartThreads(NULL,0,0,1);
  //printf("endthreads over\n");
}
int bindcpu(int id){
#ifdef __APPLE__
  // macOS does not support sched_setaffinity, so we just return 0
  return 0;
#else
  cpu_set_t mask;  //CPU核的集合
  CPU_ZERO(&mask);    //置空
  CPU_SET(id,&mask);   //设置亲和力值
  if (sched_setaffinity(0, sizeof(mask), &mask) == -1){//设置线程CPU亲和力
    printf("warning: could not set CPU affinity, continuing...\n");
    return -1;
  }
  return 0;
#endif
}
void bindthread(){
  if(ti) bindcpu(ti->indg);
}
void bindthread_(){
  bindthread();
}
threadProc md={0};
void ntdelay(int n){
  usleep(n/20);
  //static int vt=0; for(int i=0;i<n;i++) vt=(vt*1357+2581);
}
void ntdelay_(int n){
   //usleep(n/20);
  static int vt=0; for(int i=0;i<n;i++) vt=(vt*1357+2581);
}

#define MAXCHECK 10000000
#define VLINE ,__LINE__
#define PLINE ,int nl
#define FWAIT(NM,COP) int Wait_##NM(volatile int *pv,int cv PLINE){ \
  for(int i=0;GetVInt(pv) COP cv;i++){ \
  if(i>=MAXCHECK){ printf("wait state timeout %d %s %d\n",*pv,#COP,cv);exit(0); }\
  ntdelay(1); } \
  /* printf("Group=%d, Thread=%d: wait func runed!, %d %s %d \n", ti->igrp, ti->ind, GetVInt(pv), #COP, cv); */ \
  return 0;\
}
FWAIT(LNE,==);
FWAIT(LEQ,!=);
FWAIT(LGE,<);
FWAIT(LLE,>);

#ifdef DBGSYNC
#define THLOG
#endif
void opentf(){
#ifdef THLOG
  if(!ti->fo){
    char fn[256];
    if(ThreadG) sprintf(fn,"MTW_%2.2d_%2.2d_%2.2d.log",mpi_id,ti->igrp,ti->ind);
    else sprintf(fn,"MTW_%d.log",ti->ind);
    ti->fo=fopen(fn,"wt");
  }
#endif    
}
//SM 
void sWaitState(int state){
#ifdef DBGSYNC
  if(ti->fo) {
    fprintf(ti->fo,"TWT:3.2d %3.2d %3.2d ... ",md.nstep,state,ti->ind);
    fflush(ti->fo);
  }
#endif
  Wait_LGE(&md.state,state VLINE);
#ifdef DBGSYNC
  fprintf(ti->fo,"twt end\n");fflush(ti->fo);
#endif
}
void sWaitStater(int state){
  Wait_LLE(&md.state,state VLINE);
}
void sSetState(int state){
#ifdef DBGSYNC
  if(ti->fo) {
    fprintf(ti->fo,"TST:%3.2d %3.2d %3.2d %3.2d ...",md.nstep,ti->state,state,ti->ind);
    fflush(ti->fo);
  }
#endif
  ti->state=state;
#ifdef DBGSYNC
  fprintf(ti->fo,"tst end\n");fflush(ti->fo);
#endif
}
//MS
void mWaitSubs(int state){
#ifdef DBGSYNC
  if(ti->fo) {
    fprintf(ti->fo,"MWT:%3.2d %3.2d:",md.nstep,state);
    for(int i=0;i<md.Nthreads-1;i++){
      fprintf(ti->fo,"%3.2d ",md.threads[i].state);
    }
    fprintf(ti->fo," ...  "); fflush(ti->fo);
  }
#endif
  for(int i=0;i<md.Nthreads-1;i++){
#ifdef DBGSYNC
    fprintf(ti->fo,"V %d:%3.2d ",i,md.threads[i].state); fflush(ti->fo);
#endif
    Wait_LGE(&md.threads[i].state,state VLINE);
  }
#ifdef DBGSYNC
  fprintf(ti->fo,"mws end\n");fflush(ti->fo);
#endif
}
void mWaitSubsr(int state){
  for(int i=0;i<md.Nthreads-1;i++){
    Wait_LLE(&md.threads[i].state,state VLINE);
  }
}
void mSetSubs(int state){
#ifdef DBGSYNC
  if(ti->fo) {
    fprintf(ti->fo,"MST:%3.2d %3.2d %3.2d...",md.nstep,md.state,state);
    fflush(ti->fo);
  }
#endif
  md.state=state;
#ifdef DBGSYNC
  fprintf(ti->fo,"mss end\n");fflush(ti->fo);
#endif
}
//SG thread wait gmain state,pg->state
void sWaitGrp(int state){
#ifdef DBGSYNC
  fprintf(ti->fo,"SWG:%3.2d %3.2d %3.2d ...",md.nstep,state,ti->ind);
  fflush(ti->fo);
#endif
  Wait_LGE(&gi->state,state VLINE);
#ifdef DBGSYNC
  fprintf(ti->fo,"swg end\n");fflush(ti->fo);
#endif
}
void sWaitGrpr(int state){
  Wait_LLE(&gi->state,state VLINE);
}
//thread set state for gmain ,ti->state
#define USECSTATE
#ifdef USECSTATE
#define PSSTATE(i) &gi->sstate[i*MSBG]
#define SSSTATE(i)  gi->sstate[i*MSBG]
#else
#define PSSTATE(i) &ti->threads[i].state
#define SSSTATE(i)  ti->state
#endif
void sSetGrp(int state){
#ifdef DBGSYNC
  fprintf(ti->fo,"SSG:%3.2d %3.2d %3.2d %3.2d :",md.nstep,ti->pg->sstate[ti->ind*MSBG],state,ti->ind);
#endif
  if(ti->sib){
#ifdef DBGSYNC
    fprintf(ti->fo,"gWS S:%3.2d %3.2d:",md.nstep,state);
    for(int i= ti->sib;i<ti->sie;i++){
      fprintf(ti->fo,"%3.2d ",SSSTATE(i));
    }
    fprintf(ti->fo," ..."); fflush(ti->fo);
#endif
    for(int i= ti->sib;i<ti->sie;i++){
      Wait_LGE(PSSTATE(i),state,__LINE__);
    }
  }
  SSSTATE(ti->ind)=state;
#ifdef DBGSYNC
  fprintf(ti->fo,"ssg end\n");fflush(ti->fo);
#endif
}
//GS gmain thread wait sub threads,threads[x].state
void gWaitSubs(int state){//mt
#ifdef DBGSYNC
  fprintf(ti->fo,"GWS S:%3.2d %3.2d:",md.nstep,state);
  for(int i= ti->sib;i<ti->sie;i++){
    fprintf(ti->fo,"%3.2d ",SSSTATE(i));
  }
  fprintf(ti->fo," ... "); fflush(ti->fo);
#endif
  // printf("[GroupMain] Group=%d, pid %d, %d %d\n", ti->igrp, ti->ind, ti->sib, ti->sie);
  for(int i= ti->sib;i<ti->sie;i++){
    Wait_LGE(PSSTATE(i),state,__LINE__);
  }
#ifdef DBGSYNC
  fprintf(ti->fo,"gws s end\nGWS:%3.2d %3.2d:",md.nstep,state);
  fflush(ti->fo);
  for(int i=NGG;i<ti->Nthreads;i+=NGG){
    fprintf(ti->fo,"%3.2d ",SSSTATE(i));
  }
  fprintf(ti->fo,"... "); fflush(ti->fo);
#endif
  for(int i=NGG;i<ti->Nthreads;i+=NGG){
    Wait_LGE(PSSTATE(i),state,__LINE__);
  }
#ifdef DBGSYNC
  fprintf(ti->fo,"gws end\n");fflush(ti->fo);
#endif
}
void gWaitSubsr(int state){//mt
  for(int i= ti->sib;i<ti->sie;i++){
    Wait_LLE(PSSTATE(i),state,__LINE__);
  }
  for(int i=NGG;i<ti->Nthreads;i+=NGG){
    Wait_LLE(PSSTATE(i),state,__LINE__);
  }
}
//gmain thread set for subthread ,pg->state
void gSetSubs (int state){
#ifdef DBGSYNC
  if(ti->fo) {
    fprintf(ti->fo,"GSS:%3.2d %3.2d %3.2d %3.2d\n",md.nstep,gi->state,state,ti->ind);
    fflush(ti->fo);
  }
#endif
  gi->state=state;
}

//GM gmain thread wait mmthreads,pg->gstate
void gWaitMain(int state){//mt
  int i;
#ifdef DBGSYNC
  if(ti->fo) {
    fprintf(ti->fo,"GWM:%3.2d %3.2d %3.2d:... ",md.nstep,state,ti->igrp); fflush(ti->fo);
  }
#endif
  Wait_LGE(&gi->gstate,state VLINE);
#ifdef DBGSYNC
  if(ti->fo) {
    fprintf(ti->fo,"gwm:%3.2d %3.2d %d %3.2d end\n",md.nstep,state,ti->pg->gstate,ti->igrp); fflush(ti->fo);
  }
#endif
}
void gWaitMainr(int state){//mt
  Wait_LLE(&gi->gstate,state VLINE);
}
//gmain thread set for mmthread ,md.gstate[(x)*MSB]
void gSetMain(int state){
#ifdef DBGSYNC
  if(ti->fo) {
    fprintf(ti->fo,"GSM:%3.2d %3.2d %3.2d %3.2d ...",md.nstep,md.gstate[(ti->igrp)*MSB],state,ti->igrp);
    fflush(ti->fo);
  }
#endif
  md.gstate[(ti->igrp)*MSB]=state;
#ifdef DBGSYNC
  if(ti->fo) {
    fprintf(ti->fo,"gsm:%3.2d %3.2d %3.2d %3.2d end\n",md.nstep,md.gstate[(ti->igrp)*MSB],state,ti->igrp);
    fflush(ti->fo);
  }
#endif
}
//called by main main threads
//mmthread wait gmain thread,md.gstate[(x)*MSB]
void mWaitGrps(int state){//mmt
#ifdef DBGSYNC
  if(ti->fo) {
    fprintf(ti->fo,"MWG:%d %d %d:",md.nstep,state,md.ngrp);
    for(int i=0;i<md.ngrp;i++){ fprintf(ti->fo,"%d ",md.gstate[(i)*MSB]); }
    fprintf(ti->fo,"\n ");
    fflush(ti->fo);
  }
#endif
  for(int i=0;i<md.ngrp;i++){
    Wait_LGE(&md.gstate[(i)*MSB],state VLINE);
  }
}
void mWaitGrpsr(int state){//mmt
  for(int i=0;i<md.ngrp;i++){
    Wait_LLE(&md.gstate[(i)*MSB],state VLINE);
  }
}
//mmthread set for gmain threads,md.grps[x].gstate
void mSetGrps(int state){// mmt
#ifdef DBGSYNC
  if(ti->fo) {
    fprintf(ti->fo,"MSG:%3.2d %3.2d %3.2d ...",md.nstep,md.state,state);
    fflush(ti->fo);
  }
#endif
  md.state=state;
  for(int i=0;i<md.ngrp;i++){
    md.grps[i]->gstate=state;
  }
#ifdef DBGSYNC
  if(ti->fo) {
    fprintf(ti->fo,"MSG:%3.2d %3.2d ...",md.nstep,state);
    for(int i=0;i<md.ngrp;i++){
      fprintf(ti->fo,"%d:%3.2d ",i,md.grps[i]->gstate); 
    }
    fprintf(ti->fo,"msg end\n");fflush(ti->fo);
  }
#endif
}

// Fortran interface
//SM
void swaitstate_ (int *state){ sWaitState (*state); }
void swaitstater_(int *state){ sWaitStater(*state); }
void ssetstate_  (int *state){ sSetState  (*state);   }
//MS
void mwaitsubs_  (int *state)    { mWaitSubs (*state);}
void mwaitsubsr_ (int *state)    { mWaitSubsr(*state);}
void msetsubs_  (int *state)    { mSetSubs (*state);}

//GS
void gwaitsubs_  (int *state){ gWaitSubs  (*state); }
void gwaitsubsr_ (int *state){ gWaitSubsr (*state); }
void gsetsubs_   (int *state){ gSetSubs   (*state); }
//GM
void gwaitmain_  (int *state){ gWaitMain  (*state); }
void gwaitmainr_ (int *state){ gWaitMainr (*state); }
void gsetmain_   (int *state){ gSetMain   (*state); }
//MG
void mwaitgrps_  (int *state){  mWaitGrps (*state); }
void mwaitgrpsr_ (int *state){  mWaitGrpsr(*state); }
void msetgrps_   (int *state){  mSetGrps  (*state); }

static inline uint64_t atsc(void) {
#ifdef __ARM_ARCH
  uint64_t tsc;
  asm volatile("mrs %0, cntvct_el0" : "=r" (tsc));    //读取系统时间戳
  return tsc&0xFFFFFFFFFFFFFFL;//most tsc is 56bits
#else
  uint64_t a,d;
  asm volatile("rdtsc " : "=a" (a),"=d"(d));    //读取系统时间戳
  return (d<<32)|a;
#endif
}
void tscinit(){
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
void prtsc_(){
  prtsc("MT");
}
#define VFREE(p) if(p)free(p)
static void ClearThread(HTHREADINFO pti){
  if(pti->threadid){
    //printf("mpi_id:%d threadid:%d indg:%d end stage1\n",mpi_id,pti->threadid,pti->indg);
    void*vres;
    //printf("mpi_id:%d threadid:%d indg:%d end stage12\n",mpi_id,pti->threadid,pti->indg);

    pthread_cancel(pti->threadid);
    //printf("mpi_id:%d threadid:%d indg:%d end stage13\n",mpi_id,pti->threadid,pti->indg);

    pthread_join(pti->threadid,&vres);
    //printf("mpi_id:%d threadid:%d indg:%d end stage14\n",mpi_id,pti->threadid,pti->indg);
    
    //printf("mpi_id:%d threadid:%d ind:%d end stage16\n",mpi_id,pti->threadid,pti->ind);
    if(ThreadG){
      pti->pg->sstate[pti->ind*MSBG]=0;
    }
    //printf("mpi_id:%d threadid:%d indg:%d end stage15\n",mpi_id,pti->threadid,pti->indg);

  }
  //printf("mpi_id:%d threadid:%d indg:%d end stage2\n",mpi_id,pti->threadid,pti->indg);

  pti->threadid=0;
  //printf("mpi_id:%d indg:%d end stage2\n",mpi_id,pti->threadid);

  if(pti->nlocv){
    //printf("mpi_id:%d threadid:%d indg:%d end stage3\n",mpi_id,pti->threadid,pti->indg);

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
  //printf("mpi_id:%d threadid:%d indg:%d end stage4\n",mpi_id,pti->threadid,pti->indg);

  pti->nlocv=0;
}


void initmd(){
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
      mg=-1;//fixme: only for this machine
      mt=-1;
      ManageCoreId = (ipr+NCorePClu*NCluPNode)
          ;
    }else{
      if(ManageCoreId<NCorePClu){
        mg=0;mt=ManageCoreId;
        ManageCoreId+=cbase;
      }else{
        mg=-1;mt=-1;//fixme: only for this machine
        ManageCoreId = (ipr+NCorePClu*NCluPNode);
      }
    }
    for(i=0;i<NGrpPProc;i++){
      bindcpu(cbase+i*NCorePClu);// bind to core in grp for hbm
      threadGroup *pgi=md.grps[i]=HNEW(threadGroup);
      memset(pgi,0,sizeof(threadGroup));
      pgi->pmd=&md;
      pgi->ngrp=NGrpPProc;
      pgi->Nthreads=NThPGrp;
      //printf("Nthread1:%d\n",pgi->Nthreads);
      pgi->idMainThread=0;
      int ib=0;
      if(mg==i || mt==-1){
        if(!mt){ ib++; pgi->Nthreads--; }
        gi=ti->pg=pgi;
      }
      //printf("Nthread2:%d\n",pgi->Nthreads);

      pgi->igrp=i;
      pgi->NSpecial=0;
      int nsize=_getgdsize_();
      if(nsize)pgi->gd=hmalloc(nsize);
      for(j=0;j<pgi->Nthreads;j++){
        //printf("Nthread3:%d\n",pgi->Nthreads);

        struct _THREADINFO *pti=pgi->threads+j;
        int nsize=_gettdsize_();
        if(nsize)pti->td=hmalloc(nsize);
        pti->pg=pgi;
        pti->threads=pgi->threads;
        //printf("pti->ind:%d \n",j);
        pti->ind=j;
        pti->sync_flag = sync_flag_init;
        pti->Nthreads=pgi->Nthreads;
        //printf("Nthread1:%d\n",pgi->Nthreads);

        pti->MainThread=(j==0);//pgi->idMainThread);
        pti->igrp=pgi->igrp;
        pti->indg=cbase+pti->igrp*NCorePClu +j+ib;
        //printf("pti->indg:%d\n",pti->indg);
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
    //printf("mpi_id:%d managecoreid:%d\n",mpi_id,ManageCoreId);
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
    //md.mthreads=HNEW(threadGroup)
    //printf("mpi_id:%d managecoreid:%d\n",mpi_id,cbase);
    md.Nthreads=NThPGrp;
    for(j=0;j<NThPGrp-1;j++){
      struct _THREADINFO *pti=md.threads+j;
      int nsize=_gettdsize_();
      if(nsize){
        //printf("nsize:%d\n",nsize);
        pti->td=hmalloc(nsize);
      }
      pti->ind=j;
      pti->sync_flag = sync_flag_init;
      
      pti->indg=cbase+j+1;
      pti->state = 0;
      //printf("kkkkk:%d\n",pti->indg);
      pti->Nthreads=NThPGrp;
      pti->MainThread=(j==0);
    }
    int nsize=_gettdsize_();
    if(nsize)ti->td=hmalloc(nsize);
  }
}
static void InitThread(HTHREADINFO pti,TFunc tfun,void*para,int detach,int nlocv){
  ClearThread(pti);
  pti->tfun=tfun;
  pti->para=para;
  pti->detach=detach;
  pti->nlocv=nlocv;
  if(tfun==NULL)return;
}
void zStartThreads(TFunc tfun,void*para,int detach,int clear){
  int ib,ie;
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
  //initthreads(mpi_id,38,16,38,8,4,4,2); 
  //InitThreads(int mpi_id_,int NCorePClu_ ,int NCluPNode_,int NCorePGrp_,int NThPGrp_,int NGrpPProc_,int NProcPNode_,int *ManageCoreId_)
  if(tfun){
    InitThread(&md.tm,0,0,0,0);
    //printf("ti->threadid %d\n",ti->indg);
    if(ThreadG){

      for(i=0;i<NGrpPProc;i++){
        threadGroup *pgi=md.grps[i];
        tis=pgi->threads;
        for(j=0;j<pgi->Nthreads;j++){
          //printf("grp:%d pgi->indg:%d\n",i,pgi->threads[j].indg);
         
          InitThread(tis+j,tfun,para,detach,0);
          //printf("tis->indg:%d\n",(tis+i)->indg);
          pthread_create(&tis[j].threadid,NULL,(TSFunc)threadMain,&tis[j]);
          if(detach)pthread_detach(tis[j].threadid);
        }
      }
    }else{
      tis=md.threads;
      //printf("stage3\n");

      md.NSpecial=0;

      for(int i=0;i<md.Nthreads-1;i++){
        md.threads[i].Nthreads=md.Nthreads;
      }
      for(int i=0;i<md.Nthreads-1;i++){
        InitThread(tis+i,tfun,para,detach,0);
        pthread_create(&tis[i].threadid,NULL,(TSFunc)threadMain,&tis[i]);
        if(detach)pthread_detach(tis[i].threadid);
      }
      //for(int i=0;i<md.Nthreads-1;i++){
        //InitThread(tis+i,tfun,para,detach,0);
        //pthread_join(tis[i].threadid,NULL);

        //if(detach)pthread_detach(tis[i].threadid);
      //}

    }    
  }else{
    //printf("1\n");
    //printf("end threads mpi_id:%d ti->indg:%d\n",mpi_id,ti->indg);
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

void padr_(int *id,double*p){
  printf("ZZZZZZZZZZZZZZz : %d %p\n",*id,p);
}

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

static inline int mt_taskpool_typ(){
  if(ThreadG && gi && ti && ti->igrp>=0) return 1;
  return 2;
}

static inline int mt_taskpool_nworkers(){
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
