#include "mythread_sync.h"
#include "mythread_util.h"
#include <stdlib.h>

#ifdef DEBUG
/* #define DBGSYNC */
#endif

__attribute__((weak))
int GetVInt(int volatile *volatile p){
  return *p;
}

#define MAXCHECK 10000000
#define VLINE ,__LINE__
#define PLINE ,int nl
#define FWAIT(NM,COP) int Wait_##NM(volatile int *pv,int cv PLINE){ \
  for(int i=0;GetVInt(pv) COP cv;i++){ \
  if(i>=MAXCHECK){ printf("wait state timeout %d %s %d\n",*pv,#COP,cv);exit(0); }\
  ntdelay(1); } \
  return 0;\
}
FWAIT(LNE,==);
FWAIT(LEQ,!=);
FWAIT(LGE,<);
FWAIT(LLE,>);

/* SM */
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

/* MS */
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

/* SG */
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

#define NGG 8
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

/* GS */
void gWaitSubs(int state){
#ifdef DBGSYNC
  fprintf(ti->fo,"GWS S:%3.2d %3.2d:",md.nstep,state);
  for(int i= ti->sib;i<ti->sie;i++){
    fprintf(ti->fo,"%3.2d ",SSSTATE(i));
  }
  fprintf(ti->fo," ... "); fflush(ti->fo);
#endif
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
void gWaitSubsr(int state){
  for(int i= ti->sib;i<ti->sie;i++){
    Wait_LLE(PSSTATE(i),state,__LINE__);
  }
  for(int i=NGG;i<ti->Nthreads;i+=NGG){
    Wait_LLE(PSSTATE(i),state,__LINE__);
  }
}
void gSetSubs(int state){
#ifdef DBGSYNC
  if(ti->fo) {
    fprintf(ti->fo,"GSS:%3.2d %3.2d %3.2d %3.2d\n",md.nstep,gi->state,state,ti->ind);
    fflush(ti->fo);
  }
#endif
  gi->state=state;
}

/* GM */
void gWaitMain(int state){
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
void gWaitMainr(int state){
  Wait_LLE(&gi->gstate,state VLINE);
}
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

/* MG */
void mWaitGrps(int state){
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
void mWaitGrpsr(int state){
  for(int i=0;i<md.ngrp;i++){
    Wait_LLE(&md.gstate[(i)*MSB],state VLINE);
  }
}
void mSetGrps(int state){
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
