#include "mythread_fortran.h"
#include "mythread_thread.h"
#include "mythread_sync.h"
#include "mythread_timer.h"

/* ── 初始化 ── */
void inithreads_(int *mpi_id_,int *NCorePClu_,int *NCluPNode_,int *NCorePGrp_,int *NThPGrp_,int *NGrpPProc_,int *NProcPNode_,int *ManageCoreId_,int *err){
  *err=InitThreads(*mpi_id_,*NCorePClu_,* NCluPNode_,* NCorePGrp_,*NThPGrp_,*NGrpPProc_,*NProcPNode_,ManageCoreId_);
}

void startthreads_(void){
  StartThreads(thread_run);
}

void endthreads_(void){
  EndThreads();
}

/* ── CPU 绑定 ── */
void bindthread_(void){
  bindthread();
}

/* ── 同步 — SM ── */
void swaitstate_(int *state){ sWaitState(*state); }
void swaitstater_(int *state){ sWaitStater(*state); }
void ssetstate_(int *state){ sSetState(*state); }

/* ── 同步 — MS ── */
void mwaitsubs_(int *state){ mWaitSubs(*state); }
void mwaitsubsr_(int *state){ mWaitSubsr(*state); }
void msetsubs_(int *state){ mSetSubs(*state); }

/* ── 同步 — GS ── */
void gwaitsubs_(int *state){ gWaitSubs(*state); }
void gwaitsubsr_(int *state){ gWaitSubsr(*state); }
void gsetsubs_(int *state){ gSetSubs(*state); }

/* ── 同步 — GM ── */
void gwaitmain_(int *state){ gWaitMain(*state); }
void gwaitmainr_(int *state){ gWaitMainr(*state); }
void gsetmain_(int *state){ gSetMain(*state); }

/* ── 同步 — MG ── */
void mwaitgrps_(int *state){ mWaitGrps(*state); }
void mwaitgrpsr_(int *state){ mWaitGrpsr(*state); }
void msetgrps_(int *state){ mSetGrps(*state); }
