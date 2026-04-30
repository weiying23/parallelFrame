#ifndef MTHREAD_THREAD_H_INCLUDED
#define MTHREAD_THREAD_H_INCLUDED
#include "mythread_types.h"

__BEGIN_DECLS
void initmd(void);
int getcpuid(void);
void thread_run(void);
int gettid(void);
int getnt(void);
void setthread(int NCorePClu_,int NThPGrp_,int NGrpPProc_,int NProcPNode_,int ManageCoreId_);
int InitThreads(int mpi_id_,int NCorePClu_,int NCluPNode_,int NCorePGrp_,int NThPGrp_,int NGrpPProc_,int NProcPNode_,int *ManageCoreId_);
void StartThreads(TFunc tfun);
void EndThreads(void);
void bindthread(void);
int bindcpu(int id);
void _threadmain_(HTHREADINFO ti);
__END_DECLS

#endif
