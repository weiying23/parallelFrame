#ifndef MTHREAD_FORTRAN_H_INCLUDED
#define MTHREAD_FORTRAN_H_INCLUDED
#include "mythread_types.h"

__BEGIN_DECLS
/* 初始化 */
void inithreads_(int *mpi_id_,int *NCorePClu_,int *NCluPNode_,int *NCorePGrp_,int *NThPGrp_,int *NGrpPProc_,int *NProcPNode_,int *ManageCoreId_,int *err);
void startthreads_(void);
void endthreads_(void);

/* CPU 绑定 */
void bindthread_(void);

/* 同步 — SM */
void swaitstate_(int *state);
void swaitstater_(int *state);
void ssetstate_(int *state);
/* 同步 — MS */
void mwaitsubs_(int *state);
void mwaitsubsr_(int *state);
void msetsubs_(int *state);
/* 同步 — SG */
void swaitgrp_(int *state);
void swaitgrpr_(int *state);
void ssetgrp_(int *state);
/* 同步 — GS */
void gwaitsubs_(int *state);
void gwaitsubsr_(int *state);
void gsetsubs_(int *state);
/* 同步 — GM */
void gwaitmain_(int *state);
void gwaitmainr_(int *state);
void gsetmain_(int *state);
/* 同步 — MG */
void mwaitgrps_(int *state);
void mwaitgrpsr_(int *state);
void msetgrps_(int *state);

/* 计时器 */
void tscb_(int *id_);
void tsce_(int *id_);
void tsceb_(int *id_);
void prtsc_(void);

/* 工具 */
void ntdelay_(int n);
__END_DECLS

#endif
