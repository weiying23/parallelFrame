#ifndef MTHREAD_TYPES_H_INCLUDED
#define MTHREAD_TYPES_H_INCLUDED
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
#define __BEGIN_DECLS extern "C" {
#define __END_DECLS }
#else
#define __BEGIN_DECLS
#define __END_DECLS
#endif

#ifndef RFB
#define RFB 1
#endif
#ifndef MCOREPC
#define MCOREPC 64
#endif
#ifndef MCLUST
#define MCLUST  20
#endif
#ifndef MSB
#define MSB     16
#endif
#ifndef MSBG
#define MSBG    16
#endif

typedef struct _THREADINFO THREADINFO,*HTHREADINFO;
typedef void (*TFunc)(void);
typedef void (*Func)();
typedef void* (*TSFunc)(void*);

typedef struct _THREADINFO{
  short Nthreads,MainThread,igrp;
  short nlocv,mlocv,detach,indg,sib,sie;
  int ind;
  int sync_flag;
  int nstep;
  int volatile state;
  uint64_t tscbs[16];
  uint64_t tscds[16];
  struct _threadGroup *pg;
  struct _THREADINFO  *threads;
  TFunc tfun;
  void *para;
  pthread_t threadid;
  FILE *fo;
  void*td;
  void*locv[8];
}THREADINFO;

typedef struct _threadGroup{
  short igrp,ngrp,mainGroup,Nthreads,idMainThread,NSpecial;
  short nlocv,mlocv,PThreadInited;
  int nstep;
  int volatile state;
  int volatile gstate;
  struct _threadProc *pmd;
  int volatile sstate[(MCOREPC)*MSBG];
  THREADINFO threads[MCOREPC];
  void*gd;
  void*locv[8];
}threadGroup;

typedef struct _threadProc{
  short Nthreads,NSpecial,idMainThread,ngrp;
  short nlocv,mlocv,PThreadInited;
  int nstep;
  int mpi_id,mpi_npe;
  THREADINFO tm;
  THREADINFO threads[MCOREPC];
  threadGroup *grps[MCLUST];
  int volatile state;
  int volatile gstate[(MCLUST)*MSB];
  void*locv[8];
}threadProc;

typedef void (*mt_task_fn)(void *ctx);
typedef struct {
  mt_task_fn fn;
  void *ctx;
} mt_task;
typedef struct _mt_taskpool mt_taskpool;

#define MT_TASKPOOL_BLOCK 1
#define MT_TASKPOOL_SPIN  2
#define MT_TASKPOOL_TRY   4

__BEGIN_DECLS
extern threadProc md;
extern int mpi_id, NCorePClu,NCluPNode, NCorePGrp,NGrpPProc,NThPGrp,NProcPNode,ManageCoreId,NThreads;
extern __thread THREADINFO *ti;
extern __thread threadGroup *gi;
extern int ThreadG;
__END_DECLS

#define PTI HTHREADINFO ti
#define PTIM HTHREADINFO ti,
#define VTI ti
#define VTIM ti,

#ifndef USEHBM
#define hmalloc malloc
#define hrealloc realloc
#else
__BEGIN_DECLS
extern void*hmalloc(size_t size);
extern void*hrealloc(void*p,size_t size);
__END_DECLS
#endif

#endif
