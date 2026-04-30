#ifndef MTHREAD_POOL_H_INCLUDED
#define MTHREAD_POOL_H_INCLUDED
#include "mythread_types.h"

__BEGIN_DECLS
int mt_taskpool_attach(int slot,int capacity,int flags);
int mt_taskpool_detach(int slot);
int mt_taskpool_begin(int slot);
int mt_taskpool_submit(int slot,mt_task_fn fn,void *ctx);
int mt_taskpool_close(int slot);
int mt_taskpool_wait(int slot);
int mt_taskpool_shutdown(int slot);
int mt_taskpool_worker_loop(int slot);
__END_DECLS

#endif
