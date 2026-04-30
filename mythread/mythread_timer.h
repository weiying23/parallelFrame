#ifndef MTHREAD_TIMER_H_INCLUDED
#define MTHREAD_TIMER_H_INCLUDED
#include "mythread_types.h"

__BEGIN_DECLS
void tscinit(void);
void tscb(int id);
void tsce(int id);
void tsceb(int id);
void prtsc(const char*tag);
void tscb_(int *id_);
void tsce_(int *id_);
void tsceb_(int *id_);
void prtsc_(void);
__END_DECLS

#endif
