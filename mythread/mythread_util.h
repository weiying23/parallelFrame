#ifndef MTHREAD_UTIL_H_INCLUDED
#define MTHREAD_UTIL_H_INCLUDED
#include "mythread_types.h"

__BEGIN_DECLS
void ntdelay(int n);
void ntdelay_(int n);
void opentf(void);
__END_DECLS

#ifndef HNEW
#define HNEW(T) ((T*)hmalloc(sizeof(T)))
#endif

#endif
