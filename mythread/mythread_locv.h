#ifndef MTHREAD_LOCV_H_INCLUDED
#define MTHREAD_LOCV_H_INCLUDED
#include "mythread_types.h"

__BEGIN_DECLS
void SetLocV(int typ,int ind,void *p);
void *GetLocV(int typ,int ind,void *p);
__END_DECLS

#endif
