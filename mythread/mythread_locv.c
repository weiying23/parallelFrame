#include "mythread_locv.h"
#include <string.h>

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
