#include "mythread_util.h"
#include <unistd.h>

void ntdelay(int n){
  usleep(n/20);
}

void ntdelay_(int n){
  static int vt=0; for(int i=0;i<n;i++) vt=(vt*1357+2581);
}

#ifdef THLOG
void opentf(void){
  if(!ti->fo){
    char fn[256];
    if(ThreadG) sprintf(fn,"MTW_%2.2d_%2.2d_%2.2d.log",mpi_id,ti->igrp,ti->ind);
    else sprintf(fn,"MTW_%d.log",ti->ind);
    ti->fo=fopen(fn,"wt");
  }
}
#else
void opentf(void){}
#endif
