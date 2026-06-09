#include <math.h>
#include <mpi.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mythread/mythread.h"

static int    cfg_NX=14000,cfg_NY=14000,cfg_NT=480;
static double cfg_A=10.0,cfg_DT=0.001,cfg_C0=0.1;
static int    cfg_USE_FIXED_DOMAIN=0;
static double cfg_DX,cfg_DY,cfg_LX,cfg_LY,cfg_DT2,cfg_CFL_X,cfg_CFL_Y,cfg_CFL_SUM2;
static int    cfg_HALO=1,cfg_N_GROUPS=2,cfg_N_WORKERS=2,cfg_GROUP_DECOMP=0;
static int    cfg_NCorePClu=5,cfg_NCluPNode=2,cfg_NCorePGrp=4,cfg_ManageCoreId=4;
static int    cfg_CoreOffset=0,cfg_ClustOffset=0;

static void cfg_compute_derived(void){
  if(cfg_USE_FIXED_DOMAIN){cfg_LX=1.0;cfg_LY=1.0;cfg_DX=cfg_LX/(cfg_NX-1);cfg_DY=cfg_LY/(cfg_NY-1);}
  else{cfg_DX=0.01;cfg_DY=0.01;cfg_LX=(cfg_NX-1)*cfg_DX;cfg_LY=(cfg_NY-1)*cfg_DY;}
  cfg_DT2=cfg_DT*cfg_DT;cfg_CFL_X=cfg_C0*cfg_DT/cfg_DX;cfg_CFL_Y=cfg_C0*cfg_DT/cfg_DY;
  cfg_CFL_SUM2=cfg_CFL_X*cfg_CFL_X+cfg_CFL_Y*cfg_CFL_Y;
}

#define NX cfg_NX
#define NY cfg_NY
#define NT cfg_NT
#define AA cfg_A
#define DT cfg_DT
#define C0 cfg_C0
#define DX cfg_DX
#define DY cfg_DY
#define DT2 cfg_DT2
#define HALO cfg_HALO
#define CFL_SUM2 cfg_CFL_SUM2

typedef struct {
  int mx,my,nx,ny,mr,ms,px,py,px_id,py_id;
  int nb_l,nb_r,nb_d,nb_u;
} Dom;

static Dom g_d={0};
static GroupField *g_gfields=NULL;
static const mythread_decomp *g_dc=NULL;
static mythread_mpi_ctx g_mpi_ctx={0};

typedef struct {
  int gid,tid;
  int y_begin,y_end,x_begin,x_end;
  GroupField *gf;
  double partial_energy;
} ThreadTask;

int _gettdsize_(void){return (int)sizeof(ThreadTask);}
int _getgdsize_(void){return 0;}

static inline int gx(int gid,int lx){return g_d.mx+(g_dc->group_tiles[gid].x_begin+lx-HALO);}
static inline int gy(int gid,int ly){return g_d.my+(g_dc->group_tiles[gid].y_begin+ly-HALO);}

static int init_state(void){return 1;}
static int copy_state(void){return 2;}
static int compute_state(int step){return 2*step+3;}
static int boundary_state(int step){return 2*step+4;}
static int energy_state(void){return 2*NT+5;}

static double init_val(int gy_,int gx_){
  double cx=0.5*(cfg_USE_FIXED_DOMAIN?1.0:((cfg_NX-1)*cfg_DX)),
         cy=0.5*(cfg_USE_FIXED_DOMAIN?1.0:((cfg_NY-1)*cfg_DY)),
         sigma=0.06*((cfg_USE_FIXED_DOMAIN?1.0:(cfg_DX*(cfg_NX-1)))<(cfg_USE_FIXED_DOMAIN?1.0:(cfg_DY*(cfg_NY-1)))?
           (cfg_USE_FIXED_DOMAIN?1.0:(cfg_DX*(cfg_NX-1))):(cfg_USE_FIXED_DOMAIN?1.0:(cfg_DY*(cfg_NY-1))));
  return AA*exp(-((gx_*DX-cx)*(gx_*DX-cx)+(gy_*DY-cy)*(gy_*DY-cy))/(2.*sigma*sigma));
}

static void gf_apply_d(GroupField*gf,int gid){
  const mythread_tile*t=&g_dc->group_tiles[gid];
  int s=gf->stride,h=gf->ny_padded,nx=t->nx,ny=t->ny;
  int global_x0=g_d.mx+t->x_begin,global_x1=g_d.mx+t->x_end;
  int global_y0=g_d.my+t->y_begin,global_y1=g_d.my+t->y_end;

  if(global_x0==0)for(int y=0;y<h;y++)gf->u_curr[GFIDX(gf,y,HALO)]=0.0;
  if(global_x1==NX){int xr=HALO+nx-1;for(int y=0;y<h;y++)gf->u_curr[GFIDX(gf,y,xr)]=0.0;}

  if(g_d.nb_l<0&&mythread_decomp_is_domain_boundary(g_dc,gid,MYTHREAD_NEIGHBOR_LEFT))
    for(int y=0;y<h;y++)gf->u_curr[GFIDX(gf,y,0)]=0.0;
  if(g_d.nb_r<0&&mythread_decomp_is_domain_boundary(g_dc,gid,MYTHREAD_NEIGHBOR_RIGHT)){
    int xh=HALO+nx;for(int y=0;y<h;y++)gf->u_curr[GFIDX(gf,y,xh)]=0.0;
  }

  if(g_d.nb_d<0&&mythread_decomp_is_domain_boundary(g_dc,gid,MYTHREAD_NEIGHBOR_DOWN)){
    memset(&gf->u_curr[0],0,(size_t)s*sizeof(double));
    if(global_y0==0)memset(&gf->u_curr[GFIDX(gf,HALO,0)],0,(size_t)s*sizeof(double));
  }
  if(g_d.nb_u<0&&mythread_decomp_is_domain_boundary(g_dc,gid,MYTHREAD_NEIGHBOR_UP)){
    memset(&gf->u_curr[(ny+HALO)*s],0,(size_t)s*sizeof(double));
    if(global_y1==NY)memset(&gf->u_curr[GFIDX(gf,ny,0)],0,(size_t)s*sizeof(double));
  }
}

static void apply_d_all(void){
  for(int g=0;g<g_dc->n_groups;g++)gf_apply_d(&g_gfields[g],g);
}

static void halo_all(void){
  mythread_halo_exchange_intra(g_gfields,g_dc);
  mythread_halo_exchange_mpi(g_gfields,g_dc,&g_mpi_ctx,(uintptr_t)MPI_COMM_WORLD);
}

static void copy_prev_block(GroupField*gf,int gid,int y_begin,int y_end){
  const mythread_tile*t=&g_dc->group_tiles[gid];
  if(y_begin==HALO)
    memcpy(&gf->u_prev[GFIDX(gf,0,0)],&gf->u_curr[GFIDX(gf,0,0)],
           (size_t)gf->stride*sizeof(double));
  for(int y=y_begin;y<y_end;y++)
    memcpy(&gf->u_prev[GFIDX(gf,y,0)],&gf->u_curr[GFIDX(gf,y,0)],
           (size_t)gf->stride*sizeof(double));
  if(y_end==HALO+t->ny)
    memcpy(&gf->u_prev[GFIDX(gf,HALO+t->ny,0)],&gf->u_curr[GFIDX(gf,HALO+t->ny,0)],
           (size_t)gf->stride*sizeof(double));
}

static void gf_swap_all(void){
  for(int g=0;g<g_dc->n_groups;g++)group_field_swap(&g_gfields[g]);
}

static void init_field_block(GroupField*gf,int gid,int y_begin,int y_end,int x_begin,int x_end){
  for(int y=y_begin;y<y_end;y++){int gy_=gy(gid,y);
    for(int x=x_begin;x<x_end;x++){int gx_=gx(gid,x);
      double v=(gx_!=0&&gx_!=NX-1&&gy_!=0&&gy_!=NY-1)?init_val(gy_,gx_):0.0;
      gf->u_curr[GFIDX(gf,y,x)]=v;gf->u_prev[GFIDX(gf,y,x)]=v;gf->u_next[GFIDX(gf,y,x)]=0.0;}}
}

static void comp_region(GroupField*gf,int gid,int y_begin,int y_end,int x_begin,int x_end){
  const double inv_dx2=1.0/(DX*DX),inv_dy2=1.0/(DY*DY),c2dt2=C0*C0*DT2;
  for(int y=y_begin;y<y_end;y++){int gy_=gy(gid,y);if(gy_==0||gy_==NY-1)continue;
    for(int x=x_begin;x<x_end;x++){int gx_=gx(gid,x);if(gx_==0||gx_==NX-1)continue;
      double u=gf->u_curr[GFIDX(gf,y,x)];
      double d2x=(gf->u_curr[GFIDX(gf,y,x-1)]-2*u+gf->u_curr[GFIDX(gf,y,x+1)])*inv_dx2;
      double d2y=(gf->u_curr[GFIDX(gf,y-1,x)]-2*u+gf->u_curr[GFIDX(gf,y+1,x)])*inv_dy2;
      gf->u_next[GFIDX(gf,y,x)]=2*u-gf->u_prev[GFIDX(gf,y,x)]+c2dt2*(d2x+d2y);}}
}

static void comp_interior_block(GroupField*gf,int gid,int y_begin,int y_end,int x_begin,int x_end){
  const mythread_tile*t=&g_dc->group_tiles[gid];int ny=t->ny,nx=t->nx;
  if(g_d.nb_d>=0&&mythread_decomp_is_domain_boundary(g_dc,gid,MYTHREAD_NEIGHBOR_DOWN)&&y_begin<HALO+1)y_begin=HALO+1;
  if(g_d.nb_u>=0&&mythread_decomp_is_domain_boundary(g_dc,gid,MYTHREAD_NEIGHBOR_UP)&&y_end>ny)y_end=ny;
  if(g_d.nb_l>=0&&mythread_decomp_is_domain_boundary(g_dc,gid,MYTHREAD_NEIGHBOR_LEFT)&&x_begin<HALO+1)x_begin=HALO+1;
  if(g_d.nb_r>=0&&mythread_decomp_is_domain_boundary(g_dc,gid,MYTHREAD_NEIGHBOR_RIGHT)&&x_end>nx)x_end=nx;
  if(y_begin<y_end&&x_begin<x_end)comp_region(gf,gid,y_begin,y_end,x_begin,x_end);
}

static void comp_boundary_block(GroupField*gf,int gid,int y_begin,int y_end,int x_begin,int x_end){
  const mythread_tile*t=&g_dc->group_tiles[gid];int ny=t->ny,nx=t->nx;
  int lr=HALO,ur=ny,lc=HALO,rc=nx;
  if(x_begin<HALO)x_begin=HALO;if(x_end>HALO+nx)x_end=HALO+nx;
  if(y_begin<HALO)y_begin=HALO;if(y_end>HALO+ny)y_end=HALO+ny;
  if(g_d.nb_d>=0&&mythread_decomp_is_domain_boundary(g_dc,gid,MYTHREAD_NEIGHBOR_DOWN)&&y_begin<=lr&&lr<y_end)
    comp_region(gf,gid,lr,lr+1,x_begin,x_end);
  if(g_d.nb_u>=0&&mythread_decomp_is_domain_boundary(g_dc,gid,MYTHREAD_NEIGHBOR_UP)&&y_begin<=ur&&ur<y_end)
    comp_region(gf,gid,ur,ur+1,x_begin,x_end);
  if(g_d.nb_l>=0&&mythread_decomp_is_domain_boundary(g_dc,gid,MYTHREAD_NEIGHBOR_LEFT)&&x_begin<=lc&&lc<x_end)
    comp_region(gf,gid,y_begin,y_end,lc,lc+1);
  if(g_d.nb_r>=0&&mythread_decomp_is_domain_boundary(g_dc,gid,MYTHREAD_NEIGHBOR_RIGHT)&&x_begin<=rc&&rc<x_end)
    comp_region(gf,gid,y_begin,y_end,rc,rc+1);
}

static double energy_block(GroupField*gf,int gid,int y_begin,int y_end,int x_begin,int x_end){
  double my_e=0.0;
  for(int y=y_begin;y<y_end;y++){int gy_=gy(gid,y);if(gy_==0||gy_==NY-1)continue;
    for(int x=x_begin;x<x_end;x++){int gx_=gx(gid,x);if(gx_==0||gx_==NX-1)continue;
      double u=gf->u_curr[GFIDX(gf,y,x)],up=gf->u_prev[GFIDX(gf,y,x)];
      double ut=(u-up)/DT,dxc=(gf->u_curr[GFIDX(gf,y,x+1)]-gf->u_curr[GFIDX(gf,y,x-1)])/(2*DX);
      double dyc=(gf->u_curr[GFIDX(gf,y+1,x)]-gf->u_curr[GFIDX(gf,y-1,x)])/(2*DY);
      double dxp=(gf->u_prev[GFIDX(gf,y,x+1)]-gf->u_prev[GFIDX(gf,y,x-1)])/(2*DX);
      double dyp=(gf->u_prev[GFIDX(gf,y+1,x)]-gf->u_prev[GFIDX(gf,y-1,x)])/(2*DY);
      my_e+=0.5*(ut*ut+C0*C0*(dxc*dxp+dyc*dyp))*DX*DY;
    }}
  return my_e;
}

static void setup_domain(int mr,int ms){
  mythread_decomp*tmp=mythread_decomp_create(NX,NY,HALO,ms,0,MYTHREAD_DECOMP_XY_2D);
  int px=tmp->gx,py=tmp->gy;mythread_decomp_free(tmp);
  g_d.px=px;g_d.py=py;g_d.px_id=mr%px;g_d.py_id=mr/px;
  {int tn=NX/px,rn=NX%px,tm=NY/py,rm=NY%py;
   int ox=g_d.px_id*tn+(g_d.px_id<rn?g_d.px_id:rn),oy=g_d.py_id*tm+(g_d.py_id<rm?g_d.py_id:rm);
   g_d.mx=ox;g_d.nx=tn+(g_d.px_id<rn?1:0);
   g_d.my=oy;g_d.ny=tm+(g_d.py_id<rm?1:0);}
  g_d.nb_l=(g_d.px_id>0)?(mr-1):-1;g_d.nb_r=(g_d.px_id+1<px)?(mr+1):-1;
  g_d.nb_d=(g_d.py_id>0)?(mr-px):-1;g_d.nb_u=(g_d.py_id+1<py)?(mr+px):-1;
  g_d.mr=mr;g_d.ms=ms;
  g_mpi_ctx.mpirank_left=g_d.nb_l;g_mpi_ctx.mpirank_right=g_d.nb_r;
  g_mpi_ctx.mpirank_down=g_d.nb_d;g_mpi_ctx.mpirank_up=g_d.nb_u;
  g_mpi_ctx.mpi_tag_base=100;
}

static void cfg_load(const char*cp,const char*hp){
  mythread_cfg*hw=mythread_cfg_load(hp),*cs=mythread_cfg_load(cp);
  cfg_NX=cs?mythread_cfg_get_int(cs,"","NX",cfg_NX):cfg_NX;cfg_NY=cs?mythread_cfg_get_int(cs,"","NY",cfg_NY):cfg_NY;
  cfg_NT=cs?mythread_cfg_get_int(cs,"","NT",cfg_NT):cfg_NT;cfg_A=cs?mythread_cfg_get_double(cs,"","A",cfg_A):cfg_A;
  cfg_DT=cs?mythread_cfg_get_double(cs,"","DT",cfg_DT):cfg_DT;cfg_C0=cs?mythread_cfg_get_double(cs,"","C0",cfg_C0):cfg_C0;
  cfg_HALO=cs?mythread_cfg_get_int(cs,"","HALO",cfg_HALO):cfg_HALO;
  cfg_USE_FIXED_DOMAIN=cs?mythread_cfg_get_int(cs,"","USE_FIXED_DOMAIN",cfg_USE_FIXED_DOMAIN):cfg_USE_FIXED_DOMAIN;
  cfg_N_GROUPS=hw?mythread_cfg_get_int(hw,"","N_GROUPS",cfg_N_GROUPS):cfg_N_GROUPS;
  cfg_N_WORKERS=hw?mythread_cfg_get_int(hw,"","N_WORKERS",cfg_N_WORKERS):cfg_N_WORKERS;
  cfg_GROUP_DECOMP=hw?mythread_cfg_get_int(hw,"","GROUP_DECOMP",cfg_GROUP_DECOMP):cfg_GROUP_DECOMP;
  cfg_NCorePClu=hw?mythread_cfg_get_int(hw,"","NCorePClu",cfg_NCorePClu):cfg_NCorePClu;
  cfg_NCluPNode=hw?mythread_cfg_get_int(hw,"","NCluPNode",cfg_NCluPNode):cfg_NCluPNode;
  cfg_NCorePGrp=hw?mythread_cfg_get_int(hw,"","NCorePGrp",cfg_NCorePGrp):cfg_NCorePGrp;
  cfg_ManageCoreId=hw?mythread_cfg_get_int(hw,"","ManageCoreId",cfg_ManageCoreId):cfg_ManageCoreId;
  cfg_CoreOffset=hw?mythread_cfg_get_int(hw,"","CoreOffset",cfg_CoreOffset):cfg_CoreOffset;
  cfg_ClustOffset=hw?mythread_cfg_get_int(hw,"","ClustOffset",cfg_ClustOffset):cfg_ClustOffset;
  mythread_cfg_free(cs);mythread_cfg_free(hw);cfg_compute_derived();
}

static void setup_thread_tasks(void){
  for(int g=0;g<md.ngrp;g++){
    threadGroup*pg=md.grps[g];
    const mythread_tile*t=&g_dc->group_tiles[g];
    int nworkers=pg->Nthreads-1;
    for(int i=0;i<pg->Nthreads;i++){
      ThreadTask*task=(ThreadTask*)pg->threads[i].td;
      task->gid=g;task->tid=i;task->gf=&g_gfields[g];
      task->x_begin=HALO;task->x_end=HALO+t->nx;
      task->y_begin=HALO;task->y_end=HALO;
      task->partial_energy=0.0;
    }
    if(nworkers>0){
      int base=t->ny/nworkers,rem=t->ny%nworkers,y=HALO;
      for(int i=1;i<pg->Nthreads;i++){
        ThreadTask*task=(ThreadTask*)pg->threads[i].td;
        int rows=base+((i-1)<rem?1:0);
        task->y_begin=y;task->y_end=y+rows;y+=rows;
      }
    }
  }
}

static void worker_thread(void){
  ThreadTask*task=(ThreadTask*)ti->td;
  GroupField*gf=task->gf;int gid=task->gid;

  sWaitGrp(init_state());
  init_field_block(gf,gid,task->y_begin,task->y_end,task->x_begin,task->x_end);
  sSetGrp(init_state());

  sWaitGrp(copy_state());
  copy_prev_block(gf,gid,task->y_begin,task->y_end);
  sSetGrp(copy_state());

  for(int step=0;step<NT;step++){
    int cs=compute_state(step),bs=boundary_state(step);
    sWaitGrp(cs);
    comp_interior_block(gf,gid,task->y_begin,task->y_end,task->x_begin,task->x_end);
    sSetGrp(cs);

    sWaitGrp(bs);
    comp_boundary_block(gf,gid,task->y_begin,task->y_end,task->x_begin,task->x_end);
    sSetGrp(bs);
  }

  sWaitGrp(energy_state());
  task->partial_energy=energy_block(gf,gid,task->y_begin,task->y_end,task->x_begin,task->x_end);
  sSetGrp(energy_state());
}

static double collect_worker_energy(void){
  double e=0.0;
  for(int g=0;g<md.ngrp;g++){
    threadGroup*pg=md.grps[g];
    for(int i=1;i<pg->Nthreads;i++){
      ThreadTask*task=(ThreadTask*)pg->threads[i].td;
      e+=task->partial_energy;
    }
  }
  return e;
}

static void group_main_thread(void){
  int gid=ti->igrp;
  ThreadTask*task=(ThreadTask*)ti->td;
  if(group_alloc_field(&g_gfields[gid],gid,g_dc)!=0){
    fprintf(stderr,"[E] rank %d gid %d group_alloc_field failed\n",mpi_id,gid);
    MPI_Abort(MPI_COMM_WORLD,1);
  }
  task->gid=gid;task->tid=ti->ind;task->gf=&g_gfields[gid];

  gWaitMain(init_state());gSetSubs(init_state());gWaitSubs(init_state());gSetMain(init_state());
  gWaitMain(copy_state());gSetSubs(copy_state());gWaitSubs(copy_state());gSetMain(copy_state());

  for(int step=0;step<NT;step++){
    int cs=compute_state(step),bs=boundary_state(step);
    gWaitMain(cs);gSetSubs(cs);gWaitSubs(cs);gSetMain(cs);
    gWaitMain(bs);gSetSubs(bs);gWaitSubs(bs);gSetMain(bs);
  }

  gWaitMain(energy_state());gSetSubs(energy_state());gWaitSubs(energy_state());gSetMain(energy_state());
}

static void main_thread_run(int mr){
  mSetGrps(init_state());
  mWaitGrps(init_state());

  group_field_link_buffers(g_gfields,g_dc);
  apply_d_all();halo_all();apply_d_all();

  mSetGrps(copy_state());
  mWaitGrps(copy_state());

  double t0=MPI_Wtime();
  for(int step=0;step<NT;step++){
    MPI_Barrier(MPI_COMM_WORLD);
    apply_d_all();halo_all();apply_d_all();

    mSetGrps(compute_state(step));
    mWaitGrps(compute_state(step));

    mSetGrps(boundary_state(step));
    mWaitGrps(boundary_state(step));

    gf_swap_all();
    if(mr==0&&(step+1)%((NT>10?NT/10:1))==0)printf("[Main] Step %d/%d\n",step+1,NT);
  }
  double elapsed=MPI_Wtime()-t0;

  mSetGrps(energy_state());
  mWaitGrps(energy_state());
  double my_e=collect_worker_energy();
  double global_e;MPI_Allreduce(&my_e,&global_e,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  if(mr==0){printf("[Main] Done: %.3f s, %.1f Mpts/s\n",elapsed,(double)NX*NY*NT/elapsed/1e6);
    printf("[Main] Final energy: %.6f (expected %.6f) %s\n",global_e,1.570796,fabs(global_e-1.570796)<0.01?"OK":"?");}
}

void thread_run(void){
  if(ti->igrp<0)main_thread_run(mpi_id);
  else if(ti->ind==0)group_main_thread();
  else worker_thread();
}

int main(int argc,char**argv){
  int mr,ms,lr=1,gr=1,err;
  MPI_Init(&argc,&argv);MPI_Comm_rank(MPI_COMM_WORLD,&mr);MPI_Comm_size(MPI_COMM_WORLD,&ms);
  {const char*cp=mythread_env_get("WAVE_CASE_CFG","config/case.cfg"),
        *hp=mythread_env_get("WAVE_HARDWARE_CFG","config/hardware.cfg");cfg_load(cp,hp);}
  if(cfg_N_GROUPS<2)cfg_N_GROUPS=2;
  if(cfg_N_WORKERS<1)cfg_N_WORKERS=1;
  if(cfg_GROUP_DECOMP!=MYTHREAD_DECOMP_Y_ONLY&&mr==0)
    fprintf(stderr,"[W] wave_propagation_ghost_mpi_sync.c GroupField path is optimized for Y_ONLY; forcing GROUP_DECOMP=0\n");
  cfg_GROUP_DECOMP=MYTHREAD_DECOMP_Y_ONLY;
  if(NX<3||NY<3){if(mr==0)fprintf(stderr,"[E] NX/NY>=3\n");MPI_Finalize();return 1;}
  if(CFL_SUM2>1.0){if(mr==0)fprintf(stderr,"[E] CFL>1\n");MPI_Finalize();return 1;}

  setup_domain(mr,ms);
  if(g_d.nx<=0||g_d.ny<=0)lr=0;
  MPI_Allreduce(&lr,&gr,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);if(!gr){MPI_Finalize();return 1;}

  g_dc=mythread_decomp_create(g_d.nx,g_d.ny,HALO,cfg_N_GROUPS,cfg_N_WORKERS,MYTHREAD_DECOMP_Y_ONLY);
  if(!g_dc){fprintf(stderr,"[E] rank %d decomp\n",mr);MPI_Finalize();return 1;}
  g_gfields=(GroupField*)calloc((size_t)g_dc->n_groups,sizeof(GroupField));
  if(!g_gfields){fprintf(stderr,"[E] rank %d group fields\n",mr);mythread_decomp_free((mythread_decomp*)g_dc);MPI_Finalize();return 1;}

  MPI_Comm node_comm;int node_size=1,node_rank=0;
  MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,0,MPI_INFO_NULL,&node_comm);
  MPI_Comm_size(node_comm,&node_size);
  MPI_Comm_rank(node_comm,&node_rank);

  int manage_core=cfg_ManageCoreId;
  CoreOffset=cfg_CoreOffset;ClustOffset=cfg_ClustOffset;
  err=InitThreads(node_rank,cfg_NCorePClu,cfg_NCluPNode,cfg_NCorePGrp,
                  cfg_N_WORKERS+1,cfg_N_GROUPS,node_size,&manage_core);
  if(err!=0){fprintf(stderr,"[E] rank %d InitThreads failed: %d\n",mr,err);
    free(g_gfields);mythread_decomp_free((mythread_decomp*)g_dc);MPI_Comm_free(&node_comm);MPI_Finalize();return 1;}

  setup_thread_tasks();

  if(mr==0)printf("=== Wave MPI+mythread-sync GroupField %dx%dx%d procs=%d (%dx%d) groups/rank=%d workers/group=%d DT=%.4f CFL=%.4f ===\n",
    NX,NY,NT,ms,g_d.px,g_d.py,cfg_N_GROUPS,cfg_N_WORKERS,cfg_DT,cfg_C0*cfg_DT/cfg_DX);
  printf("[R%d] domain %dx%d at (%d,%d) nbr L=%d R=%d D=%d U=%d groups=%d\n",
    mr,g_d.nx,g_d.ny,g_d.mx,g_d.my,g_d.nb_l,g_d.nb_r,g_d.nb_d,g_d.nb_u,g_dc->n_groups);

  StartThreads(thread_run);
  thread_run();
  EndThreads();

  for(int g=0;g<g_dc->n_groups;g++)group_free_field(&g_gfields[g]);
  free(g_gfields);mythread_decomp_free((mythread_decomp*)g_dc);MPI_Comm_free(&node_comm);MPI_Finalize();return 0;
}
