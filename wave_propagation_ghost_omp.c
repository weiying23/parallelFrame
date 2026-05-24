#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#if defined(__linux__)
#include <sys/sysinfo.h>
#endif

#include "mythread/mythread_config.h"
#include "mythread/mythread_decomp.h"

/*
 * MPI + OpenMP — NUMA 分组版
 *
 * 线程数 = N_GROUPS × N_WORKERS
 * 每个组一个 GroupField，绑定到对应 NUMA 节点
 * 组间 halo 通过 memcpy 交换（intra），域边界通过 MPI
 */

static int    cfg_NX=14000, cfg_NY=14000, cfg_NT=480;
static double cfg_A=10.0, cfg_DT=0.001, cfg_C0=0.1;
static int    cfg_USE_FIXED_DOMAIN=0;
static double cfg_DX, cfg_DY, cfg_LX, cfg_LY, cfg_DT2, cfg_CFL_X, cfg_CFL_Y, cfg_CFL_SUM2;
static int    cfg_HALO=1, cfg_N_GROUPS=2, cfg_N_WORKERS=3;
static int    cfg_ENERGY_REPORT_INTERVAL=60, cfg_GROUP_DECOMP=0;

static void cfg_compute_derived(void) {
  if(cfg_USE_FIXED_DOMAIN){cfg_LX=1.0;cfg_LY=1.0;
    cfg_DX=cfg_LX/(double)(cfg_NX-1);cfg_DY=cfg_LY/(double)(cfg_NY-1);}
  else{cfg_DX=0.01;cfg_DY=0.01;
    cfg_LX=(cfg_NX-1)*cfg_DX;cfg_LY=(cfg_NY-1)*cfg_DY;}
  cfg_DT2=cfg_DT*cfg_DT;
  cfg_CFL_X=cfg_C0*cfg_DT/cfg_DX;cfg_CFL_Y=cfg_C0*cfg_DT/cfg_DY;
  cfg_CFL_SUM2=cfg_CFL_X*cfg_CFL_X+cfg_CFL_Y*cfg_CFL_Y;
}

#define NX    cfg_NX
#define NY    cfg_NY
#define NT    cfg_NT
#define A     cfg_A
#define DT    cfg_DT
#define C0    cfg_C0
#define DX    cfg_DX
#define DY    cfg_DY
#define LX    cfg_LX
#define LY    cfg_LY
#define DT2   cfg_DT2
#define CFL_X cfg_CFL_X
#define CFL_Y cfg_CFL_Y
#define CFL_SUM2 cfg_CFL_SUM2
#define HALO  cfg_HALO
#define ENERGY_REPORT_INTERVAL cfg_ENERGY_REPORT_INTERVAL

typedef struct GroupField {
  double *u_prev,*u_curr,*u_next;
  int ny_padded,nx_padded,stride;
  size_t plane_bytes;
  double *send_to_up,*recv_from_up,*send_to_down,*recv_from_down;
  double *send_to_left,*recv_from_left,*send_to_right,*recv_from_right;
  double *send_up,*recv_up,*send_down,*recv_down;
  double *send_left,*recv_left,*send_right,*recv_right;
  double group_energy;
} GroupField;

#define GFIDX(gf,y,x) ((size_t)(y)*(size_t)(gf)->stride+(size_t)(x))

typedef struct {
  int mpirank_up,mpirank_down,mpirank_left,mpirank_right;
  int mpi_tag_base;
} MpiCtx;

typedef struct {
  int local_x_begin,local_x_end,local_nx;
  int local_y_begin,local_y_end,local_ny;
  int mpi_rank,mpi_size;
  int proc_x,proc_y,proc_px,proc_py;
  int neighbor_left,neighbor_right,neighbor_up,neighbor_down;
  double initial_energy;
} SimData;

enum {TM_COMP=0,TM_ENERGY=1,TM_HALO=2,TM_N=8}; /* 8 doubles = 64B cache-line pad */
#define TM_SLOT(tid) ((tid)*TM_N)

static SimData g_sim={0};
static GroupField *g_gfields=NULL;
static const mythread_decomp *g_decomp=NULL;
static MpiCtx g_mpi_ctx={0};
static double *g_l2_acc=NULL,*g_max_acc=NULL;
static int    *g_max_x=NULL,*g_max_y=NULL;
static double *g_times=NULL;
static int g_nthreads,g_n_groups,g_n_workers;

static inline double wall_time(void){return MPI_Wtime();}

static inline int global_x_from_local(int gid,int local_x){
  return g_sim.local_x_begin+(g_decomp->group_tiles[gid].x_begin+local_x-HALO);
}
static inline int global_y_from_local(int gid,int local_y){
  return g_sim.local_y_begin+(g_decomp->group_tiles[gid].y_begin+local_y-HALO);
}

static void*xcalloc(size_t n,size_t s){
  void*p=calloc(n,s);if(!p){fprintf(stderr,"[Error] alloc failed\n");MPI_Abort(MPI_COMM_WORLD,1);}
  return p;
}

static int gf_alloc(GroupField*gf,int gid,const mythread_decomp*dc){
  if(!gf||!dc||gid<0||gid>=dc->n_groups)return-1;
  const mythread_tile*tile=&dc->group_tiles[gid];
  int h=dc->halo,nx=tile->nx,ny=tile->ny;
  memset(gf,0,sizeof(*gf));
  gf->ny_padded=ny+2*h;gf->nx_padded=nx+2*h;gf->stride=gf->nx_padded;
  gf->plane_bytes=(size_t)gf->ny_padded*(size_t)gf->nx_padded*sizeof(double);
#define ALLOC(p,sz) do{(p)=(double*)calloc((sz),sizeof(double));if(!(p))goto fail;}while(0)
  ALLOC(gf->u_prev,(size_t)gf->ny_padded*gf->nx_padded);
  ALLOC(gf->u_curr,(size_t)gf->ny_padded*gf->nx_padded);
  ALLOC(gf->u_next,(size_t)gf->ny_padded*gf->nx_padded);
  if(mythread_decomp_neighbor(dc,gid,MYTHREAD_NEIGHBOR_UP)>=0){ALLOC(gf->recv_from_up,nx);}
  if(mythread_decomp_neighbor(dc,gid,MYTHREAD_NEIGHBOR_DOWN)>=0){ALLOC(gf->recv_from_down,nx);}
  if(mythread_decomp_neighbor(dc,gid,MYTHREAD_NEIGHBOR_LEFT)>=0){ALLOC(gf->recv_from_left,ny);}
  if(mythread_decomp_neighbor(dc,gid,MYTHREAD_NEIGHBOR_RIGHT)>=0){ALLOC(gf->recv_from_right,ny);}
  if(mythread_decomp_is_domain_boundary(dc,gid,MYTHREAD_NEIGHBOR_UP)){ALLOC(gf->send_up,nx);ALLOC(gf->recv_up,nx);}
  if(mythread_decomp_is_domain_boundary(dc,gid,MYTHREAD_NEIGHBOR_DOWN)){ALLOC(gf->send_down,nx);ALLOC(gf->recv_down,nx);}
  if(mythread_decomp_is_domain_boundary(dc,gid,MYTHREAD_NEIGHBOR_LEFT)){ALLOC(gf->send_left,ny);ALLOC(gf->recv_left,ny);}
  if(mythread_decomp_is_domain_boundary(dc,gid,MYTHREAD_NEIGHBOR_RIGHT)){ALLOC(gf->send_right,ny);ALLOC(gf->recv_right,ny);}
#undef ALLOC
  return 0;
fail:
  free(gf->u_prev);free(gf->u_curr);free(gf->u_next);
  free(gf->recv_from_up);free(gf->recv_from_down);free(gf->recv_from_left);free(gf->recv_from_right);
  free(gf->send_up);free(gf->recv_up);free(gf->send_down);free(gf->recv_down);
  free(gf->send_left);free(gf->recv_left);free(gf->send_right);free(gf->recv_right);
  memset(gf,0,sizeof(*gf));return-1;
}
static void gf_free(GroupField*gf){
  if(!gf)return;
  free(gf->u_prev);free(gf->u_curr);free(gf->u_next);
  free(gf->recv_from_up);free(gf->recv_from_down);free(gf->recv_from_left);free(gf->recv_from_right);
  free(gf->send_up);free(gf->recv_up);free(gf->send_down);free(gf->recv_down);
  free(gf->send_left);free(gf->recv_left);free(gf->send_right);free(gf->recv_right);
  memset(gf,0,sizeof(*gf));
}
static void gf_swap(GroupField*gf){
  if(!gf)return;
  double*tmp=gf->u_prev;gf->u_prev=gf->u_curr;gf->u_curr=gf->u_next;gf->u_next=tmp;
}

/* 组间缓冲指针互连（同 fix.c 的 group_field_link_buffers） */
static void gf_link_buffers(GroupField*gfs,const mythread_decomp*dc){
  if(!gfs||!dc)return;
  for(int g=0;g<dc->n_groups;g++){
    int nb=mythread_decomp_neighbor(dc,g,MYTHREAD_NEIGHBOR_UP);
    if(nb>=0)gfs[g].send_to_up=gfs[nb].recv_from_down;
    nb=mythread_decomp_neighbor(dc,g,MYTHREAD_NEIGHBOR_DOWN);
    if(nb>=0)gfs[g].send_to_down=gfs[nb].recv_from_up;
    nb=mythread_decomp_neighbor(dc,g,MYTHREAD_NEIGHBOR_LEFT);
    if(nb>=0)gfs[g].send_to_left=gfs[nb].recv_from_right;
    nb=mythread_decomp_neighbor(dc,g,MYTHREAD_NEIGHBOR_RIGHT);
    if(nb>=0)gfs[g].send_to_right=gfs[nb].recv_from_left;
  }
}

/* 组间 halo 交换（纯 memcpy，无 MPI） */
static void halo_exchange_intra(GroupField*gfs,const mythread_decomp*dc){
  if(!gfs||!dc)return;
  int h=dc->halo,ng=dc->n_groups;
  /* Phase 1: 拷贝源组边界到目标组 recv 缓冲 */
  for(int g=0;g<ng;g++){
    GroupField*gf=&gfs[g];if(!gf->u_curr)continue;
    const mythread_tile*t=&dc->group_tiles[g];
    int nx=t->nx,ny=t->ny,nb;
    nb=mythread_decomp_neighbor(dc,g,MYTHREAD_NEIGHBOR_UP);
    if(nb>=0&&gfs[nb].recv_from_down)
      memcpy(gfs[nb].recv_from_down,&gf->u_curr[ny*gf->stride+h],(size_t)nx*sizeof(double));
    nb=mythread_decomp_neighbor(dc,g,MYTHREAD_NEIGHBOR_DOWN);
    if(nb>=0&&gfs[nb].recv_from_up)
      memcpy(gfs[nb].recv_from_up,&gf->u_curr[h*gf->stride+h],(size_t)nx*sizeof(double));
    nb=mythread_decomp_neighbor(dc,g,MYTHREAD_NEIGHBOR_LEFT);
    if(nb>=0&&gfs[nb].recv_from_right)
      for(int r=0;r<ny;r++)gfs[nb].recv_from_right[r]=gf->u_curr[(h+r)*gf->stride+h];
    nb=mythread_decomp_neighbor(dc,g,MYTHREAD_NEIGHBOR_RIGHT);
    if(nb>=0&&gfs[nb].recv_from_left){
      int sc=h+nx-1;
      for(int r=0;r<ny;r++)gfs[nb].recv_from_left[r]=gf->u_curr[(h+r)*gf->stride+sc];
    }
  }
  /* Phase 2: 应用到目标组 halo 区域 */
  for(int g=0;g<ng;g++){
    GroupField*gf=&gfs[g];if(!gf->u_curr)continue;
    const mythread_tile*t=&dc->group_tiles[g];int nx=t->nx,ny=t->ny;
    if(gf->recv_from_up)
      memcpy(&gf->u_curr[(ny+h)*gf->stride+h],gf->recv_from_up,(size_t)nx*sizeof(double));
    if(gf->recv_from_down)
      memcpy(&gf->u_curr[0*gf->stride+h],gf->recv_from_down,(size_t)nx*sizeof(double));
    if(gf->recv_from_left)
      for(int r=0;r<ny;r++)gf->u_curr[(h+r)*gf->stride+0]=gf->recv_from_left[r];
    if(gf->recv_from_right){
      int dc2=h+nx;
      for(int r=0;r<ny;r++)gf->u_curr[(h+r)*gf->stride+dc2]=gf->recv_from_right[r];
    }
  }
}

/* MPI 域边界 halo 交换 */
static void halo_exchange_mpi(GroupField*gfs,const mythread_decomp*dc){
  if(!gfs||!dc)return;
  int h=dc->halo,tag=g_mpi_ctx.mpi_tag_base;
  MPI_Request reqs[8];int nr=0;
  for(int g=0;g<dc->n_groups;g++){
    GroupField*gf=&gfs[g];if(!gf->u_curr)continue;
    const mythread_tile*t=&dc->group_tiles[g];int nx=t->nx,ny=t->ny,s=gf->stride;
    /* Y up */
    if(g_sim.neighbor_up>=0&&gf->send_up&&gf->recv_up&&
       mythread_decomp_is_domain_boundary(dc,g,MYTHREAD_NEIGHBOR_UP)){
      memcpy(gf->send_up,&gf->u_curr[ny*s+h],(size_t)nx*sizeof(double));
      MPI_Isend(gf->send_up,nx,MPI_DOUBLE,g_sim.neighbor_up,tag,MPI_COMM_WORLD,&reqs[nr++]);
      MPI_Irecv(gf->recv_up,nx,MPI_DOUBLE,g_sim.neighbor_up,tag+1,MPI_COMM_WORLD,&reqs[nr++]);
    }
    /* Y down */
    if(g_sim.neighbor_down>=0&&gf->send_down&&gf->recv_down&&
       mythread_decomp_is_domain_boundary(dc,g,MYTHREAD_NEIGHBOR_DOWN)){
      memcpy(gf->send_down,&gf->u_curr[h*s+h],(size_t)nx*sizeof(double));
      MPI_Isend(gf->send_down,nx,MPI_DOUBLE,g_sim.neighbor_down,tag+1,MPI_COMM_WORLD,&reqs[nr++]);
      MPI_Irecv(gf->recv_down,nx,MPI_DOUBLE,g_sim.neighbor_down,tag,MPI_COMM_WORLD,&reqs[nr++]);
    }
  }
  MPI_Waitall(nr,reqs,MPI_STATUSES_IGNORE);
  /* Apply recv to halo */
  for(int g=0;g<dc->n_groups;g++){
    GroupField*gf=&gfs[g];if(!gf->u_curr)continue;
    const mythread_tile*t=&dc->group_tiles[g];int nx=t->nx,ny=t->ny,s=gf->stride;
    if(gf->recv_up) memcpy(&gf->u_curr[(ny+h)*s+h],gf->recv_up,(size_t)nx*sizeof(double));
    if(gf->recv_down) memcpy(&gf->u_curr[0*s+h],gf->recv_down,(size_t)nx*sizeof(double));
  }
}

static void apply_dirichlet(GroupField*gf,int gid){
  const mythread_tile*tile=&g_decomp->group_tiles[gid];
  int s=gf->stride,h=gf->ny_padded,nx=tile->nx,ny=tile->ny;
  if(g_sim.local_x_begin==0)
    for(int y=0;y<h;y++)gf->u_curr[GFIDX(gf,y,HALO)]=0.0;
  if(g_sim.local_x_end==NX){int xr=HALO+nx-1;
    for(int y=0;y<h;y++)gf->u_curr[GFIDX(gf,y,xr)]=0.0;}
  if(g_sim.neighbor_left<0&&mythread_decomp_is_domain_boundary(g_decomp,gid,MYTHREAD_NEIGHBOR_LEFT))
    for(int y=0;y<h;y++)gf->u_curr[GFIDX(gf,y,0)]=0.0;
  if(g_sim.neighbor_right<0&&mythread_decomp_is_domain_boundary(g_decomp,gid,MYTHREAD_NEIGHBOR_RIGHT))
  {int xh=HALO+nx;for(int y=0;y<h;y++)gf->u_curr[GFIDX(gf,y,xh)]=0.0;}
  if(g_sim.neighbor_down<0&&mythread_decomp_is_domain_boundary(g_decomp,gid,MYTHREAD_NEIGHBOR_DOWN)){
    memset(&gf->u_curr[0],0,(size_t)s*sizeof(double));
    if(g_sim.local_y_begin==0)memset(&gf->u_curr[GFIDX(gf,HALO,0)],0,(size_t)s*sizeof(double));
  }
  if(g_sim.neighbor_up<0&&mythread_decomp_is_domain_boundary(g_decomp,gid,MYTHREAD_NEIGHBOR_UP)){
    memset(&gf->u_curr[(ny+HALO)*s],0,(size_t)s*sizeof(double));
    if(g_sim.local_y_end==NY)memset(&gf->u_curr[GFIDX(gf,ny,0)],0,(size_t)s*sizeof(double));
  }
}
static void apply_dirichlet_all(void){
  for(int g=0;g<g_decomp->n_groups;g++)apply_dirichlet(&g_gfields[g],g);
}
static void copy_curr_to_prev_all(void){
  for(int g=0;g<g_decomp->n_groups;g++)
    memcpy(g_gfields[g].u_prev,g_gfields[g].u_curr,g_gfields[g].plane_bytes);
}

static void compute_region(GroupField*gf,int gid,int yb,int ye,int xb,int xe){
  for(int y=yb;y<ye;y++){int gy=global_y_from_local(gid,y);if(gy==0||gy==NY-1)continue;
    for(int x=xb;x<xe;x++){int gx=global_x_from_local(gid,x);if(gx==0||gx==NX-1)continue;
      double u=gf->u_curr[GFIDX(gf,y,x)];
      double d2x=(gf->u_curr[GFIDX(gf,y,x-1)]-2*u+gf->u_curr[GFIDX(gf,y,x+1)])/(DX*DX);
      double d2y=(gf->u_curr[GFIDX(gf,y-1,x)]-2*u+gf->u_curr[GFIDX(gf,y+1,x)])/(DY*DY);
      gf->u_next[GFIDX(gf,y,x)]=2*u-gf->u_prev[GFIDX(gf,y,x)]+C0*C0*DT2*(d2x+d2y);
    }
  }
}
static void compute_interior(GroupField*gf,int gid,int yb,int ye,int xb,int xe){
  const mythread_tile*t=&g_decomp->group_tiles[gid];int ny=t->ny,nx=t->nx;
  if(g_sim.neighbor_down>=0&&mythread_decomp_is_domain_boundary(g_decomp,gid,MYTHREAD_NEIGHBOR_DOWN)&&yb<HALO+1)yb=HALO+1;
  if(g_sim.neighbor_up>=0&&mythread_decomp_is_domain_boundary(g_decomp,gid,MYTHREAD_NEIGHBOR_UP)&&ye>ny)ye=ny;
  if(g_sim.neighbor_left>=0&&mythread_decomp_is_domain_boundary(g_decomp,gid,MYTHREAD_NEIGHBOR_LEFT)&&xb<HALO+1)xb=HALO+1;
  if(g_sim.neighbor_right>=0&&mythread_decomp_is_domain_boundary(g_decomp,gid,MYTHREAD_NEIGHBOR_RIGHT)&&xe>nx)xe=nx;
  if(yb<ye&&xb<xe)compute_region(gf,gid,yb,ye,xb,xe);
}
static void compute_boundary(GroupField*gf,int gid,int yb,int ye,int xb,int xe){
  const mythread_tile*t=&g_decomp->group_tiles[gid];int ny=t->ny,nx=t->nx;
  int lr=HALO,ur=ny,lc=HALO,rc=nx;
  if(xb<HALO)xb=HALO;if(xe>HALO+nx)xe=HALO+nx;
  if(yb<HALO)yb=HALO;if(ye>HALO+ny)ye=HALO+ny;
  int nd=g_sim.neighbor_down>=0&&mythread_decomp_is_domain_boundary(g_decomp,gid,MYTHREAD_NEIGHBOR_DOWN);
  int nu=g_sim.neighbor_up>=0&&mythread_decomp_is_domain_boundary(g_decomp,gid,MYTHREAD_NEIGHBOR_UP);
  int nl=g_sim.neighbor_left>=0&&mythread_decomp_is_domain_boundary(g_decomp,gid,MYTHREAD_NEIGHBOR_LEFT);
  int nr=g_sim.neighbor_right>=0&&mythread_decomp_is_domain_boundary(g_decomp,gid,MYTHREAD_NEIGHBOR_RIGHT);
  if(nd&&yb<=lr&&lr<ye)compute_region(gf,gid,lr,lr+1,xb,xe);
  if(nu&&ur!=lr&&yb<=ur&&ur<ye)compute_region(gf,gid,ur,ur+1,xb,xe);
  if(nl&&xb<=lc&&lc<xe)compute_region(gf,gid,yb,ye,lc,lc+1);
  if(nr&&rc!=lc&&xb<=rc&&rc<xe)compute_region(gf,gid,yb,ye,rc,rc+1);
}

static double compute_energy_block(GroupField*gf,int gid,int yb,int ye,int xb,int xe,
    double*l2_out,double*max_out,int*max_x,int*max_y){
  if(xb>=xe||yb>=ye){if(l2_out)*l2_out=0.0;if(max_out)*max_out=0.0;return 0.0;}
  const mythread_tile*t=&g_decomp->group_tiles[gid];
  double ke=0,px=0,py=0,ca=DX*DY,l2=0,ma=0;int mx=0,my=0;
  int xeb=(xb==HALO)?(HALO-1):xb,xee=xe;if(xee>HALO+t->nx)xee=HALO+t->nx;
  for(int y=yb;y<ye;y++){int gy=global_y_from_local(gid,y);
    for(int x=xb;x<xe;x++){double u=gf->u_curr[GFIDX(gf,y,x)];l2+=u*u;
      double au=fabs(u);if(au>ma){ma=au;mx=global_x_from_local(gid,x);my=gy;}}
    if(gy>0&&gy<NY-1)for(int x=xb;x<xe;x++)
      {double ut=(gf->u_curr[GFIDX(gf,y,x)]-gf->u_prev[GFIDX(gf,y,x)])/DT;ke+=ut*ut;}
    for(int x=xeb;x<xee;x++)
      {double dc=(gf->u_curr[GFIDX(gf,y,x+1)]-gf->u_curr[GFIDX(gf,y,x)])/DX;
       double dp=(gf->u_prev[GFIDX(gf,y,x+1)]-gf->u_prev[GFIDX(gf,y,x)])/DX;px+=dc*dp;}
    if(gy<NY-1)for(int x=xb;x<xe;x++)
      {double dc=(gf->u_curr[GFIDX(gf,y+1,x)]-gf->u_curr[GFIDX(gf,y,x)])/DY;
       double dp=(gf->u_prev[GFIDX(gf,y+1,x)]-gf->u_prev[GFIDX(gf,y,x)])/DY;py+=dc*dp;}
  }
  if(l2_out)*l2_out=l2;if(max_out)*max_out=ma;if(max_x)*max_x=mx;if(max_y)*max_y=my;
  return 0.5*(ke+C0*C0*(px+py))*ca;
}
static int should_measure_energy(int step){
  if(step==0||step==NT-1)return 1;
  if(ENERGY_REPORT_INTERVAL>0&&((step+1)%ENERGY_REPORT_INTERVAL)==0)return 1;
  return 0;
}

static double initial_condition_value(int gy,int gx);

/* ═══════════════════════════════════════════════════════════════
 *  NUMA 分组版 run_simulation
 *  线程数 = n_groups × n_workers
 *  每个线程: gid = tid / n_workers, 负责该组的 y 子区间
 * ═══════════════════════════════════════════════════════════════ */
static void run_simulation(void) {
  for(int g=0;g<g_n_groups;g++)
    if(gf_alloc(&g_gfields[g],g,g_decomp)!=0)
    {fprintf(stderr,"[Error] gf_alloc group %d\n",g);MPI_Abort(MPI_COMM_WORLD,1);}
  gf_link_buffers(g_gfields,g_decomp);

  double *thread_pe=(double*)calloc((size_t)g_nthreads,sizeof(double));
  double *thread_l2=(double*)calloc((size_t)g_nthreads,sizeof(double));
  double *thread_ma=(double*)calloc((size_t)g_nthreads,sizeof(double));
  int    *thread_mx=(int*)calloc((size_t)g_nthreads,sizeof(int));
  int    *thread_my=(int*)calloc((size_t)g_nthreads,sizeof(int));

#pragma omp parallel num_threads(g_nthreads)
  {
    int tid=omp_get_thread_num();
    int gid=tid/g_n_workers;         /* 组 */
    int wid=tid%g_n_workers;         /* 组内 worker */
    double*tm=&g_times[TM_SLOT(tid)];
    double t0;
    GroupField*gf=&g_gfields[gid];
    const mythread_tile*tile=&g_decomp->group_tiles[gid];
    int ny=tile->ny, nx=tile->nx;

    /* 组内 Y 切分：n_workers 个线程均分 ny 行 */
    int base=ny/g_n_workers, rem=ny%g_n_workers;
    int y_start=HALO+wid*base+(wid<rem?wid:rem);
    int y_end=y_start+base+(wid<rem?1:0);

    /* ── Phase 1: 初始化波场（各线程用自己的 y 范围）── */
#pragma omp barrier
#pragma omp single
    t0=wall_time();
    for(int y=y_start;y<y_end;y++){
      int gy=global_y_from_local(gid,y);
      for(int x=HALO;x<HALO+nx;x++){
        int gx=global_x_from_local(gid,x);
        size_t p=GFIDX(gf,y,x);
        double v=(gx!=0&&gx!=NX-1&&gy!=0&&gy!=NY-1)?initial_condition_value(gy,gx):0.0;
        gf->u_curr[p]=v;gf->u_prev[p]=v;gf->u_next[p]=0.0;
      }
    }
#pragma omp single
    tm[TM_COMP]+=wall_time()-t0;

    /* Dirichlet + Halo */
#pragma omp single
    { apply_dirichlet_all();
      t0=wall_time();halo_exchange_intra(g_gfields,g_decomp);
      halo_exchange_mpi(g_gfields,g_decomp);tm[TM_HALO]+=wall_time()-t0;
      apply_dirichlet_all();copy_curr_to_prev_all(); }

    /* ── Phase 2: 初始能量 ── */
#pragma omp single
    { t0=wall_time();g_l2_acc[0]=0.0;g_max_acc[0]=0.0;
      for(int g=0;g<g_n_groups;g++)g_gfields[g].group_energy=0.0; }
#pragma omp barrier
    { thread_pe[tid]=0.0;thread_l2[tid]=0.0;double ma=0.0;int mx=0,my=0;
      for(int y=y_start;y<y_end;y++){
        double pe_,l2_,ma_;int mx_,my_;
        pe_=compute_energy_block(gf,gid,y,y+1,HALO,HALO+nx,&l2_,&ma_,&mx_,&my_);
        thread_pe[tid]+=pe_;thread_l2[tid]+=l2_;
        if(ma_>ma){ma=ma_;mx=mx_;my=my_;}
      }
      thread_ma[tid]=ma;thread_mx[tid]=mx;thread_my[tid]=my; }
#pragma omp barrier
#pragma omp single
    { double pe=0,l2=0,ma=0;int mx=0,my=0;
      for(int t=0;t<g_nthreads;t++){pe+=thread_pe[t];l2+=thread_l2[t];
       if(thread_ma[t]>ma){ma=thread_ma[t];mx=thread_mx[t];my=thread_my[t];}}
      for(int g=0;g<g_n_groups;g++)g_gfields[g].group_energy=0.0;
      g_gfields[0].group_energy=pe;g_l2_acc[0]=l2;
      g_max_acc[0]=ma;g_max_x[0]=mx;g_max_y[0]=my;
      tm[TM_ENERGY]+=wall_time()-t0;
      double local_e=pe,global_e;
      MPI_Allreduce(&local_e,&global_e,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      g_sim.initial_energy=global_e;
      if(g_sim.mpi_rank==0)printf("[Main] Initial: E=%.6f L2=%.6f max|u|=%.6f@(%d,%d)\n",
        global_e,sqrt(g_l2_acc[0]*DX*DY),g_max_acc[0],g_max_x[0],g_max_y[0]); }

    /* ════════════ 时间步循环 ════════════ */
    double prev_time=wall_time();
    for(int step=0;step<NT;step++){
      MPI_Barrier(MPI_COMM_WORLD);
      int need_e=should_measure_energy(step);

      /* (a) Dirichlet + Halo */
#pragma omp single
      { apply_dirichlet_all();
        t0=wall_time();halo_exchange_intra(g_gfields,g_decomp);
        halo_exchange_mpi(g_gfields,g_decomp);tm[TM_HALO]+=wall_time()-t0;
        apply_dirichlet_all(); }

      /* (b) Compute interior（每线程用自己的 y 范围）*/
#pragma omp single
      t0=wall_time();
      for(int y=y_start;y<y_end;y++)
        compute_interior(gf,gid,y,y+1,HALO,HALO+nx);
#pragma omp single
      tm[TM_COMP]+=wall_time()-t0;

      /* (c) Compute boundary */
#pragma omp single
      t0=wall_time();
      for(int y=y_start;y<y_end;y++)
        compute_boundary(gf,gid,y,y+1,HALO,HALO+nx);
#pragma omp single
      tm[TM_COMP]+=wall_time()-t0;

      /* (d) Swap */
#pragma omp single
      for(int g=0;g<g_n_groups;g++)gf_swap(&g_gfields[g]);

      /* (e) Energy */
      if(need_e){
#pragma omp single
        { apply_dirichlet_all();
          t0=wall_time();halo_exchange_intra(g_gfields,g_decomp);
          halo_exchange_mpi(g_gfields,g_decomp);tm[TM_HALO]+=wall_time()-t0;
          apply_dirichlet_all(); }
#pragma omp barrier
        { thread_pe[tid]=0.0;thread_l2[tid]=0.0;double ma=0.0;int mx=0,my=0;
          for(int y=y_start;y<y_end;y++){
            double pe_,l2_,ma_;int mx_,my_;
            pe_=compute_energy_block(gf,gid,y,y+1,HALO,HALO+nx,&l2_,&ma_,&mx_,&my_);
            thread_pe[tid]+=pe_;thread_l2[tid]+=l2_;
            if(ma_>ma){ma=ma_;mx=mx_;my=my_;}
          }
          thread_ma[tid]=ma;thread_mx[tid]=mx;thread_my[tid]=my; }
#pragma omp barrier
#pragma omp single
        { double pe=0,l2=0,ma=0;int mx=0,my=0;
          for(int t=0;t<g_nthreads;t++){pe+=thread_pe[t];l2+=thread_l2[t];
           if(thread_ma[t]>ma){ma=thread_ma[t];mx=thread_mx[t];my=thread_my[t];}}
          g_l2_acc[0]=l2;g_max_acc[0]=ma;g_max_x[0]=mx;g_max_y[0]=my;
          tm[TM_ENERGY]+=wall_time()-t0;
          double local_e=pe,global_e;
          MPI_Allreduce(&local_e,&global_e,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
          if(g_sim.mpi_rank==0){double cur_time=MPI_Wtime();
            printf("[Main] Step %4d/%d, time %.3f,  E=%.6f L2=%.6f max|u|=%.6f@(%d,%d)\n",
              step+1,NT,cur_time-prev_time,global_e,sqrt(g_l2_acc[0]*DX*DY),
              g_max_acc[0],g_max_x[0],g_max_y[0]);prev_time=cur_time;}
        }
      }
    }
  }

  free(thread_pe);free(thread_l2);free(thread_ma);free(thread_mx);free(thread_my);
  double elapsed=wall_time();
  if(g_sim.mpi_rank==0){
    printf("[Main] Simulation completed in %.3f seconds\n",elapsed);
    printf("[Main] Throughput: %.2f Mpoint-updates/s\n",(double)NX*NY*NT/elapsed/1.0e6);
  }
}

static void setup_process_domain(int mpi_rank,int mpi_size){
  mythread_decomp*tmp=mythread_decomp_create(NX,NY,HALO,mpi_size,0,MYTHREAD_DECOMP_XY_2D);
  int px=tmp->gx,py=tmp->gy;mythread_decomp_free(tmp);
  g_sim.proc_px=px;g_sim.proc_py=py;
  g_sim.proc_x=mpi_rank%px;g_sim.proc_y=mpi_rank/px;
  {int t_nx=NX/px,r_nx=NX%px,t_ny=NY/py,r_ny=NY%py;
   int ox=g_sim.proc_x*t_nx+(g_sim.proc_x<r_nx?g_sim.proc_x:r_nx);
   int oy=g_sim.proc_y*t_ny+(g_sim.proc_y<r_ny?g_sim.proc_y:r_ny);
   g_sim.local_x_begin=ox;g_sim.local_x_end=ox+t_nx+(g_sim.proc_x<r_nx?1:0);
   g_sim.local_y_begin=oy;g_sim.local_y_end=oy+t_ny+(g_sim.proc_y<r_ny?1:0);}
  g_sim.local_nx=g_sim.local_x_end-g_sim.local_x_begin;
  g_sim.local_ny=g_sim.local_y_end-g_sim.local_y_begin;
  g_sim.neighbor_left=(g_sim.proc_x>0)?(mpi_rank-1):-1;
  g_sim.neighbor_right=(g_sim.proc_x+1<px)?(mpi_rank+1):-1;
  g_sim.neighbor_down=(g_sim.proc_y>0)?(mpi_rank-px):-1;
  g_sim.neighbor_up=(g_sim.proc_y+1<py)?(mpi_rank+px):-1;
  g_mpi_ctx.mpirank_up=g_sim.neighbor_up;g_mpi_ctx.mpirank_down=g_sim.neighbor_down;
  g_mpi_ctx.mpirank_left=g_sim.neighbor_left;g_mpi_ctx.mpirank_right=g_sim.neighbor_right;
  g_mpi_ctx.mpi_tag_base=100;
}
static double initial_condition_value(int gy,int gx){
  double cx=0.5*LX,cy=0.5*LY,sigma=0.06*((LX<LY)?LX:LY);
  return A*exp(-((gx*DX-cx)*(gx*DX-cx)+(gy*DY-cy)*(gy*DY-cy))/(2.*sigma*sigma));
}
static void cfg_load(const char*cp,const char*hp){
  mythread_cfg*hw=mythread_cfg_load(hp),*cs=mythread_cfg_load(cp);
  cfg_NX=cs?mythread_cfg_get_int(cs,"","NX",cfg_NX):cfg_NX;
  cfg_NY=cs?mythread_cfg_get_int(cs,"","NY",cfg_NY):cfg_NY;
  cfg_NT=cs?mythread_cfg_get_int(cs,"","NT",cfg_NT):cfg_NT;
  cfg_A=cs?mythread_cfg_get_double(cs,"","A",cfg_A):cfg_A;
  cfg_DT=cs?mythread_cfg_get_double(cs,"","DT",cfg_DT):cfg_DT;
  cfg_C0=cs?mythread_cfg_get_double(cs,"","C0",cfg_C0):cfg_C0;
  cfg_USE_FIXED_DOMAIN=cs?mythread_cfg_get_int(cs,"","USE_FIXED_DOMAIN",cfg_USE_FIXED_DOMAIN):cfg_USE_FIXED_DOMAIN;
  cfg_HALO=cs?mythread_cfg_get_int(cs,"","HALO",cfg_HALO):cfg_HALO;
  cfg_ENERGY_REPORT_INTERVAL=cs?mythread_cfg_get_int(cs,"","ENERGY_REPORT_INTERVAL",cfg_ENERGY_REPORT_INTERVAL):cfg_ENERGY_REPORT_INTERVAL;
  cfg_N_GROUPS=hw?mythread_cfg_get_int(hw,"","N_GROUPS",cfg_N_GROUPS):cfg_N_GROUPS;
  cfg_N_WORKERS=hw?mythread_cfg_get_int(hw,"","N_WORKERS",cfg_N_WORKERS):cfg_N_WORKERS;
  cfg_GROUP_DECOMP=hw?mythread_cfg_get_int(hw,"","GROUP_DECOMP",cfg_GROUP_DECOMP):cfg_GROUP_DECOMP;
  mythread_cfg_free(cs);mythread_cfg_free(hw);
  cfg_compute_derived();
}
static void print_timing_report(int mr,int ms,int nt){
  static const char*nm[]={"comp","energy","halo"};
  for(int r=0;r<ms;r++){MPI_Barrier(MPI_COMM_WORLD);
    if(mr==r){for(int t=0;t<nt;t++){double*tm=&g_times[TM_SLOT(t)];
      printf("[Timing] rank %d/%d tid=%d gid=%d ",mr,ms,t,t/g_n_workers);
      for(int s=0;s<3;s++)printf("%s=%.3f ",nm[s],tm[s]);
      printf("\n");}fflush(stdout);}}
  MPI_Barrier(MPI_COMM_WORLD);
}

int main(int argc,char**argv){
  int mpi_rank,mpi_size,node_size,local_ready=1,global_ready=1;
  int provided;
  MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE,&provided);
  MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
  MPI_Comm nc;MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,0,MPI_INFO_NULL,&nc);
  MPI_Comm_size(nc,&node_size);

  {const char*cp=mythread_env_get("WAVE_CASE_CFG","config/case.cfg"),
        *hp=mythread_env_get("WAVE_HARDWARE_CFG","config/hardware.cfg");
   cfg_load(cp,hp);}
  g_n_groups=cfg_N_GROUPS>0?cfg_N_GROUPS:1;
  g_n_workers=cfg_N_WORKERS>0?cfg_N_WORKERS:omp_get_max_threads();
  g_nthreads=g_n_groups*g_n_workers;
  omp_set_num_threads(g_nthreads);

  if(NX<3||NY<3){if(mpi_rank==0)fprintf(stderr,"[Error] NX/NY>=3\n");MPI_Finalize();return 1;}
  if(CFL_SUM2>1.0){if(mpi_rank==0)fprintf(stderr,"[Error] CFL>1\n");MPI_Finalize();return 1;}

  memset(&g_sim,0,sizeof(g_sim));g_sim.mpi_rank=mpi_rank;g_sim.mpi_size=mpi_size;
  setup_process_domain(mpi_rank,mpi_size);
  if(g_sim.local_nx<=0||g_sim.local_ny<=0)local_ready=0;
  MPI_Allreduce(&local_ready,&global_ready,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  if(!global_ready){MPI_Finalize();return 1;}

  g_decomp=mythread_decomp_create(g_sim.local_nx,g_sim.local_ny,HALO,g_n_groups,g_n_workers,cfg_GROUP_DECOMP);
  if(!g_decomp){fprintf(stderr,"[Error] decomp_create\n");MPI_Finalize();return 1;}

  g_gfields=(GroupField*)xcalloc((size_t)g_n_groups,sizeof(GroupField));
  g_l2_acc=(double*)xcalloc((size_t)g_n_groups,sizeof(double));
  g_max_acc=(double*)xcalloc((size_t)g_n_groups,sizeof(double));
  g_max_x=(int*)xcalloc((size_t)g_n_groups,sizeof(int));
  g_max_y=(int*)xcalloc((size_t)g_n_groups,sizeof(int));
  g_times=(double*)xcalloc((size_t)g_nthreads*TM_N,sizeof(double));

  if(mpi_rank==0){
    printf("============================================\n");
    printf("  Wave Equation (MPI+OpenMP, NUMA groups)\n");
    printf("============================================\n");
    printf("Global: %dx%d  steps=%d  DT=%.4f CFL=%.4f\n",NX,NY,NT,cfg_DT,CFL_X);
    printf("MPI: %d  groups: %d  workers/group: %d  total-threads: %d\n",
           mpi_size,g_n_groups,g_n_workers,g_nthreads);
    printf("Decomp: %s (%dx%d groups)\n",
           g_decomp->policy==MYTHREAD_DECOMP_Y_ONLY?"Y_ONLY":"XY_2D",g_decomp->gx,g_decomp->gy);
    printf("============================================\n");
  }
  MPI_Barrier(MPI_COMM_WORLD);
  printf("[Domain] rank %d/%d proc=(%d,%d)/(%d,%d) local=%dx%d\n",
    mpi_rank,mpi_size,g_sim.proc_x,g_sim.proc_y,g_sim.proc_px,g_sim.proc_py,
    g_sim.local_nx,g_sim.local_ny);
  MPI_Barrier(MPI_COMM_WORLD);

  double t_start=MPI_Wtime();
  run_simulation();
  double t_elapsed=MPI_Wtime()-t_start;
  print_timing_report(mpi_rank,mpi_size,g_nthreads);
  if(mpi_rank==0)printf("[Main] Total: %.3f s\n",t_elapsed);

  for(int g=0;g<g_n_groups;g++)gf_free(&g_gfields[g]);
  free(g_gfields);mythread_decomp_free((mythread_decomp*)g_decomp);
  free(g_l2_acc);free(g_max_acc);free(g_max_x);free(g_max_y);free(g_times);
  MPI_Finalize();return 0;
}
