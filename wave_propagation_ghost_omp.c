#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mythread/mythread_config.h"
#include "mythread/mythread_decomp.h"

static int    cfg_NX=14000,cfg_NY=14000,cfg_NT=480;
static double cfg_A=10.0,cfg_DT=0.001,cfg_C0=0.1;
static int    cfg_USE_FIXED_DOMAIN=0;
static double cfg_DX,cfg_DY,cfg_LX,cfg_LY,cfg_DT2,cfg_CFL_X,cfg_CFL_Y,cfg_CFL_SUM2;
static int    cfg_HALO=1,cfg_N_GROUPS=2,cfg_N_WORKERS=3,cfg_GROUP_DECOMP=0;

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

typedef struct{ double *u_prev,*u_curr,*u_next; int ny_padded,nx_padded,stride; size_t plane_bytes;
  double *send_to_up,*recv_from_up,*send_to_down,*recv_from_down;
  double *send_to_left,*recv_from_left,*send_to_right,*recv_from_right;
  double *send_up,*recv_up,*send_down,*recv_down,*send_left,*recv_left,*send_right,*recv_right;
} GroupField;
#define GFIDX(gf,y,x) ((size_t)(y)*(size_t)(gf)->stride+(size_t)(x))

typedef struct{ int mpirank_up,mpirank_down,mpirank_left,mpirank_right,mpi_tag_base; } MpiCtx;
typedef struct{ int local_x_begin,local_x_end,local_nx,local_y_begin,local_y_end,local_ny;
  int mpi_rank,mpi_size,proc_x,proc_y,proc_px,proc_py;
  int neighbor_left,neighbor_right,neighbor_up,neighbor_down; } SimData;

static SimData g_sim={0};
static GroupField *g_gfields=NULL;
static const mythread_decomp *g_decomp=NULL;
static MpiCtx g_mpi_ctx={0};
static int g_nthreads;

static inline double wall_time(void){return MPI_Wtime();}
static inline int gx_loc(int gid,int lx){return g_sim.local_x_begin+(g_decomp->group_tiles[gid].x_begin+lx-HALO);}
static inline int gy_loc(int gid,int ly){return g_sim.local_y_begin+(g_decomp->group_tiles[gid].y_begin+ly-HALO);}

static void*xc(size_t n,size_t s){void*p=calloc(n,s);if(!p){fprintf(stderr,"[Error] alloc\n");MPI_Abort(MPI_COMM_WORLD,1);}return p;}

static int gf_alloc(GroupField*gf,int gid,const mythread_decomp*dc){
  if(!gf||!dc||gid<0||gid>=dc->n_groups)return-1;
  const mythread_tile*t=&dc->group_tiles[gid];int h=dc->halo,nx=t->nx,ny=t->ny;
  memset(gf,0,sizeof(*gf));gf->ny_padded=ny+2*h;gf->nx_padded=nx+2*h;gf->stride=gf->nx_padded;
  gf->plane_bytes=(size_t)gf->ny_padded*(size_t)gf->nx_padded*sizeof(double);
#define AL(p,sz) do{(p)=(double*)calloc((sz),sizeof(double));if(!(p))goto fail;}while(0)
  AL(gf->u_prev,(size_t)gf->ny_padded*gf->nx_padded);AL(gf->u_curr,(size_t)gf->ny_padded*gf->nx_padded);
  AL(gf->u_next,(size_t)gf->ny_padded*gf->nx_padded);
  if(mythread_decomp_neighbor(dc,gid,MYTHREAD_NEIGHBOR_UP)>=0){AL(gf->recv_from_up,nx);}
  if(mythread_decomp_neighbor(dc,gid,MYTHREAD_NEIGHBOR_DOWN)>=0){AL(gf->recv_from_down,nx);}
  if(mythread_decomp_neighbor(dc,gid,MYTHREAD_NEIGHBOR_LEFT)>=0){AL(gf->recv_from_left,ny);}
  if(mythread_decomp_neighbor(dc,gid,MYTHREAD_NEIGHBOR_RIGHT)>=0){AL(gf->recv_from_right,ny);}
  if(mythread_decomp_is_domain_boundary(dc,gid,MYTHREAD_NEIGHBOR_UP)){AL(gf->send_up,nx);AL(gf->recv_up,nx);}
  if(mythread_decomp_is_domain_boundary(dc,gid,MYTHREAD_NEIGHBOR_DOWN)){AL(gf->send_down,nx);AL(gf->recv_down,nx);}
  if(mythread_decomp_is_domain_boundary(dc,gid,MYTHREAD_NEIGHBOR_LEFT)){AL(gf->send_left,ny);AL(gf->recv_left,ny);}
  if(mythread_decomp_is_domain_boundary(dc,gid,MYTHREAD_NEIGHBOR_RIGHT)){AL(gf->send_right,ny);AL(gf->recv_right,ny);}
#undef AL
  return 0;
fail:
  free(gf->u_prev);free(gf->u_curr);free(gf->u_next);
  free(gf->recv_from_up);free(gf->recv_from_down);free(gf->recv_from_left);free(gf->recv_from_right);
  free(gf->send_up);free(gf->recv_up);free(gf->send_down);free(gf->recv_down);
  free(gf->send_left);free(gf->recv_left);free(gf->send_right);free(gf->recv_right);
  memset(gf,0,sizeof(*gf));return-1;
}
static void gf_free(GroupField*gf){if(!gf)return;
  free(gf->u_prev);free(gf->u_curr);free(gf->u_next);
  free(gf->recv_from_up);free(gf->recv_from_down);free(gf->recv_from_left);free(gf->recv_from_right);
  free(gf->send_up);free(gf->recv_up);free(gf->send_down);free(gf->recv_down);
  free(gf->send_left);free(gf->recv_left);free(gf->send_right);free(gf->recv_right);memset(gf,0,sizeof(*gf));}
static void gf_swap(GroupField*gf){if(!gf)return;
  double*tmp=gf->u_prev;gf->u_prev=gf->u_curr;gf->u_curr=gf->u_next;gf->u_next=tmp;}

static void gf_link(GroupField*gfs,const mythread_decomp*dc){
  if(!gfs||!dc)return;
  for(int g=0;g<dc->n_groups;g++){int nb;
    nb=mythread_decomp_neighbor(dc,g,MYTHREAD_NEIGHBOR_UP);if(nb>=0)gfs[g].send_to_up=gfs[nb].recv_from_down;
    nb=mythread_decomp_neighbor(dc,g,MYTHREAD_NEIGHBOR_DOWN);if(nb>=0)gfs[g].send_to_down=gfs[nb].recv_from_up;
    nb=mythread_decomp_neighbor(dc,g,MYTHREAD_NEIGHBOR_LEFT);if(nb>=0)gfs[g].send_to_left=gfs[nb].recv_from_right;
    nb=mythread_decomp_neighbor(dc,g,MYTHREAD_NEIGHBOR_RIGHT);if(nb>=0)gfs[g].send_to_right=gfs[nb].recv_from_left;
  }
}

static void halo_intra(GroupField*gfs,const mythread_decomp*dc){
  if(!gfs||!dc)return;int h=dc->halo,ng=dc->n_groups;
  for(int g=0;g<ng;g++){GroupField*gf=&gfs[g];if(!gf->u_curr)continue;
    const mythread_tile*t=&dc->group_tiles[g];int nx=t->nx,ny=t->ny,nb;
    nb=mythread_decomp_neighbor(dc,g,MYTHREAD_NEIGHBOR_UP);
    if(nb>=0&&gfs[nb].recv_from_down)memcpy(gfs[nb].recv_from_down,&gf->u_curr[ny*gf->stride+h],(size_t)nx*sizeof(double));
    nb=mythread_decomp_neighbor(dc,g,MYTHREAD_NEIGHBOR_DOWN);
    if(nb>=0&&gfs[nb].recv_from_up)memcpy(gfs[nb].recv_from_up,&gf->u_curr[h*gf->stride+h],(size_t)nx*sizeof(double));
    nb=mythread_decomp_neighbor(dc,g,MYTHREAD_NEIGHBOR_LEFT);
    if(nb>=0&&gfs[nb].recv_from_right)for(int r=0;r<ny;r++)gfs[nb].recv_from_right[r]=gf->u_curr[(h+r)*gf->stride+h];
    nb=mythread_decomp_neighbor(dc,g,MYTHREAD_NEIGHBOR_RIGHT);
    if(nb>=0&&gfs[nb].recv_from_left){int sc=h+nx-1;for(int r=0;r<ny;r++)gfs[nb].recv_from_left[r]=gf->u_curr[(h+r)*gf->stride+sc];}
  }
  for(int g=0;g<ng;g++){GroupField*gf=&gfs[g];if(!gf->u_curr)continue;
    const mythread_tile*t=&dc->group_tiles[g];int nx=t->nx,ny=t->ny;
    if(gf->recv_from_up)memcpy(&gf->u_curr[(ny+h)*gf->stride+h],gf->recv_from_up,(size_t)nx*sizeof(double));
    if(gf->recv_from_down)memcpy(&gf->u_curr[0*gf->stride+h],gf->recv_from_down,(size_t)nx*sizeof(double));
    if(gf->recv_from_left)for(int r=0;r<ny;r++)gf->u_curr[(h+r)*gf->stride+0]=gf->recv_from_left[r];
    if(gf->recv_from_right){int d2=h+nx;for(int r=0;r<ny;r++)gf->u_curr[(h+r)*gf->stride+d2]=gf->recv_from_right[r];}
  }
}

static void halo_mpi(GroupField*gfs,const mythread_decomp*dc){
  if(!gfs||!dc)return;int h=dc->halo,tag=g_mpi_ctx.mpi_tag_base;MPI_Request reqs[8];int nr=0;
  for(int g=0;g<dc->n_groups;g++){GroupField*gf=&gfs[g];if(!gf->u_curr)continue;
    const mythread_tile*t=&dc->group_tiles[g];int nx=t->nx,ny=t->ny,s=gf->stride;
    if(g_sim.neighbor_up>=0&&gf->send_up&&gf->recv_up&&mythread_decomp_is_domain_boundary(dc,g,MYTHREAD_NEIGHBOR_UP)){
      memcpy(gf->send_up,&gf->u_curr[ny*s+h],(size_t)nx*sizeof(double));
      MPI_Isend(gf->send_up,nx,MPI_DOUBLE,g_sim.neighbor_up,tag,MPI_COMM_WORLD,&reqs[nr++]);
      MPI_Irecv(gf->recv_up,nx,MPI_DOUBLE,g_sim.neighbor_up,tag+1,MPI_COMM_WORLD,&reqs[nr++]);}
    if(g_sim.neighbor_down>=0&&gf->send_down&&gf->recv_down&&mythread_decomp_is_domain_boundary(dc,g,MYTHREAD_NEIGHBOR_DOWN)){
      memcpy(gf->send_down,&gf->u_curr[h*s+h],(size_t)nx*sizeof(double));
      MPI_Isend(gf->send_down,nx,MPI_DOUBLE,g_sim.neighbor_down,tag+1,MPI_COMM_WORLD,&reqs[nr++]);
      MPI_Irecv(gf->recv_down,nx,MPI_DOUBLE,g_sim.neighbor_down,tag,MPI_COMM_WORLD,&reqs[nr++]);}
  }
  MPI_Waitall(nr,reqs,MPI_STATUSES_IGNORE);
  for(int g=0;g<dc->n_groups;g++){GroupField*gf=&gfs[g];if(!gf->u_curr)continue;
    const mythread_tile*t=&dc->group_tiles[g];int nx=t->nx,ny=t->ny,s=gf->stride;
    if(gf->recv_up)memcpy(&gf->u_curr[(ny+h)*s+h],gf->recv_up,(size_t)nx*sizeof(double));
    if(gf->recv_down)memcpy(&gf->u_curr[0*s+h],gf->recv_down,(size_t)nx*sizeof(double));}
}

static void apply_d(GroupField*gf,int gid){
  const mythread_tile*t=&g_decomp->group_tiles[gid];int s=gf->stride,h=gf->ny_padded,nx=t->nx,ny=t->ny;
  if(g_sim.local_x_begin==0)for(int y=0;y<h;y++)gf->u_curr[GFIDX(gf,y,HALO)]=0.0;
  if(g_sim.local_x_end==NX){int xr=HALO+nx-1;for(int y=0;y<h;y++)gf->u_curr[GFIDX(gf,y,xr)]=0.0;}
  if(g_sim.neighbor_left<0&&mythread_decomp_is_domain_boundary(g_decomp,gid,MYTHREAD_NEIGHBOR_LEFT))
    for(int y=0;y<h;y++)gf->u_curr[GFIDX(gf,y,0)]=0.0;
  if(g_sim.neighbor_right<0&&mythread_decomp_is_domain_boundary(g_decomp,gid,MYTHREAD_NEIGHBOR_RIGHT))
  {int xh=HALO+nx;for(int y=0;y<h;y++)gf->u_curr[GFIDX(gf,y,xh)]=0.0;}
  if(g_sim.neighbor_down<0&&mythread_decomp_is_domain_boundary(g_decomp,gid,MYTHREAD_NEIGHBOR_DOWN)){
    memset(&gf->u_curr[0],0,(size_t)s*sizeof(double));
    if(g_sim.local_y_begin==0)memset(&gf->u_curr[GFIDX(gf,HALO,0)],0,(size_t)s*sizeof(double));}
  if(g_sim.neighbor_up<0&&mythread_decomp_is_domain_boundary(g_decomp,gid,MYTHREAD_NEIGHBOR_UP)){
    memset(&gf->u_curr[(ny+HALO)*s],0,(size_t)s*sizeof(double));
    if(g_sim.local_y_end==NY)memset(&gf->u_curr[GFIDX(gf,ny,0)],0,(size_t)s*sizeof(double));}
}
static void apply_d_all(void){for(int g=0;g<g_decomp->n_groups;g++)apply_d(&g_gfields[g],g);}
static void copy_prev(void){for(int g=0;g<g_decomp->n_groups;g++)memcpy(g_gfields[g].u_prev,g_gfields[g].u_curr,g_gfields[g].plane_bytes);}

static void comp_region(GroupField*gf,int gid,int yb,int ye,int xb,int xe){
  for(int y=yb;y<ye;y++){int gy=gy_loc(gid,y);if(gy==0||gy==NY-1)continue;
    for(int x=xb;x<xe;x++){int gx=gx_loc(gid,x);if(gx==0||gx==NX-1)continue;
      double u=gf->u_curr[GFIDX(gf,y,x)];
      double d2x=(gf->u_curr[GFIDX(gf,y,x-1)]-2*u+gf->u_curr[GFIDX(gf,y,x+1)])/(DX*DX);
      double d2y=(gf->u_curr[GFIDX(gf,y-1,x)]-2*u+gf->u_curr[GFIDX(gf,y+1,x)])/(DY*DY);
      gf->u_next[GFIDX(gf,y,x)]=2*u-gf->u_prev[GFIDX(gf,y,x)]+C0*C0*DT2*(d2x+d2y);}}}
static void comp_interior(GroupField*gf,int gid,int yb,int ye,int xb,int xe){
  const mythread_tile*t=&g_decomp->group_tiles[gid];int ny=t->ny,nx=t->nx;
  if(g_sim.neighbor_down>=0&&mythread_decomp_is_domain_boundary(g_decomp,gid,MYTHREAD_NEIGHBOR_DOWN)&&yb<HALO+1)yb=HALO+1;
  if(g_sim.neighbor_up>=0&&mythread_decomp_is_domain_boundary(g_decomp,gid,MYTHREAD_NEIGHBOR_UP)&&ye>ny)ye=ny;
  if(g_sim.neighbor_left>=0&&mythread_decomp_is_domain_boundary(g_decomp,gid,MYTHREAD_NEIGHBOR_LEFT)&&xb<HALO+1)xb=HALO+1;
  if(g_sim.neighbor_right>=0&&mythread_decomp_is_domain_boundary(g_decomp,gid,MYTHREAD_NEIGHBOR_RIGHT)&&xe>nx)xe=nx;
  if(yb<ye&&xb<xe)comp_region(gf,gid,yb,ye,xb,xe);}
static void comp_boundary(GroupField*gf,int gid,int yb,int ye,int xb,int xe){
  const mythread_tile*t=&g_decomp->group_tiles[gid];int ny=t->ny,nx=t->nx,lr=HALO,ur=ny,lc=HALO,rc=nx;
  if(xb<HALO)xb=HALO;if(xe>HALO+nx)xe=HALO+nx;if(yb<HALO)yb=HALO;if(ye>HALO+ny)ye=HALO+ny;
  int nd=g_sim.neighbor_down>=0&&mythread_decomp_is_domain_boundary(g_decomp,gid,MYTHREAD_NEIGHBOR_DOWN);
  int nu=g_sim.neighbor_up>=0&&mythread_decomp_is_domain_boundary(g_decomp,gid,MYTHREAD_NEIGHBOR_UP);
  int nl=g_sim.neighbor_left>=0&&mythread_decomp_is_domain_boundary(g_decomp,gid,MYTHREAD_NEIGHBOR_LEFT);
  int nr=g_sim.neighbor_right>=0&&mythread_decomp_is_domain_boundary(g_decomp,gid,MYTHREAD_NEIGHBOR_RIGHT);
  if(nd&&yb<=lr&&lr<ye)comp_region(gf,gid,lr,lr+1,xb,xe);
  if(nu&&ur!=lr&&yb<=ur&&ur<ye)comp_region(gf,gid,ur,ur+1,xb,xe);
  if(nl&&xb<=lc&&lc<xe)comp_region(gf,gid,yb,ye,lc,lc+1);
  if(nr&&rc!=lc&&xb<=rc&&rc<xe)comp_region(gf,gid,yb,ye,rc,rc+1);}

static double init_val(int gy,int gx){
  double cx=0.5*(cfg_USE_FIXED_DOMAIN?1.0:((cfg_NX-1)*cfg_DX)),
         cy=0.5*(cfg_USE_FIXED_DOMAIN?1.0:((cfg_NY-1)*cfg_DY)),
         sigma=0.06*((cfg_USE_FIXED_DOMAIN?1.0:(cfg_DX*(cfg_NX-1)))<(cfg_USE_FIXED_DOMAIN?1.0:(cfg_DY*(cfg_NY-1)))?
                     (cfg_USE_FIXED_DOMAIN?1.0:(cfg_DX*(cfg_NX-1))):(cfg_USE_FIXED_DOMAIN?1.0:(cfg_DY*(cfg_NY-1))));
  return AA*exp(-((gx*DX-cx)*(gx*DX-cx)+(gy*DY-cy)*(gy*DY-cy))/(2.*sigma*sigma));
}

static void run_simulation(void){
  if(gf_alloc(&g_gfields[0],0,g_decomp)!=0)
    {fprintf(stderr,"[Error] gf_alloc\n");MPI_Abort(MPI_COMM_WORLD,1);}

  double *ebuf=(double*)calloc((size_t)g_nthreads,sizeof(double));
#pragma omp parallel num_threads(g_nthreads)
  { int tid=omp_get_thread_num();
    GroupField*gf=&g_gfields[0];
    const mythread_tile*tile=&g_decomp->group_tiles[0];
    int ny=tile->ny,nx=tile->nx;
    int base=ny/g_nthreads,rem=ny%g_nthreads;
    int y0=HALO+tid*base+(tid<rem?tid:rem),y1=y0+base+(tid<rem?1:0);

#pragma omp barrier
    for(int y=y0;y<y1;y++){int gy=gy_loc(0,y);
      for(int x=HALO;x<HALO+nx;x++){int gx=gx_loc(0,x);size_t p=GFIDX(gf,y,x);
        double v=(gx!=0&&gx!=NX-1&&gy!=0&&gy!=NY-1)?init_val(gy,gx):0.0;
        gf->u_curr[p]=v;gf->u_prev[p]=v;gf->u_next[p]=0.0;}}

#pragma omp single
    {apply_d_all();halo_mpi(g_gfields,g_decomp);apply_d_all();copy_prev();}

    double t0=MPI_Wtime();
    for(int step=0;step<NT;step++){MPI_Barrier(MPI_COMM_WORLD);
#pragma omp single
      {apply_d_all();halo_mpi(g_gfields,g_decomp);apply_d_all();}
      for(int y=y0;y<y1;y++)comp_interior(gf,0,y,y+1,HALO,HALO+nx);
#pragma omp barrier
      for(int y=y0;y<y1;y++)comp_boundary(gf,0,y,y+1,HALO,HALO+nx);
#pragma omp barrier
#pragma omp single
      gf_swap(gf);
      /* single has implicit barrier */
      if(tid==0&&g_sim.mpi_rank==0&&(step+1)%((NT>10?NT/10:1))==0)printf("[Main] Step %d/%d\n",step+1,NT);
    }
    /* 计算最终能量（所有线程，验证守恒） */
    { double my_e=0.0;
    for(int y=y0;y<y1;y++){int gy=gy_loc(0,y);if(gy==0||gy==NY-1)continue;
      for(int x=HALO;x<HALO+nx;x++){int gx=gx_loc(0,x);if(gx==0||gx==NX-1)continue;
        double u=gf->u_curr[GFIDX(gf,y,x)],up=gf->u_prev[GFIDX(gf,y,x)];
        double ut=(u-up)/DT,dxc=(gf->u_curr[GFIDX(gf,y,x+1)]-gf->u_curr[GFIDX(gf,y,x-1)])/(2*DX);
        double dyc=(gf->u_curr[GFIDX(gf,y+1,x)]-gf->u_curr[GFIDX(gf,y-1,x)])/(2*DY);
        double dxp=(gf->u_prev[GFIDX(gf,y,x+1)]-gf->u_prev[GFIDX(gf,y,x-1)])/(2*DX);
        double dyp=(gf->u_prev[GFIDX(gf,y+1,x)]-gf->u_prev[GFIDX(gf,y-1,x)])/(2*DY);
        my_e+=0.5*(ut*ut+C0*C0*(dxc*dxp+dyc*dyp))*DX*DY;
      }
    }
    ebuf[tid]=my_e; }
#pragma omp barrier
#pragma omp single
    { double local_e=0;for(int t=0;t<g_nthreads;t++)local_e+=ebuf[t];
      double global_e;MPI_Allreduce(&local_e,&global_e,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      if(g_sim.mpi_rank==0)printf("[Main] Final energy: %.6f (expected %.6f) %s\n",
        global_e,1.570796,fabs(global_e-1.570796)<0.01?"OK":"?"); }

    if(tid==0&&g_sim.mpi_rank==0)printf("[Main] Done: %.3f s, %.1f Mpts/s\n",
      MPI_Wtime()-t0,(double)NX*NY*NT/(MPI_Wtime()-t0)/1e6);
  }
  free(ebuf);
}

static void setup_domain(int mr,int ms){
  mythread_decomp*tmp=mythread_decomp_create(NX,NY,HALO,ms,0,MYTHREAD_DECOMP_XY_2D);
  int px=tmp->gx,py=tmp->gy;mythread_decomp_free(tmp);
  g_sim.proc_px=px;g_sim.proc_py=py;g_sim.proc_x=mr%px;g_sim.proc_y=mr/px;
  {int tn=NX/px,rn=NX%px,tm=NY/py,rm=NY%py;
   int ox=g_sim.proc_x*tn+(g_sim.proc_x<rn?g_sim.proc_x:rn),oy=g_sim.proc_y*tm+(g_sim.proc_y<rm?g_sim.proc_y:rm);
   g_sim.local_x_begin=ox;g_sim.local_x_end=ox+tn+(g_sim.proc_x<rn?1:0);
   g_sim.local_y_begin=oy;g_sim.local_y_end=oy+tm+(g_sim.proc_y<rm?1:0);}
  g_sim.local_nx=g_sim.local_x_end-g_sim.local_x_begin;g_sim.local_ny=g_sim.local_y_end-g_sim.local_y_begin;
  g_sim.neighbor_left=(g_sim.proc_x>0)?(mr-1):-1;g_sim.neighbor_right=(g_sim.proc_x+1<px)?(mr+1):-1;
  g_sim.neighbor_down=(g_sim.proc_y>0)?(mr-px):-1;g_sim.neighbor_up=(g_sim.proc_y+1<py)?(mr+px):-1;
  g_mpi_ctx.mpirank_up=g_sim.neighbor_up;g_mpi_ctx.mpirank_down=g_sim.neighbor_down;
  g_mpi_ctx.mpirank_left=g_sim.neighbor_left;g_mpi_ctx.mpirank_right=g_sim.neighbor_right;g_mpi_ctx.mpi_tag_base=100;
}

static void cfg_load(const char*cp,const char*hp){
  mythread_cfg*hw=mythread_cfg_load(hp),*cs=mythread_cfg_load(cp);
  cfg_NX=cs?mythread_cfg_get_int(cs,"","NX",cfg_NX):cfg_NX;cfg_NY=cs?mythread_cfg_get_int(cs,"","NY",cfg_NY):cfg_NY;
  cfg_NT=cs?mythread_cfg_get_int(cs,"","NT",cfg_NT):cfg_NT;cfg_A=cs?mythread_cfg_get_double(cs,"","A",cfg_A):cfg_A;
  cfg_DT=cs?mythread_cfg_get_double(cs,"","DT",cfg_DT):cfg_DT;cfg_C0=cs?mythread_cfg_get_double(cs,"","C0",cfg_C0):cfg_C0;
  cfg_HALO=cs?mythread_cfg_get_int(cs,"","HALO",cfg_HALO):cfg_HALO;
  cfg_N_GROUPS=hw?mythread_cfg_get_int(hw,"","N_GROUPS",cfg_N_GROUPS):cfg_N_GROUPS;
  cfg_N_WORKERS=hw?mythread_cfg_get_int(hw,"","N_WORKERS",cfg_N_WORKERS):cfg_N_WORKERS;
  cfg_GROUP_DECOMP=hw?mythread_cfg_get_int(hw,"","GROUP_DECOMP",cfg_GROUP_DECOMP):cfg_GROUP_DECOMP;
  mythread_cfg_free(cs);mythread_cfg_free(hw);cfg_compute_derived();
}

int main(int argc,char**argv){
  int mr,ms,ns,lr=1,gr=1,prov;MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE,&prov);
  MPI_Comm_rank(MPI_COMM_WORLD,&mr);MPI_Comm_size(MPI_COMM_WORLD,&ms);
  MPI_Comm nc;MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,0,MPI_INFO_NULL,&nc);MPI_Comm_size(nc,&ns);
  {const char*cp=mythread_env_get("WAVE_CASE_CFG","config/case.cfg"),
        *hp=mythread_env_get("WAVE_HARDWARE_CFG","config/hardware.cfg");cfg_load(cp,hp);}
  g_nthreads=cfg_N_WORKERS>0?cfg_N_WORKERS:omp_get_max_threads();
  omp_set_num_threads(g_nthreads);
  if(NX<3||NY<3){if(mr==0)fprintf(stderr,"[Error] NX/NY>=3\n");MPI_Finalize();return 1;}
  if(CFL_SUM2>1.0){if(mr==0)fprintf(stderr,"[Error] CFL>1\n");MPI_Finalize();return 1;}
  memset(&g_sim,0,sizeof(g_sim));g_sim.mpi_rank=mr;g_sim.mpi_size=ms;setup_domain(mr,ms);
  if(g_sim.local_nx<=0||g_sim.local_ny<=0)lr=0;
  MPI_Allreduce(&lr,&gr,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);if(!gr){MPI_Finalize();return 1;}
  g_decomp=mythread_decomp_create(g_sim.local_nx,g_sim.local_ny,HALO,1,g_nthreads,cfg_GROUP_DECOMP);
  if(!g_decomp){fprintf(stderr,"[Error] decomp\n");MPI_Finalize();return 1;}
  g_gfields=(GroupField*)xc(1,sizeof(GroupField));
  if(mr==0)printf("=== Wave OMP %dx%dx%d MPI=%d threads=%d DT=%.4f CFL=%.4f ===\n",
    NX,NY,NT,ms,g_nthreads,cfg_DT,cfg_C0*cfg_DT/cfg_DX);
  MPI_Barrier(MPI_COMM_WORLD);
  double t0=MPI_Wtime();run_simulation();
  if(mr==0)printf("[Main] Wall: %.3f s\n",MPI_Wtime()-t0);
  gf_free(&g_gfields[0]);free(g_gfields);
  mythread_decomp_free((mythread_decomp*)g_decomp);MPI_Finalize();return 0;
}
