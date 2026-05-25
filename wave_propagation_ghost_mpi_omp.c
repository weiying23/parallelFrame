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
static int    cfg_HALO=1,cfg_N_THREADS=0;

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
  double *u_prev,*u_curr,*u_next;
  int ny_padded,nx_padded,stride; size_t plane_bytes;
  double *send_up,*recv_up,*send_down,*recv_down;
  double *send_left,*recv_left,*send_right,*recv_right;
} Field;

#define F(f,y,x) ((size_t)(y)*(size_t)(f)->stride+(size_t)(x))

typedef struct {
  int mx,my,nx,ny,mr,ms,px,py,px_id,py_id;
  int nb_l,nb_r,nb_d,nb_u;
} Dom;

static Dom g_d={0};
static Field g_f={0};
static const mythread_decomp *g_dc=NULL;

static inline double wtime(void){return MPI_Wtime();}
static inline int gx(int lx){return g_d.mx+(g_dc->group_tiles[0].x_begin+lx-HALO);}
static inline int gy(int ly){return g_d.my+(g_dc->group_tiles[0].y_begin+ly-HALO);}

static void*xc(size_t n,size_t s){void*p=calloc(n,s);if(!p){fprintf(stderr,"[E] alloc\n");MPI_Abort(MPI_COMM_WORLD,1);}return p;}

static int f_alloc(Field*f){
  const mythread_tile*t=&g_dc->group_tiles[0];int h=HALO,nx=t->nx,ny=t->ny;
  memset(f,0,sizeof(*f));f->ny_padded=ny+2*h;f->nx_padded=nx+2*h;f->stride=f->nx_padded;
  f->plane_bytes=(size_t)f->ny_padded*(size_t)f->nx_padded*sizeof(double);
#define AL(p,sz) do{(p)=(double*)calloc((sz),sizeof(double));if(!(p))goto fail;}while(0)
  AL(f->u_prev,(size_t)f->ny_padded*f->nx_padded);
  AL(f->u_curr,(size_t)f->ny_padded*f->nx_padded);
  AL(f->u_next,(size_t)f->ny_padded*f->nx_padded);
  if(g_d.nb_u>=0||g_dc->n_groups>1&&mythread_decomp_is_domain_boundary(g_dc,0,MYTHREAD_NEIGHBOR_UP)){AL(f->send_up,nx);AL(f->recv_up,nx);}
  if(g_d.nb_d>=0||g_dc->n_groups>1&&mythread_decomp_is_domain_boundary(g_dc,0,MYTHREAD_NEIGHBOR_DOWN)){AL(f->send_down,nx);AL(f->recv_down,nx);}
  if(g_d.nb_l>=0||g_dc->n_groups>1&&mythread_decomp_is_domain_boundary(g_dc,0,MYTHREAD_NEIGHBOR_LEFT)){AL(f->send_left,ny);AL(f->recv_left,ny);}
  if(g_d.nb_r>=0||g_dc->n_groups>1&&mythread_decomp_is_domain_boundary(g_dc,0,MYTHREAD_NEIGHBOR_RIGHT)){AL(f->send_right,ny);AL(f->recv_right,ny);}
#undef AL
  return 0;
fail:
  free(f->u_prev);free(f->u_curr);free(f->u_next);
  free(f->send_up);free(f->recv_up);free(f->send_down);free(f->recv_down);
  free(f->send_left);free(f->recv_left);free(f->send_right);free(f->recv_right);
  memset(f,0,sizeof(*f));return-1;
}
static void f_free(Field*f){if(!f)return;
  free(f->u_prev);free(f->u_curr);free(f->u_next);
  free(f->send_up);free(f->recv_up);free(f->send_down);free(f->recv_down);
  free(f->send_left);free(f->recv_left);free(f->send_right);free(f->recv_right);
  memset(f,0,sizeof(*f));}
static void f_swap(Field*f){if(!f)return;
  double*t=f->u_prev;f->u_prev=f->u_curr;f->u_curr=f->u_next;f->u_next=t;}

static void apply_d(Field*f){
  const mythread_tile*t=&g_dc->group_tiles[0];int s=f->stride,h=f->ny_padded,nx=t->nx,ny=t->ny;
  if(g_d.mx==0)for(int y=0;y<h;y++)f->u_curr[F(f,y,HALO)]=0.0;
  if(g_d.mx+g_d.nx==NX){int xr=HALO+nx-1;for(int y=0;y<h;y++)f->u_curr[F(f,y,xr)]=0.0;}
  if(g_d.nb_l<0)for(int y=0;y<h;y++)f->u_curr[F(f,y,0)]=0.0;
  if(g_d.nb_r<0){int xh=HALO+nx;for(int y=0;y<h;y++)f->u_curr[F(f,y,xh)]=0.0;}
  if(g_d.nb_d<0){memset(&f->u_curr[0],0,(size_t)s*sizeof(double));
    if(g_d.my==0)memset(&f->u_curr[F(f,HALO,0)],0,(size_t)s*sizeof(double));}
  if(g_d.nb_u<0){memset(&f->u_curr[(ny+HALO)*s],0,(size_t)s*sizeof(double));
    if(g_d.my+g_d.ny==NY)memset(&f->u_curr[F(f,ny,0)],0,(size_t)s*sizeof(double));}
}

static void halo_x(Field*f){
  int nx=g_dc->group_tiles[0].nx,ny=g_dc->group_tiles[0].ny,s=f->stride,tag=100,tag_lr=200;
  MPI_Request r[8];int nr=0;
  if(g_d.nb_u>=0&&f->send_up){
    memcpy(f->send_up,&f->u_curr[ny*s+HALO],(size_t)nx*sizeof(double));
    MPI_Isend(f->send_up,nx,MPI_DOUBLE,g_d.nb_u,tag,MPI_COMM_WORLD,&r[nr++]);
    MPI_Irecv(f->recv_up,nx,MPI_DOUBLE,g_d.nb_u,tag+1,MPI_COMM_WORLD,&r[nr++]);}
  if(g_d.nb_d>=0&&f->send_down){
    memcpy(f->send_down,&f->u_curr[HALO*s+HALO],(size_t)nx*sizeof(double));
    MPI_Isend(f->send_down,nx,MPI_DOUBLE,g_d.nb_d,tag+1,MPI_COMM_WORLD,&r[nr++]);
    MPI_Irecv(f->recv_down,nx,MPI_DOUBLE,g_d.nb_d,tag,MPI_COMM_WORLD,&r[nr++]);}
  if(g_d.nb_r>=0&&f->send_right){
    int xr=HALO+nx-1;
    for(int y=0;y<ny;y++)f->send_right[y]=f->u_curr[(HALO+y)*s+xr];
    MPI_Isend(f->send_right,ny,MPI_DOUBLE,g_d.nb_r,tag_lr,MPI_COMM_WORLD,&r[nr++]);
    MPI_Irecv(f->recv_right,ny,MPI_DOUBLE,g_d.nb_r,tag_lr+1,MPI_COMM_WORLD,&r[nr++]);}
  if(g_d.nb_l>=0&&f->send_left){
    for(int y=0;y<ny;y++)f->send_left[y]=f->u_curr[(HALO+y)*s+HALO];
    MPI_Isend(f->send_left,ny,MPI_DOUBLE,g_d.nb_l,tag_lr+1,MPI_COMM_WORLD,&r[nr++]);
    MPI_Irecv(f->recv_left,ny,MPI_DOUBLE,g_d.nb_l,tag_lr,MPI_COMM_WORLD,&r[nr++]);}
  MPI_Waitall(nr,r,MPI_STATUSES_IGNORE);
  if(f->recv_up)memcpy(&f->u_curr[(ny+HALO)*s+HALO],f->recv_up,(size_t)nx*sizeof(double));
  if(f->recv_down)memcpy(&f->u_curr[0*s+HALO],f->recv_down,(size_t)nx*sizeof(double));
  if(f->recv_right){int xh=HALO+nx;for(int y=0;y<ny;y++)f->u_curr[(HALO+y)*s+xh]=f->recv_right[y];}
  if(f->recv_left)for(int y=0;y<ny;y++)f->u_curr[(HALO+y)*s+0]=f->recv_left[y];
}

static double init_val(int gy,int gx){
  double cx=0.5*(cfg_USE_FIXED_DOMAIN?1.0:((cfg_NX-1)*cfg_DX)),
         cy=0.5*(cfg_USE_FIXED_DOMAIN?1.0:((cfg_NY-1)*cfg_DY)),
         sigma=0.06*((cfg_USE_FIXED_DOMAIN?1.0:(cfg_DX*(cfg_NX-1)))<(cfg_USE_FIXED_DOMAIN?1.0:(cfg_DY*(cfg_NY-1)))?
           (cfg_USE_FIXED_DOMAIN?1.0:(cfg_DX*(cfg_NX-1))):(cfg_USE_FIXED_DOMAIN?1.0:(cfg_DY*(cfg_NY-1))));
  return AA*exp(-((gx*DX-cx)*(gx*DX-cx)+(gy*DY-cy)*(gy*DY-cy))/(2.*sigma*sigma));
}

static void comp_interior(Field*f){
  const mythread_tile*t=&g_dc->group_tiles[0];int ny=t->ny,nx=t->nx;
  int yb=HALO,ye=HALO+ny,xb=HALO,xe=HALO+nx;
  const int gx_base=g_d.mx+t->x_begin-HALO,gy_base=g_d.my+t->y_begin-HALO;
  const double inv_dx2=1.0/(DX*DX),inv_dy2=1.0/(DY*DY),c2dt2=C0*C0*DT2;
  if(g_d.nb_d>=0)yb=HALO+1;
  if(g_d.nb_u>=0)ye=ny;
  #pragma omp parallel for schedule(static)
  for(int y=yb;y<ye;y++){int gy_=gy_base+y;if(gy_==0||gy_==NY-1)continue;
    for(int x=xb;x<xe;x++){int gx_=gx_base+x;if(gx_==0||gx_==NX-1)continue;
      double u=f->u_curr[F(f,y,x)];
      double d2x=(f->u_curr[F(f,y,x-1)]-2*u+f->u_curr[F(f,y,x+1)])*inv_dx2;
      double d2y=(f->u_curr[F(f,y-1,x)]-2*u+f->u_curr[F(f,y+1,x)])*inv_dy2;
      f->u_next[F(f,y,x)]=2*u-f->u_prev[F(f,y,x)]+c2dt2*(d2x+d2y);}}
}

static void comp_boundary(Field*f){
  const mythread_tile*t=&g_dc->group_tiles[0];int ny=t->ny,nx=t->nx;
  int lr=HALO,ur=ny;
  const int gx_base=g_d.mx+t->x_begin-HALO;
  const double inv_dx2=1.0/(DX*DX),inv_dy2=1.0/(DY*DY),c2dt2=C0*C0*DT2;
  if(g_d.nb_d>=0){
    #pragma omp parallel for schedule(static)
    for(int x=HALO;x<HALO+nx;x++){
      int gx_=gx_base+x;if(gx_==0||gx_==NX-1)continue;
      double u=f->u_curr[F(f,lr,x)];
      double d2x=(f->u_curr[F(f,lr,x-1)]-2*u+f->u_curr[F(f,lr,x+1)])*inv_dx2;
      double d2y=(f->u_curr[F(f,lr-1,x)]-2*u+f->u_curr[F(f,lr+1,x)])*inv_dy2;
      f->u_next[F(f,lr,x)]=2*u-f->u_prev[F(f,lr,x)]+c2dt2*(d2x+d2y);}
  }
  if(g_d.nb_u>=0){
    #pragma omp parallel for schedule(static)
    for(int x=HALO;x<HALO+nx;x++){
      int gx_=gx_base+x;if(gx_==0||gx_==NX-1)continue;
      double u=f->u_curr[F(f,ur,x)];
      double d2x=(f->u_curr[F(f,ur,x-1)]-2*u+f->u_curr[F(f,ur,x+1)])*inv_dx2;
      double d2y=(f->u_curr[F(f,ur-1,x)]-2*u+f->u_curr[F(f,ur+1,x)])*inv_dy2;
      f->u_next[F(f,ur,x)]=2*u-f->u_prev[F(f,ur,x)]+c2dt2*(d2x+d2y);}
  }
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
}

static void cfg_load(const char*cp,const char*hp){
  mythread_cfg*hw=mythread_cfg_load(hp),*cs=mythread_cfg_load(cp);
  cfg_NX=cs?mythread_cfg_get_int(cs,"","NX",cfg_NX):cfg_NX;cfg_NY=cs?mythread_cfg_get_int(cs,"","NY",cfg_NY):cfg_NY;
  cfg_NT=cs?mythread_cfg_get_int(cs,"","NT",cfg_NT):cfg_NT;cfg_A=cs?mythread_cfg_get_double(cs,"","A",cfg_A):cfg_A;
  cfg_DT=cs?mythread_cfg_get_double(cs,"","DT",cfg_DT):cfg_DT;cfg_C0=cs?mythread_cfg_get_double(cs,"","C0",cfg_C0):cfg_C0;
  cfg_HALO=cs?mythread_cfg_get_int(cs,"","HALO",cfg_HALO):cfg_HALO;
  cfg_N_THREADS=hw?mythread_cfg_get_int(hw,"","N_WORKERS",cfg_N_THREADS):cfg_N_THREADS;
  mythread_cfg_free(cs);mythread_cfg_free(hw);cfg_compute_derived();
}

int main(int argc,char**argv){
  int mr,ms,lr=1,gr=1;
  MPI_Init(&argc,&argv);MPI_Comm_rank(MPI_COMM_WORLD,&mr);MPI_Comm_size(MPI_COMM_WORLD,&ms);
  {const char*cp=mythread_env_get("WAVE_CASE_CFG","config/case.cfg"),
        *hp=mythread_env_get("WAVE_HARDWARE_CFG","config/hardware.cfg");cfg_load(cp,hp);}
  if(cfg_N_THREADS>0)omp_set_num_threads(cfg_N_THREADS);
  int nth=omp_get_max_threads();
  if(NX<3||NY<3){if(mr==0)fprintf(stderr,"[E] NX/NY>=3\n");MPI_Finalize();return 1;}
  if(CFL_SUM2>1.0){if(mr==0)fprintf(stderr,"[E] CFL>1\n");MPI_Finalize();return 1;}

  setup_domain(mr,ms);
  if(g_d.nx<=0||g_d.ny<=0)lr=0;
  MPI_Allreduce(&lr,&gr,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);if(!gr){MPI_Finalize();return 1;}

  g_dc=mythread_decomp_create(g_d.nx,g_d.ny,HALO,1,0,0);
  if(!g_dc||f_alloc(&g_f)!=0){fprintf(stderr,"[E] init\n");MPI_Finalize();return 1;}

  if(mr==0)printf("=== Wave MPI+OMP %dx%dx%d procs=%d (%dx%d) threads/rank=%d DT=%.4f CFL=%.4f ===\n",
    NX,NY,NT,ms,g_d.px,g_d.py,nth,cfg_DT,cfg_C0*cfg_DT/cfg_DX);
  printf("[R%d] domain %dx%d at (%d,%d) nbr L=%d R=%d D=%d U=%d\n",
    mr,g_d.nx,g_d.ny,g_d.mx,g_d.my,g_d.nb_l,g_d.nb_r,g_d.nb_d,g_d.nb_u);

  /* Init field */
  Field*f=&g_f;const mythread_tile*t=&g_dc->group_tiles[0];
  #pragma omp parallel for schedule(static)
  for(int y=HALO;y<HALO+t->ny;y++){int gy_=gy(y);
    for(int x=HALO;x<HALO+t->nx;x++){int gx_=gx(x);
      double v=(gx_!=0&&gx_!=NX-1&&gy_!=0&&gy_!=NY-1)?init_val(gy_,gx_):0.0;
      f->u_curr[F(f,y,x)]=v;f->u_prev[F(f,y,x)]=v;f->u_next[F(f,y,x)]=0.0;}}

  apply_d(f);halo_x(f);apply_d(f);
  memcpy(f->u_prev,f->u_curr,f->plane_bytes);

  double t0=MPI_Wtime();
  for(int step=0;step<NT;step++){
    MPI_Barrier(MPI_COMM_WORLD);
    apply_d(f);halo_x(f);apply_d(f);
    comp_interior(f);
    comp_boundary(f);
    f_swap(f);
    if(mr==0&&(step+1)%((NT>10?NT/10:1))==0)printf("[Main] Step %d/%d\n",step+1,NT);
  }
  double elapsed=MPI_Wtime()-t0;

  /* Final energy */
  double my_e=0.0;
  #pragma omp parallel for reduction(+:my_e) schedule(static)
  for(int y=HALO;y<HALO+t->ny;y++){int gy_=gy(y);if(gy_==0||gy_==NY-1)continue;
    for(int x=HALO;x<HALO+t->nx;x++){int gx_=gx(x);if(gx_==0||gx_==NX-1)continue;
      double u=f->u_curr[F(f,y,x)],up=f->u_prev[F(f,y,x)];
      double ut=(u-up)/DT,dxc=(f->u_curr[F(f,y,x+1)]-f->u_curr[F(f,y,x-1)])/(2*DX);
      double dyc=(f->u_curr[F(f,y+1,x)]-f->u_curr[F(f,y-1,x)])/(2*DY);
      double dxp=(f->u_prev[F(f,y,x+1)]-f->u_prev[F(f,y,x-1)])/(2*DX);
      double dyp=(f->u_prev[F(f,y+1,x)]-f->u_prev[F(f,y-1,x)])/(2*DY);
      my_e+=0.5*(ut*ut+C0*C0*(dxc*dxp+dyc*dyp))*DX*DY;
    }}
  double global_e;MPI_Allreduce(&my_e,&global_e,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  if(mr==0){printf("[Main] Done: %.3f s, %.1f Mpts/s\n",elapsed,(double)NX*NY*NT/elapsed/1e6);
    printf("[Main] Final energy: %.6f (expected %.6f) %s\n",global_e,1.570796,fabs(global_e-1.570796)<0.01?"OK":"?");}

  f_free(&g_f);mythread_decomp_free((mythread_decomp*)g_dc);MPI_Finalize();return 0;
}
