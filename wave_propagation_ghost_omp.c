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

/* ── 运行时参数 ── */
static int    cfg_NX=14000, cfg_NY=14000, cfg_NT=480;
static double cfg_A=10.0, cfg_DT=0.001, cfg_C0=0.1;
static int    cfg_USE_FIXED_DOMAIN=0;
static double cfg_DX, cfg_DY, cfg_LX, cfg_LY, cfg_DT2, cfg_CFL_X, cfg_CFL_Y, cfg_CFL_SUM2;
static int    cfg_HALO=1, cfg_N_WORKERS=3;
static int    cfg_ENERGY_REPORT_INTERVAL=60, cfg_GROUP_DECOMP=0;
static int    cfg_NCorePClu=5, cfg_NCluPNode=2, cfg_NCorePGrp=4, cfg_ManageCoreId=4;

static void cfg_compute_derived(void) {
  if (cfg_USE_FIXED_DOMAIN) {
    cfg_LX = 1.0; cfg_LY = 1.0;
    cfg_DX = cfg_LX / (double)(cfg_NX - 1);
    cfg_DY = cfg_LY / (double)(cfg_NY - 1);
  } else {
    cfg_DX = 0.01; cfg_DY = 0.01;
    cfg_LX = (cfg_NX - 1) * cfg_DX;
    cfg_LY = (cfg_NY - 1) * cfg_DY;
  }
  cfg_DT2   = cfg_DT * cfg_DT;
  cfg_CFL_X = cfg_C0 * cfg_DT / cfg_DX;
  cfg_CFL_Y = cfg_C0 * cfg_DT / cfg_DY;
  cfg_CFL_SUM2 = cfg_CFL_X * cfg_CFL_X + cfg_CFL_Y * cfg_CFL_Y;
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

/* ── 简化版 GroupField ── */
typedef struct GroupField {
  double *u_prev, *u_curr, *u_next;
  int ny_padded, nx_padded, stride;
  size_t plane_bytes;
  double *send_up,    *recv_up;
  double *send_down,  *recv_down;
  double *send_left,  *recv_left;
  double *send_right, *recv_right;
  double group_energy;
} GroupField;

#define GFIDX(gf, y, x) ((size_t)(y) * (size_t)(gf)->stride + (size_t)(x))

/* ── 简化版 MPI 上下文 ── */
typedef struct {
  int mpirank_up, mpirank_down, mpirank_left, mpirank_right;
  int mpi_tag_base;
} MpiCtx;

/* ── 进程域信息 ── */
typedef struct {
  int local_x_begin, local_x_end, local_nx;
  int local_y_begin, local_y_end, local_ny;
  int mpi_rank, mpi_size;
  int proc_x, proc_y, proc_px, proc_py;
  int neighbor_left, neighbor_right, neighbor_up, neighbor_down;
  double initial_energy;
} SimData;

/* 计时：pad 到 8 doubles（64B cache line）避免 false sharing */
enum { TM_COMP=0, TM_ENERGY=1, TM_HALO=2, TM_N_SLOTS_PAD=8 };
#define TM_SLOT(tid) ((tid) * TM_N_SLOTS_PAD)

static SimData g_sim = {0};
static GroupField *g_gfields = NULL;
static const mythread_decomp *g_decomp = NULL;
static MpiCtx g_mpi_ctx = {0};
static double *g_l2_acc   = NULL;
static double *g_max_acc  = NULL;
static int    *g_max_x    = NULL, *g_max_y = NULL;
static double *g_times    = NULL;  /* [max_threads * TM_N_SLOTS_PAD] — cache-line padded */
static int g_nthreads;

static inline double wall_time(void) { return MPI_Wtime(); }

static inline int global_x_from_local(int gid, int local_x) {
  return g_sim.local_x_begin + (g_decomp->group_tiles[gid].x_begin + local_x - HALO);
}
static inline int global_y_from_local(int gid, int local_y) {
  return g_sim.local_y_begin + (g_decomp->group_tiles[gid].y_begin + local_y - HALO);
}

static void *xcalloc(size_t count, size_t size) {
  void *p = calloc(count, size);
  if (!p) { fprintf(stderr,"[Error] MPI=%d: alloc failed\n",g_sim.mpi_rank); MPI_Abort(MPI_COMM_WORLD,1); }
  return p;
}
static int add_size_checked(size_t a, size_t b, size_t *out) {
  if (a > SIZE_MAX - b) return 0; *out = a + b; return 1;
}
static int multiply_size_checked(size_t a, size_t b, size_t *out) {
  if (a != 0 && b > SIZE_MAX / a) return 0; *out = a * b; return 1;
}

static int gf_alloc(GroupField *gf, int gid, const mythread_decomp *dc) {
  if (!gf || !dc || gid < 0 || gid >= dc->n_groups) return -1;
  const mythread_tile *tile = &dc->group_tiles[gid];
  int halo = dc->halo, nx = tile->nx, ny = tile->ny;
  memset(gf, 0, sizeof(*gf));
  gf->ny_padded = ny + 2*halo; gf->nx_padded = nx + 2*halo; gf->stride = gf->nx_padded;
  gf->plane_bytes = (size_t)gf->ny_padded * (size_t)gf->nx_padded * sizeof(double);
#define ALLOC(p,sz) do{(p)=(double*)calloc((sz),sizeof(double));if(!(p))goto fail;}while(0)
  ALLOC(gf->u_prev,(size_t)gf->ny_padded*gf->nx_padded);
  ALLOC(gf->u_curr,(size_t)gf->ny_padded*gf->nx_padded);
  ALLOC(gf->u_next,(size_t)gf->ny_padded*gf->nx_padded);
  if(mythread_decomp_is_domain_boundary(dc,gid,MYTHREAD_NEIGHBOR_UP)){ALLOC(gf->send_up,nx);ALLOC(gf->recv_up,nx);}
  if(mythread_decomp_is_domain_boundary(dc,gid,MYTHREAD_NEIGHBOR_DOWN)){ALLOC(gf->send_down,nx);ALLOC(gf->recv_down,nx);}
  if(mythread_decomp_is_domain_boundary(dc,gid,MYTHREAD_NEIGHBOR_LEFT)){ALLOC(gf->send_left,ny);ALLOC(gf->recv_left,ny);}
  if(mythread_decomp_is_domain_boundary(dc,gid,MYTHREAD_NEIGHBOR_RIGHT)){ALLOC(gf->send_right,ny);ALLOC(gf->recv_right,ny);}
#undef ALLOC
  return 0;
fail:
  free(gf->u_prev);free(gf->u_curr);free(gf->u_next);
  free(gf->send_up);free(gf->recv_up);free(gf->send_down);free(gf->recv_down);
  free(gf->send_left);free(gf->recv_left);free(gf->send_right);free(gf->recv_right);
  memset(gf,0,sizeof(*gf)); return -1;
}
static void gf_free(GroupField *gf) {
  if(!gf)return;
  free(gf->u_prev);free(gf->u_curr);free(gf->u_next);
  free(gf->send_up);free(gf->recv_up);free(gf->send_down);free(gf->recv_down);
  free(gf->send_left);free(gf->recv_left);free(gf->send_right);free(gf->recv_right);
  memset(gf,0,sizeof(*gf));
}
static void gf_swap(GroupField *gf) {
  if(!gf)return;
  double *tmp=gf->u_prev; gf->u_prev=gf->u_curr; gf->u_curr=gf->u_next; gf->u_next=tmp;
}

static int validate_memory(int mpi_rank, int node_size, const mythread_decomp *dc) {
  size_t total = 0;
  for (int g=0; g<dc->n_groups; g++) {
    const mythread_tile *t = &dc->group_tiles[g];
    size_t pw, ph, plane_elem, planes, halo_y, halo_x, mpi_bufs, gf_total;
    if(!add_size_checked((size_t)t->nx,2u*(size_t)dc->halo,&pw))return 0;
    if(!add_size_checked((size_t)t->ny,2u*(size_t)dc->halo,&ph))return 0;
    if(!multiply_size_checked(ph,pw,&plane_elem))return 0;
    if(!multiply_size_checked(plane_elem,sizeof(double)*3,&planes))return 0;
    if(!multiply_size_checked((size_t)t->nx,sizeof(double)*2,&halo_y))return 0;
    if(!multiply_size_checked((size_t)t->ny,sizeof(double)*2,&halo_x))return 0;
    mpi_bufs=0;
    if(mythread_decomp_is_domain_boundary(dc,g,MYTHREAD_NEIGHBOR_UP))
      {size_t b;multiply_size_checked((size_t)t->nx,sizeof(double)*2,&b);mpi_bufs+=b;}
    if(mythread_decomp_is_domain_boundary(dc,g,MYTHREAD_NEIGHBOR_DOWN))
      {size_t b;multiply_size_checked((size_t)t->nx,sizeof(double)*2,&b);mpi_bufs+=b;}
    if(mythread_decomp_is_domain_boundary(dc,g,MYTHREAD_NEIGHBOR_LEFT))
      {size_t b;multiply_size_checked((size_t)t->ny,sizeof(double)*2,&b);mpi_bufs+=b;}
    if(mythread_decomp_is_domain_boundary(dc,g,MYTHREAD_NEIGHBOR_RIGHT))
      {size_t b;multiply_size_checked((size_t)t->ny,sizeof(double)*2,&b);mpi_bufs+=b;}
    if(!add_size_checked(planes,halo_y,&gf_total))return 0;
    if(!add_size_checked(gf_total,halo_x,&gf_total))return 0;
    if(!add_size_checked(gf_total,mpi_bufs,&gf_total))return 0;
    if(!add_size_checked(total,gf_total,&total))return 0;
  }
  if(mpi_rank==0)printf("Estimated rank memory: %.3f GiB (%d groups)\n",(double)total/(1024.*1024.*1024.),dc->n_groups);
  if(node_size<1)node_size=1;
  size_t nb; if(!multiply_size_checked(total,(size_t)node_size,&nb))return 0;
  if(mpi_rank==0)printf("Estimated node memory: %.3f GiB (%d ranks/node)\n",(double)nb/(1024.*1024.*1024.),node_size);
#if defined(__linux__)
  {struct sysinfo info;
   if(sysinfo(&info)==0){
     unsigned long long vis=(unsigned long long)info.totalram*(unsigned long long)info.mem_unit;
     if(mpi_rank==0)printf("Visible node memory : %.3f GiB\n",(double)vis/(1024.*1024.*1024.));
     if((unsigned long long)nb>vis){if(mpi_rank==0)fprintf(stderr,"[Error] memory exceeds visible\n");return 0;}}}
#endif
  return 1;
}

static void setup_process_domain(int mpi_rank, int mpi_size) {
  mythread_decomp *tmp = mythread_decomp_create(NX,NY,HALO,mpi_size,0,MYTHREAD_DECOMP_XY_2D);
  int px=tmp->gx, py=tmp->gy; mythread_decomp_free(tmp);
  g_sim.proc_px=px; g_sim.proc_py=py;
  g_sim.proc_x=mpi_rank%px; g_sim.proc_y=mpi_rank/px;
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

static double initial_condition_value(int gy, int gx) {
  double cx=0.5*LX,cy=0.5*LY,sigma=0.06*((LX<LY)?LX:LY);
  return A*exp(-((gx*DX-cx)*(gx*DX-cx)+(gy*DY-cy)*(gy*DY-cy))/(2.*sigma*sigma));
}

/* ── Dirichlet 边界（只清零物理边界，不做全平面扫描）── */
static void apply_dirichlet(GroupField *gf, int gid) {
  const mythread_tile *tile = &g_decomp->group_tiles[gid];
  int s=gf->stride, h=gf->ny_padded, nx=tile->nx, ny=tile->ny;
  if(g_sim.local_x_begin==0)
    for(int y=0;y<h;y++) gf->u_curr[GFIDX(gf,y,HALO)]=0.0;
  if(g_sim.local_x_end==NX){
    int xr=HALO+nx-1;
    for(int y=0;y<h;y++) gf->u_curr[GFIDX(gf,y,xr)]=0.0;
  }
  /* MPI halo 归零（仅边界进程 + 域边界组） */
  if(g_sim.neighbor_left<0&&mythread_decomp_is_domain_boundary(g_decomp,gid,MYTHREAD_NEIGHBOR_LEFT))
    for(int y=0;y<h;y++) gf->u_curr[GFIDX(gf,y,0)]=0.0;
  if(g_sim.neighbor_right<0&&mythread_decomp_is_domain_boundary(g_decomp,gid,MYTHREAD_NEIGHBOR_RIGHT)){
    int xh=HALO+nx;
    for(int y=0;y<h;y++) gf->u_curr[GFIDX(gf,y,xh)]=0.0;
  }
  if(g_sim.neighbor_down<0&&mythread_decomp_is_domain_boundary(g_decomp,gid,MYTHREAD_NEIGHBOR_DOWN)){
    memset(&gf->u_curr[0],0,(size_t)s*sizeof(double));
    /* 同时归零物理边界行（原 zero_physical_y 的逻辑） */
    if(g_sim.local_y_begin==0) memset(&gf->u_curr[GFIDX(gf,HALO,0)],0,(size_t)s*sizeof(double));
  }
  if(g_sim.neighbor_up<0&&mythread_decomp_is_domain_boundary(g_decomp,gid,MYTHREAD_NEIGHBOR_UP)){
    memset(&gf->u_curr[(ny+HALO)*s],0,(size_t)s*sizeof(double));
    if(g_sim.local_y_end==NY) memset(&gf->u_curr[GFIDX(gf,ny,0)],0,(size_t)s*sizeof(double));
  }
}

static void apply_dirichlet_all(void) {
  for(int g=0;g<g_decomp->n_groups;g++) apply_dirichlet(&g_gfields[g],g);
}

static void copy_curr_to_prev_all(void) {
  for(int g=0;g<g_decomp->n_groups;g++)
    memcpy(g_gfields[g].u_prev,g_gfields[g].u_curr,g_gfields[g].plane_bytes);
}

/* ═══════════════════════════════════════════════════════════════
 *  计算核心：内联到 parallel-for 中消除函数调用开销
 * ═══════════════════════════════════════════════════════════════ */

static int should_measure_energy(int step) {
  if(step==0||step==NT-1)return 1;
  if(ENERGY_REPORT_INTERVAL>0&&((step+1)%ENERGY_REPORT_INTERVAL)==0)return 1;
  return 0;
}

/* 前向声明 */
static double compute_energy_block(GroupField *gf, int gid,
    int y_begin, int y_end, int x_begin, int x_end,
    double *l2_out, double *max_out,
    int *max_x, int *max_y, const int *gy_cache);
static void halo_exchange_mpi(GroupField *gf, int gid);

/* 预计算不变量，避免 compute_region 中重复计算 */
static void run_simulation(int nthreads) {
  int gid=0;
  GroupField *gf=&g_gfields[gid];
  const mythread_tile *tile=&g_decomp->group_tiles[gid];
  int ny_int=tile->ny, nx_int=tile->nx, stride=gf->stride;
  double local_energy, global_energy;

  /* 预计算：每个 y 的全局坐标，避免重复调用 global_y_from_local */
  int *gy_cache = (int*)malloc((size_t)(ny_int+2*HALO)*sizeof(int));
  for(int y=0;y<ny_int+2*HALO;y++) gy_cache[y]=global_y_from_local(gid,y);
  /* 预计算：每个 x 的全局坐标 */
  int *gx_cache = (int*)malloc((size_t)(nx_int+2*HALO)*sizeof(int));
  for(int x=0;x<nx_int+2*HALO;x++) gx_cache[x]=global_x_from_local(gid,x);

  /* 域边界标记（compute_interior 判断用） */
  int nd=g_sim.neighbor_down>=0&&mythread_decomp_is_domain_boundary(g_decomp,gid,MYTHREAD_NEIGHBOR_DOWN);
  int nu=g_sim.neighbor_up>=0&&mythread_decomp_is_domain_boundary(g_decomp,gid,MYTHREAD_NEIGHBOR_UP);
  int nl=g_sim.neighbor_left>=0&&mythread_decomp_is_domain_boundary(g_decomp,gid,MYTHREAD_NEIGHBOR_LEFT);
  int nr=g_sim.neighbor_right>=0&&mythread_decomp_is_domain_boundary(g_decomp,gid,MYTHREAD_NEIGHBOR_RIGHT);
  int y_start_interior=nd?HALO+1:HALO;
  int y_end_interior=nu?ny_int:ny_int+HALO;
  int x_start_interior=nl?HALO+1:HALO;
  int x_end_interior=nr?nx_int:nx_int+HALO;
  /* 边界行 */
  int lr=HALO, ur=ny_int, lc=HALO, rc=nx_int;

  /* CFL 常量 */
  double inv_dx2=1.0/(DX*DX), inv_dy2=1.0/(DY*DY), c02_dt2=C0*C0*DT2;

  if(gf_alloc(gf,gid,g_decomp)!=0){
    fprintf(stderr,"[Error] MPI=%d: gf_alloc failed\n",g_sim.mpi_rank);
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  /* OpenMP: 嵌套 {} 块内的变量自动 private。reduction 需要的
     pe 必须在 parallel 外部声明为 shared 才能被 omp for 归约 */
  double pe;
#pragma omp parallel num_threads(nthreads) shared(pe)
  {
    int tid=omp_get_thread_num();
    double *tm=&g_times[TM_SLOT(tid)];  /* cache-line padded */
    double t0, l2, ma; int mx_, my_;

    /* ── Phase 1: 初始化波场 ── */
#pragma omp single
    t0=wall_time();
    {
#pragma omp for schedule(static,16) nowait
      for(int y=HALO;y<HALO+ny_int;y++){
        int gy=gy_cache[y];
        for(int x=HALO;x<HALO+nx_int;x++){
          int gx=gx_cache[x];
          size_t p=GFIDX(gf,y,x);
          double v=(gx!=0&&gx!=NX-1&&gy!=0&&gy!=NY-1)?initial_condition_value(gy,gx):0.0;
          gf->u_curr[p]=v;gf->u_prev[p]=v;gf->u_next[p]=0.0;
        }
      }
    }
#pragma omp single
    tm[TM_COMP]+=wall_time()-t0;

    /* ── 首次 Dirichlet + Halo ── */
#pragma omp single
    {
      apply_dirichlet_all();
      t0=wall_time();
      halo_exchange_mpi(gf,gid);
      tm[TM_HALO]+=wall_time()-t0;
      apply_dirichlet_all();
      copy_curr_to_prev_all();
    }

    /* ── Phase 2: 初始能量 ── */
#pragma omp single
    { t0=wall_time();gf->group_energy=0.0;g_l2_acc[0]=0.0;g_max_acc[0]=0.0;pe=0.0;}
    {
      l2=ma=0.0;mx_=my_=0;
#pragma omp for schedule(static,16) reduction(+:pe)
      for(int y=HALO;y<HALO+ny_int;y++){
        double l2_,ma_;int mx__,my__;
        pe+=compute_energy_block(gf,gid,y,y+1,HALO,HALO+nx_int,&l2_,&ma_,&mx__,&my__,gy_cache);
        l2+=l2_;if(ma_>ma){ma=ma_;mx_=mx__;my_=my__;}
      }
#pragma omp critical
      {gf->group_energy+=pe;g_l2_acc[0]+=l2;
       if(ma>g_max_acc[0]){g_max_acc[0]=ma;g_max_x[0]=mx_;g_max_y[0]=my_;}}
    }
#pragma omp single
    {
      tm[TM_ENERGY]+=wall_time()-t0;
      local_energy=gf->group_energy;
      MPI_Allreduce(&local_energy,&global_energy,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      g_sim.initial_energy=global_energy;
      if(g_sim.mpi_rank==0){
        double l2=sqrt(g_l2_acc[0]*DX*DY);
        printf("[Main] Initial: E=%.6f L2=%.6f max|u|=%.6f@(%d,%d)\n",
               global_energy,l2,g_max_acc[0],g_max_x[0],g_max_y[0]);
      }
    }

    /* ═══════════════════════════════════════════════════════════
     *  时间步循环（优化版：barrier 从每步 ~10 次降到 ~3 次）
     * ═══════════════════════════════════════════════════════════ */
    double prev_time=wall_time();
    for(int step=0;step<NT;step++){
      MPI_Barrier(MPI_COMM_WORLD);
      int need_e=should_measure_energy(step);

      /* (a) Dirichlet + Halo（单线程执行，omp for 隐式 barrier 保证安全） */
#pragma omp single
      {
        apply_dirichlet_all();
        t0=wall_time();halo_exchange_mpi(gf,gid);tm[TM_HALO]+=wall_time()-t0;
        apply_dirichlet_all();
      }
      /* single 末尾的隐式 barrier 确保 halo 就绪 */

      /* (b) Compute interior + boundary 合并为一个 parallel-for */
#pragma omp single
      t0=wall_time();
      {
        /* interior: y ∈ [y_start, y_end), x ∈ [x_start, x_end) */
#pragma omp for schedule(static,16) nowait
        for(int y=y_start_interior;y<y_end_interior;y++){
          int gy=gy_cache[y];
          if(gy==0||gy==NY-1)continue;
          double *u_prev=gf->u_prev, *u_curr=gf->u_curr, *u_next=gf->u_next;
          for(int x=x_start_interior;x<x_end_interior;x++){
            int gx=gx_cache[x];
            if(gx==0||gx==NX-1)continue;
            size_t p=(size_t)y*stride+(size_t)x;
            double u_ij=u_curr[p];
            double d2x=(u_curr[p-1]-2.0*u_ij+u_curr[p+1])*inv_dx2;
            double d2y=(u_curr[p-stride]-2.0*u_ij+u_curr[p+stride])*inv_dy2;
            u_next[p]=2.0*u_ij-u_prev[p]+c02_dt2*(d2x+d2y);
          }
        }

        /* boundary rows (MPI neighbor directions) */
        if(nd){
#pragma omp for schedule(static) nowait
          for(int x=HALO;x<HALO+nx_int;x++){
            int gx=gx_cache[x];if(gx==0||gx==NX-1)continue;
            double *u_prev=gf->u_prev,*u_curr=gf->u_curr,*u_next=gf->u_next;
            size_t p=(size_t)lr*stride+(size_t)x;
            double u_ij=u_curr[p];
            double d2x=(u_curr[p-1]-2.0*u_ij+u_curr[p+1])*inv_dx2;
            double d2y=(u_curr[p-stride]-2.0*u_ij+u_curr[p+stride])*inv_dy2;
            u_next[p]=2.0*u_ij-u_prev[p]+c02_dt2*(d2x+d2y);
          }
        }
        if(nu){
#pragma omp for schedule(static) nowait
          for(int x=HALO;x<HALO+nx_int;x++){
            int gx=gx_cache[x];if(gx==0||gx==NX-1)continue;
            double *u_prev=gf->u_prev,*u_curr=gf->u_curr,*u_next=gf->u_next;
            size_t p=(size_t)ur*stride+(size_t)x;
            double u_ij=u_curr[p];
            double d2x=(u_curr[p-1]-2.0*u_ij+u_curr[p+1])*inv_dx2;
            double d2y=(u_curr[p-stride]-2.0*u_ij+u_curr[p+stride])*inv_dy2;
            u_next[p]=2.0*u_ij-u_prev[p]+c02_dt2*(d2x+d2y);
          }
        }
        if(nl){
#pragma omp for schedule(static) nowait
          for(int y=HALO;y<HALO+ny_int;y++){
            int gy=gy_cache[y];if(gy==0||gy==NY-1)continue;
            double *u_prev=gf->u_prev,*u_curr=gf->u_curr,*u_next=gf->u_next;
            size_t p=(size_t)y*stride+(size_t)lc;
            double u_ij=u_curr[p];
            double d2x=(u_curr[p-1]-2.0*u_ij+u_curr[p+1])*inv_dx2;
            double d2y=(u_curr[p-stride]-2.0*u_ij+u_curr[p+stride])*inv_dy2;
            u_next[p]=2.0*u_ij-u_prev[p]+c02_dt2*(d2x+d2y);
          }
        }
        if(nr){
#pragma omp for schedule(static) nowait
          for(int y=HALO;y<HALO+ny_int;y++){
            int gy=gy_cache[y];if(gy==0||gy==NY-1)continue;
            double *u_prev=gf->u_prev,*u_curr=gf->u_curr,*u_next=gf->u_next;
            size_t p=(size_t)y*stride+(size_t)rc;
            double u_ij=u_curr[p];
            double d2x=(u_curr[p-1]-2.0*u_ij+u_curr[p+1])*inv_dx2;
            double d2y=(u_curr[p-stride]-2.0*u_ij+u_curr[p+stride])*inv_dy2;
            u_next[p]=2.0*u_ij-u_prev[p]+c02_dt2*(d2x+d2y);
          }
        }
      }
#pragma omp single
      tm[TM_COMP]+=wall_time()-t0;

      /* (c) Swap fields */
#pragma omp single
      gf_swap(gf);

      /* (d) Energy (按需) */
      if(need_e){
#pragma omp single
        {
          apply_dirichlet_all();
          t0=wall_time();halo_exchange_mpi(gf,gid);tm[TM_HALO]+=wall_time()-t0;
          apply_dirichlet_all();
          gf->group_energy=0.0;g_l2_acc[0]=0.0;g_max_acc[0]=0.0;pe=0.0;
        }
        {
          l2=ma=0.0;mx_=my_=0;
#pragma omp for schedule(static,16) reduction(+:pe)
          for(int y=HALO;y<HALO+ny_int;y++){
            double l2_,ma_;int mx__,my__;
            pe+=compute_energy_block(gf,gid,y,y+1,HALO,HALO+nx_int,&l2_,&ma_,&mx__,&my__,gy_cache);
            l2+=l2_;if(ma_>ma){ma=ma_;mx_=mx__;my_=my__;}
          }
#pragma omp critical
          {gf->group_energy+=pe;g_l2_acc[0]+=l2;
           if(ma>g_max_acc[0]){g_max_acc[0]=ma;g_max_x[0]=mx_;g_max_y[0]=my_;}}
        }
#pragma omp single
        {
          tm[TM_ENERGY]+=wall_time()-t0;
          local_energy=gf->group_energy;
          MPI_Allreduce(&local_energy,&global_energy,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
          if(g_sim.mpi_rank==0){
            double cur_time=MPI_Wtime();
            printf("[Main] Step %4d/%d, time %.3f,  E=%.6f L2=%.6f max|u|=%.6f@(%d,%d)\n",
                   step+1,NT,cur_time-prev_time,global_energy,
                   sqrt(g_l2_acc[0]*DX*DY),g_max_acc[0],g_max_x[0],g_max_y[0]);
            prev_time=cur_time;
          }
        }
      }
    } /* step loop */
  } /* omp parallel */

  free(gy_cache);free(gx_cache);

  double elapsed=wall_time();
  if(g_sim.mpi_rank==0){
    printf("[Main] Simulation completed in %.3f seconds\n",elapsed);
    printf("[Main] Throughput: %.2f Mpoint-updates/s\n",(double)NX*NY*NT/elapsed/1.0e6);
  }
}

/* ═══════════════════════════════════════════════════════════════
 *  compute_energy_block：带 gy_cache 的版本（避免重复计算 global_y）
 * ═══════════════════════════════════════════════════════════════ */
static double compute_energy_block(GroupField *gf, int gid,
                                    int y_begin, int y_end, int x_begin, int x_end,
                                    double *l2_out, double *max_out,
                                    int *max_x, int *max_y, const int *gy_cache) {
  if(x_begin>=x_end||y_begin>=y_end){
    if(l2_out)*l2_out=0.0;if(max_out)*max_out=0.0;return 0.0;
  }
  const mythread_tile *tile=&g_decomp->group_tiles[gid];
  double ke=0.0,px=0.0,py=0.0,ca=DX*DY,l2=0.0,ma=0.0;
  int mx=0,my=0;
  int xeb=(x_begin==HALO)?(HALO-1):x_begin,xee=x_end;
  if(xee>HALO+tile->nx)xee=HALO+tile->nx;
  for(int y=y_begin;y<y_end;y++){
    int gy=gy_cache[y];
    for(int x=x_begin;x<x_end;x++){
      double u=gf->u_curr[GFIDX(gf,y,x)];l2+=u*u;
      double au=fabs(u);if(au>ma){ma=au;mx=global_x_from_local(gid,x);my=gy;}
    }
    if(gy>0&&gy<NY-1)
      for(int x=x_begin;x<x_end;x++){
        double ut=(gf->u_curr[GFIDX(gf,y,x)]-gf->u_prev[GFIDX(gf,y,x)])/DT;ke+=ut*ut;
      }
    for(int x=xeb;x<xee;x++){
      double dc=(gf->u_curr[GFIDX(gf,y,x+1)]-gf->u_curr[GFIDX(gf,y,x)])/DX;
      double dp=(gf->u_prev[GFIDX(gf,y,x+1)]-gf->u_prev[GFIDX(gf,y,x)])/DX;px+=dc*dp;
    }
    if(gy<NY-1)
      for(int x=x_begin;x<x_end;x++){
        double dc=(gf->u_curr[GFIDX(gf,y+1,x)]-gf->u_curr[GFIDX(gf,y,x)])/DY;
        double dp=(gf->u_prev[GFIDX(gf,y+1,x)]-gf->u_prev[GFIDX(gf,y,x)])/DY;py+=dc*dp;
      }
  }
  if(l2_out)*l2_out=l2;if(max_out)*max_out=ma;if(max_x)*max_x=mx;if(max_y)*max_y=my;
  return 0.5*(ke+C0*C0*(px+py))*ca;
}

/* ── Halo 交换（内联 MPI）── */
static void halo_exchange_mpi(GroupField *gf, int gid) {
  const mythread_tile *tile=&g_decomp->group_tiles[gid];
  int nx=tile->nx,ny=tile->ny,s=gf->stride,tag=g_mpi_ctx.mpi_tag_base;
  MPI_Request reqs[4];int nr=0;

  if(g_sim.neighbor_up>=0&&gf->send_up&&gf->recv_up){
    memcpy(gf->send_up,&gf->u_curr[ny*s+HALO],(size_t)nx*sizeof(double));
    MPI_Isend(gf->send_up,nx,MPI_DOUBLE,g_sim.neighbor_up,tag,MPI_COMM_WORLD,&reqs[nr++]);
  }
  if(g_sim.neighbor_down>=0&&gf->send_down&&gf->recv_down){
    MPI_Irecv(gf->recv_down,nx,MPI_DOUBLE,g_sim.neighbor_down,tag,MPI_COMM_WORLD,&reqs[nr++]);
  }
  if(g_sim.neighbor_down>=0&&gf->send_down&&gf->recv_down){
    memcpy(gf->send_down,&gf->u_curr[HALO*s+HALO],(size_t)nx*sizeof(double));
    MPI_Isend(gf->send_down,nx,MPI_DOUBLE,g_sim.neighbor_down,tag+1,MPI_COMM_WORLD,&reqs[nr++]);
  }
  if(g_sim.neighbor_up>=0&&gf->send_up&&gf->recv_up){
    MPI_Irecv(gf->recv_up,nx,MPI_DOUBLE,g_sim.neighbor_up,tag+1,MPI_COMM_WORLD,&reqs[nr++]);
  }
  if(g_sim.neighbor_left>=0&&gf->send_left&&gf->recv_left){
    for(int row=0;row<ny;row++)gf->send_left[row]=gf->u_curr[(HALO+row)*s+HALO];
    MPI_Isend(gf->send_left,ny,MPI_DOUBLE,g_sim.neighbor_left,tag+2,MPI_COMM_WORLD,&reqs[nr++]);
    MPI_Irecv(gf->recv_left,ny,MPI_DOUBLE,g_sim.neighbor_left,tag+3,MPI_COMM_WORLD,&reqs[nr++]);
  }
  if(g_sim.neighbor_right>=0&&gf->send_right&&gf->recv_right){
    int src_col=HALO+nx-1;
    for(int row=0;row<ny;row++)gf->send_right[row]=gf->u_curr[(HALO+row)*s+src_col];
    MPI_Isend(gf->send_right,ny,MPI_DOUBLE,g_sim.neighbor_right,tag+3,MPI_COMM_WORLD,&reqs[nr++]);
    MPI_Irecv(gf->recv_right,ny,MPI_DOUBLE,g_sim.neighbor_right,tag+2,MPI_COMM_WORLD,&reqs[nr++]);
  }

  MPI_Waitall(nr,reqs,MPI_STATUSES_IGNORE);

  if(g_sim.neighbor_down>=0&&gf->recv_down)
    memcpy(&gf->u_curr[0*s+HALO],gf->recv_down,(size_t)nx*sizeof(double));
  if(g_sim.neighbor_up>=0&&gf->recv_up)
    memcpy(&gf->u_curr[(ny+HALO)*s+HALO],gf->recv_up,(size_t)nx*sizeof(double));
  if(g_sim.neighbor_left>=0&&gf->recv_left)
    for(int row=0;row<ny;row++)gf->u_curr[(HALO+row)*s+0]=gf->recv_left[row];
  if(g_sim.neighbor_right>=0&&gf->recv_right){
    int dst_col=HALO+nx;
    for(int row=0;row<ny;row++)gf->u_curr[(HALO+row)*s+dst_col]=gf->recv_right[row];
  }
}

/* ── 配置加载 ── */
static void cfg_load_from_files(const char *case_path, const char *hw_path) {
  mythread_cfg *hw=mythread_cfg_load(hw_path);
  mythread_cfg *cs=mythread_cfg_load(case_path);
  cfg_NX=cs?mythread_cfg_get_int(cs,"","NX",cfg_NX):cfg_NX;
  cfg_NY=cs?mythread_cfg_get_int(cs,"","NY",cfg_NY):cfg_NY;
  cfg_NT=cs?mythread_cfg_get_int(cs,"","NT",cfg_NT):cfg_NT;
  cfg_A=cs?mythread_cfg_get_double(cs,"","A",cfg_A):cfg_A;
  cfg_DT=cs?mythread_cfg_get_double(cs,"","DT",cfg_DT):cfg_DT;
  cfg_C0=cs?mythread_cfg_get_double(cs,"","C0",cfg_C0):cfg_C0;
  cfg_USE_FIXED_DOMAIN=cs?mythread_cfg_get_int(cs,"","USE_FIXED_DOMAIN",cfg_USE_FIXED_DOMAIN):cfg_USE_FIXED_DOMAIN;
  cfg_HALO=cs?mythread_cfg_get_int(cs,"","HALO",cfg_HALO):cfg_HALO;
  cfg_ENERGY_REPORT_INTERVAL=cs?mythread_cfg_get_int(cs,"","ENERGY_REPORT_INTERVAL",cfg_ENERGY_REPORT_INTERVAL):cfg_ENERGY_REPORT_INTERVAL;
  cfg_N_WORKERS=hw?mythread_cfg_get_int(hw,"","N_WORKERS",cfg_N_WORKERS):cfg_N_WORKERS;
  cfg_GROUP_DECOMP=hw?mythread_cfg_get_int(hw,"","GROUP_DECOMP",cfg_GROUP_DECOMP):cfg_GROUP_DECOMP;
  mythread_cfg_free(cs);mythread_cfg_free(hw);
  cfg_compute_derived();
}

/* ── 计时 ── */
static void print_timing_report(int mpi_rank, int mpi_size, int nthreads) {
  static const char *names[]={"comp","energy","halo"};
  for(int r=0;r<mpi_size;r++){
    MPI_Barrier(MPI_COMM_WORLD);
    if(mpi_rank==r){
      for(int t=0;t<nthreads;t++){
        double *tm=&g_times[TM_SLOT(t)];
        double tot=tm[TM_COMP]+tm[TM_ENERGY]+tm[TM_HALO];
        printf("[Timing] rank %d/%d tid=%d ",mpi_rank,mpi_size,t);
        for(int s=0;s<3;s++)printf("%s=%.3f ",names[s],tm[s]);
        printf("tot=%.3f\n",tot);
      }
      fflush(stdout);
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

/* ═══════════════════════════════════════════════════════════════
 *  main
 * ═══════════════════════════════════════════════════════════════ */
int main(int argc, char **argv) {
  int mpi_rank, mpi_size, node_size;
  int local_ready=1, global_ready=1;

  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);

  MPI_Comm node_comm;
  MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,0,MPI_INFO_NULL,&node_comm);
  MPI_Comm_size(node_comm,&node_size);

  {
    const char *cp=mythread_env_get("WAVE_CASE_CFG","config/case.cfg");
    const char *hp=mythread_env_get("WAVE_HARDWARE_CFG","config/hardware.cfg");
    cfg_load_from_files(cp,hp);
  }
  g_nthreads=cfg_N_WORKERS>0?cfg_N_WORKERS:omp_get_max_threads();
  omp_set_num_threads(g_nthreads);

  if(NX<3||NY<3){if(mpi_rank==0)fprintf(stderr,"[Error] NX/NY >=3\n");MPI_Finalize();return 1;}
  if((long long)mpi_size>(long long)NX*NY){if(mpi_rank==0)fprintf(stderr,"[Error] too many procs\n");MPI_Finalize();return 1;}
  if(CFL_SUM2>1.0){if(mpi_rank==0)fprintf(stderr,"[Error] CFL>1\n");MPI_Finalize();return 1;}

  memset(&g_sim,0,sizeof(g_sim));
  g_sim.mpi_rank=mpi_rank;g_sim.mpi_size=mpi_size;
  setup_process_domain(mpi_rank,mpi_size);
  if(g_sim.local_nx<=0||g_sim.local_ny<=0)local_ready=0;
  MPI_Allreduce(&local_ready,&global_ready,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  if(!global_ready){MPI_Finalize();return 1;}

  g_decomp=mythread_decomp_create(g_sim.local_nx,g_sim.local_ny,HALO,1,g_nthreads,cfg_GROUP_DECOMP);
  if(!g_decomp){fprintf(stderr,"[Error] decomp_create failed\n");MPI_Finalize();return 1;}
  if(!validate_memory(mpi_rank,node_size,g_decomp)){
    mythread_decomp_free((mythread_decomp*)g_decomp);MPI_Finalize();return 1;
  }

  g_gfields=(GroupField*)xcalloc((size_t)g_decomp->n_groups,sizeof(GroupField));
  g_l2_acc =(double*)xcalloc((size_t)g_decomp->n_groups,sizeof(double));
  g_max_acc=(double*)xcalloc((size_t)g_decomp->n_groups,sizeof(double));
  g_max_x  =(int*)xcalloc((size_t)g_decomp->n_groups,sizeof(int));
  g_max_y  =(int*)xcalloc((size_t)g_decomp->n_groups,sizeof(int));
  /* cache-line padded: each thread's timing slots are TM_N_SLOTS_PAD=8 doubles = 64 bytes */
  g_times  =(double*)xcalloc((size_t)g_nthreads*TM_N_SLOTS_PAD,sizeof(double));

  if(mpi_rank==0){
    printf("============================================\n");
    printf("  Wave Equation (MPI + OpenMP, optimized)\n");
    printf("============================================\n");
    printf("Global grid       : %d x %d\n",NX,NY);
    printf("Time steps        : %d\n",NT);
    printf("MPI processes     : %d (%d x %d)\n",mpi_size,g_sim.proc_px,g_sim.proc_py);
    printf("OpenMP threads    : %d\n",g_nthreads);
    printf("DT=%.4f DX=%.4f CFL: x=%.4f y=%.4f sum2=%.6f\n",cfg_DT,DX,CFL_X,CFL_Y,CFL_SUM2);
    printf("============================================\n");
  }
  MPI_Barrier(MPI_COMM_WORLD);
  printf("[Domain] rank %d/%d node_size=%d proc=(%d,%d)/(%d,%d) local=(nx=%d,ny=%d) "
         "neighbors(L=%d R=%d D=%d U=%d)\n",
         mpi_rank,mpi_size,node_size,g_sim.proc_x,g_sim.proc_y,
         g_sim.proc_px,g_sim.proc_py,g_sim.local_nx,g_sim.local_ny,
         g_sim.neighbor_left,g_sim.neighbor_right,
         g_sim.neighbor_down,g_sim.neighbor_up);
  MPI_Barrier(MPI_COMM_WORLD);

  double t_start=MPI_Wtime();
  run_simulation(g_nthreads);
  double t_elapsed=MPI_Wtime()-t_start;

  print_timing_report(mpi_rank,mpi_size,g_nthreads);
  if(mpi_rank==0)printf("[Main] Total wall time: %.3f seconds\n",t_elapsed);

  for(int g=0;g<g_decomp->n_groups;g++)gf_free(&g_gfields[g]);
  free(g_gfields);
  mythread_decomp_free((mythread_decomp*)g_decomp);
  free(g_l2_acc);free(g_max_acc);free(g_max_x);free(g_max_y);free(g_times);
  MPI_Finalize();
  return 0;
}
