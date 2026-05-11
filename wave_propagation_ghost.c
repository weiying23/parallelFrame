#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#if defined(__linux__)
#include <sys/sysinfo.h>
#endif

#include "mythread/mythread.h"

/*
 * 二维波动方程 — NUMA 感知 + 任务池调度
 *
 * - 每个线程组分配独立 GroupField，物理内存绑在组所在 NUMA 节点
 * - 组间/MPI halo 交换通过 mythread_halo 统一处理
 * - Worker 通过任务池动态竞争 CHUNK_ROWS 行粒度的 RowTaskCtx
 */

/* ── 运行时参数（编译默认值 → config 文件 → 环境变量）── */
static int    cfg_NX=14000, cfg_NY=14000, cfg_NT=480;
static double cfg_A=10.0, cfg_DT=0.01, cfg_C0=0.1;
static int    cfg_USE_FIXED_DOMAIN=0;
static double cfg_DX, cfg_DY, cfg_LX, cfg_LY, cfg_DT2, cfg_CFL_X, cfg_CFL_Y, cfg_CFL_SUM2;
static int    cfg_HALO=1, cfg_N_GROUPS=2, cfg_N_WORKERS=3, cfg_THREADS_PER_GROUP=4;
static int    cfg_ENERGY_REPORT_INTERVAL=1200, cfg_GROUP_DECOMP=0;
static int    cfg_NCorePClu=5, cfg_NCluPNode=2, cfg_NCorePGrp=4, cfg_ManageCoreId=4;
static int    cfg_CoreOffset=0, cfg_ClustOffset=0;

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
  cfg_THREADS_PER_GROUP = cfg_N_WORKERS + 1;
}

static void cfg_load_from_files(const char *case_path, const char *hw_path) {
  mythread_cfg *hw = mythread_cfg_load(hw_path);
  mythread_cfg *cs = mythread_cfg_load(case_path);
  cfg_NX  = cs ? mythread_cfg_get_int(cs, "", "NX",  cfg_NX)  : cfg_NX;
  cfg_NY  = cs ? mythread_cfg_get_int(cs, "", "NY",  cfg_NY)  : cfg_NY;
  cfg_NT  = cs ? mythread_cfg_get_int(cs, "", "NT",  cfg_NT)  : cfg_NT;
  cfg_A   = cs ? mythread_cfg_get_double(cs, "", "A",  cfg_A)   : cfg_A;
  cfg_DT  = cs ? mythread_cfg_get_double(cs, "", "DT", cfg_DT)  : cfg_DT;
  cfg_C0  = cs ? mythread_cfg_get_double(cs, "", "C0", cfg_C0)  : cfg_C0;
  cfg_USE_FIXED_DOMAIN = cs ? mythread_cfg_get_int(cs, "", "USE_FIXED_DOMAIN", cfg_USE_FIXED_DOMAIN) : cfg_USE_FIXED_DOMAIN;
  cfg_HALO = cs ? mythread_cfg_get_int(cs, "", "HALO", cfg_HALO) : cfg_HALO;
  cfg_ENERGY_REPORT_INTERVAL = cs ? mythread_cfg_get_int(cs, "", "ENERGY_REPORT_INTERVAL", cfg_ENERGY_REPORT_INTERVAL) : cfg_ENERGY_REPORT_INTERVAL;
  cfg_N_GROUPS  = hw ? mythread_cfg_get_int(hw, "", "N_GROUPS",  cfg_N_GROUPS)  : cfg_N_GROUPS;
  cfg_N_WORKERS = hw ? mythread_cfg_get_int(hw, "", "N_WORKERS", cfg_N_WORKERS) : cfg_N_WORKERS;
  cfg_NCorePClu  = hw ? mythread_cfg_get_int(hw, "", "NCorePClu",  cfg_NCorePClu)  : cfg_NCorePClu;
  cfg_NCluPNode  = hw ? mythread_cfg_get_int(hw, "", "NCluPNode",  cfg_NCluPNode)  : cfg_NCluPNode;
  cfg_NCorePGrp  = hw ? mythread_cfg_get_int(hw, "", "NCorePGrp",  cfg_NCorePGrp)  : cfg_NCorePGrp;
  cfg_ManageCoreId = hw ? mythread_cfg_get_int(hw, "", "ManageCoreId", cfg_ManageCoreId) : cfg_ManageCoreId;
  mythread_cfg_free(cs);
  mythread_cfg_free(hw);
  cfg_compute_derived();
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

#define SLOT_TP 0
#define TP_CAP 512
#define CHUNK_ROWS 16
#define CHUNK_COLS 256

/* ── 数据结构 ── */
typedef struct {
  int local_x_begin, local_x_end, local_nx;
  int local_y_begin, local_y_end, local_ny;
  int mpi_rank, mpi_size;
  int proc_x, proc_y, proc_px, proc_py;
  int neighbor_left, neighbor_right, neighbor_up, neighbor_down;
  double initial_energy;
} SimulationData;

typedef struct {
  int gid, tid, cpu_id, x_begin, x_end, y_begin, y_end;
  GroupField *gf;
  double partial_energy;
  double t_wait_init, t_work_init, t_wait_energy0, t_work_energy0;
  double t_wait_compute, t_work_compute, t_wait_boundary, t_work_boundary;
  double t_wait_energy, t_work_energy, t_comm, t_recv, t_waitr, t_allreduce;
  int energy_steps;
} ThreadTask;

typedef struct {
  int y_begin, y_end, x_begin, x_end;
  int gid;                     /* 所属组 id */
  GroupField *gf;
  double *energy_acc;
  double *l2_acc;
  double *max_amp;
  int    *max_amp_gx, *max_amp_gy;
} RowTaskCtx;

typedef struct {
  RowTaskCtx *tasks;
  size_t n_tasks;
} GroupTaskPlan;

/* ── 全局变量 ── */
static SimulationData g_sim = {0};
static GroupField *g_gfields = NULL;
static const mythread_decomp *g_decomp = NULL;
static mythread_mpi_ctx g_mpi_ctx = {0};
static GroupTaskPlan *g_group_plans = NULL;
static RowTaskCtx *g_flat_tasks = NULL;
static int g_n_flat_tasks = 0;
static double g_energy_acc  = 0.0;
static double g_l2_acc      = 0.0;
static double g_max_amp     = 0.0;
static int    g_max_amp_gx  = 0, g_max_amp_gy = 0;
static double *g_group_l2   = NULL;
static double *g_group_max  = NULL;
static int    *g_group_max_x = NULL, *g_group_max_y = NULL;

int _gettdsize_() { return (int)sizeof(ThreadTask); }
int _getgdsize_() { return 0; }

/* ── 坐标 ── */
static inline int global_x_from_local(int gid, int local_x) {
  return g_sim.local_x_begin + (g_decomp->group_tiles[gid].x_begin + local_x - HALO);
}
static inline int global_y_from_local(int gid, int local_y) {
  return g_sim.local_y_begin + (g_decomp->group_tiles[gid].y_begin + local_y - HALO);
}
static inline double wall_time(void) { return MPI_Wtime(); }

static void reset_task_timers(ThreadTask *task) {
  if (!task) return;
  task->cpu_id = -1; task->partial_energy = 0.0;
  task->t_wait_init = task->t_work_init = 0.0;
  task->t_wait_energy0 = task->t_work_energy0 = 0.0;
  task->t_wait_compute = task->t_work_compute = 0.0;
  task->t_wait_boundary = task->t_work_boundary = 0.0;
  task->t_wait_energy = task->t_work_energy = 0.0;
  task->t_comm = task->t_recv = task->t_waitr = task->t_allreduce = 0.0;
  task->energy_steps = 0;
}

/* ── 同步状态 ── */
static int init_fields_state(void)     { return 1; }
static int initial_energy_state(void)  { return 2; }
static int compute_phase_state(int s)  { return 3 * s + 3; }
static int boundary_phase_state(int s) { return 3 * s + 4; }
static int energy_phase_state(int s)   { return 3 * s + 5; }
static int should_measure_energy_step(int step) {
  if (step == 0 || step == NT - 1) return 1;
  if (ENERGY_REPORT_INTERVAL > 0 && ((step + 1) % ENERGY_REPORT_INTERVAL) == 0) return 1;
  return 0;
}
static int uses_group_threads(void) { return ThreadG != 0; }

/* ── 同步助手 ── */
static void start_phase_from_main(int state) {
  if (uses_group_threads()) mSetGrps(state); else mSetSubs(state);
}
static void wait_phase_from_main(int state) {
  if (uses_group_threads()) mWaitGrps(state); else mWaitSubs(state);
}

/* ── 内存 ── */
static void *xcalloc(size_t count, size_t size) {
  void *p = calloc(count, size);
  if (!p) {
    fprintf(stderr, "[Error] MPI=%d: allocation failed (%zu x %zu)\n",
            g_sim.mpi_rank, count, size);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  return p;
}
static int add_size_checked(size_t a, size_t b, size_t *out) {
  if (a > SIZE_MAX - b) return 0; *out = a + b; return 1;
}
static int multiply_size_checked(size_t a, size_t b, size_t *out) {
  if (a != 0 && b > SIZE_MAX / a) return 0; *out = a * b; return 1;
}
static int validate_per_group_memory(int mpi_rank, int node_size,
                                      const mythread_decomp *dc) {
  size_t total = 0;
  for (int g = 0; g < dc->n_groups; g++) {
    const mythread_tile *t = &dc->group_tiles[g];
    size_t pw, ph, plane, planes, halo_y, halo_x, mpi_bufs, gf_total;
    if (!add_size_checked((size_t)t->nx, 2u*(size_t)dc->halo, &pw)) return 0;
    if (!add_size_checked((size_t)t->ny, 2u*(size_t)dc->halo, &ph)) return 0;
    if (!multiply_size_checked(ph, pw, &plane)) return 0;
    if (!multiply_size_checked(plane, sizeof(double), &planes)) return 0;
    if (!multiply_size_checked(planes, 3u, &planes)) return 0; /* 3 平面 */
    if (!multiply_size_checked((size_t)t->nx, sizeof(double)*2, &halo_y)) return 0; /* send/recv Y */
    if (!multiply_size_checked((size_t)t->ny, sizeof(double)*2, &halo_x)) return 0;
    mpi_bufs = 0;
    if (mythread_decomp_is_domain_boundary(dc, g, MYTHREAD_NEIGHBOR_UP))
      { size_t b; multiply_size_checked((size_t)t->nx, sizeof(double)*2, &b); mpi_bufs += b; }
    if (mythread_decomp_is_domain_boundary(dc, g, MYTHREAD_NEIGHBOR_DOWN))
      { size_t b; multiply_size_checked((size_t)t->nx, sizeof(double)*2, &b); mpi_bufs += b; }
    if (mythread_decomp_is_domain_boundary(dc, g, MYTHREAD_NEIGHBOR_LEFT))
      { size_t b; multiply_size_checked((size_t)t->ny, sizeof(double)*2, &b); mpi_bufs += b; }
    if (mythread_decomp_is_domain_boundary(dc, g, MYTHREAD_NEIGHBOR_RIGHT))
      { size_t b; multiply_size_checked((size_t)t->ny, sizeof(double)*2, &b); mpi_bufs += b; }
    if (!add_size_checked(planes, halo_y, &gf_total)) return 0;
    if (!add_size_checked(gf_total, halo_x, &gf_total)) return 0;
    if (!add_size_checked(gf_total, mpi_bufs, &gf_total)) return 0;
    if (!add_size_checked(total, gf_total, &total)) return 0;
  }
  if (mpi_rank == 0) {
    printf("Estimated rank memory: %.3f GiB (%d groups)\n",
           (double)total/(1024.*1024.*1024.), dc->n_groups);
  }
  if (node_size < 1) node_size = 1;
  size_t nb;
  if (!multiply_size_checked(total, (size_t)node_size, &nb)) return 0;
  if (mpi_rank == 0)
    printf("Estimated node memory: %.3f GiB (%d ranks/node)\n",
           (double)nb/(1024.*1024.*1024.), node_size);
#if defined(__linux__)
  { struct sysinfo info;
    if (sysinfo(&info) == 0) {
      unsigned long long vis = (unsigned long long)info.totalram * (unsigned long long)info.mem_unit;
      if (mpi_rank == 0) printf("Visible node memory : %.3f GiB\n", (double)vis/(1024.*1024.*1024.));
      if ((unsigned long long)nb > vis) {
        if (mpi_rank==0) fprintf(stderr,"[Error] estimated %.3f GiB exceeds visible %.3f GiB\n",
                                 (double)nb/(1024.*1024.*1024.), (double)vis/(1024.*1024.*1024.));
        return 0;
      }
    }
  }
#endif
  return 1;
}

/* ── MPI 拓扑 ── */
static void setup_process_domain(int mpi_rank, int mpi_size) {
  mythread_decomp *tmp = mythread_decomp_create(NX, NY, HALO, mpi_size, 0, MYTHREAD_DECOMP_XY_2D);
  int px = tmp->gx, py = tmp->gy; mythread_decomp_free(tmp);
  g_sim.proc_px = px; g_sim.proc_py = py;
  g_sim.proc_x = mpi_rank % px; g_sim.proc_y = mpi_rank / px;
  { int t_nx=NX/px, r_nx=NX%px, t_ny=NY/py, r_ny=NY%py;
    int ox = g_sim.proc_x*t_nx + (g_sim.proc_x<r_nx?g_sim.proc_x:r_nx);
    int oy = g_sim.proc_y*t_ny + (g_sim.proc_y<r_ny?g_sim.proc_y:r_ny);
    g_sim.local_x_begin=ox; g_sim.local_x_end=ox+t_nx+(g_sim.proc_x<r_nx?1:0);
    g_sim.local_y_begin=oy; g_sim.local_y_end=oy+t_ny+(g_sim.proc_y<r_ny?1:0);
  }
  g_sim.local_nx = g_sim.local_x_end - g_sim.local_x_begin;
  g_sim.local_ny = g_sim.local_y_end - g_sim.local_y_begin;
  g_sim.neighbor_left  = (g_sim.proc_x>0)?(mpi_rank-1):-1;
  g_sim.neighbor_right = (g_sim.proc_x+1<px)?(mpi_rank+1):-1;
  g_sim.neighbor_down  = (g_sim.proc_y>0)?(mpi_rank-px):-1;
  g_sim.neighbor_up    = (g_sim.proc_y+1<py)?(mpi_rank+px):-1;
  g_mpi_ctx.mpirank_up=g_sim.neighbor_up; g_mpi_ctx.mpirank_down=g_sim.neighbor_down;
  g_mpi_ctx.mpirank_left=g_sim.neighbor_left; g_mpi_ctx.mpirank_right=g_sim.neighbor_right;
  g_mpi_ctx.mpi_tag_base=100;
}

static double initial_condition_value(int gy, int gx) {
  double cx=0.5*LX, cy=0.5*LY, sigma=0.06*((LX<LY)?LX:LY);
  return A*exp(-((gx*DX-cx)*(gx*DX-cx)+(gy*DY-cy)*(gy*DY-cy))/(2.*sigma*sigma));
}

/* ── 初始化 ── */
static void init_simulation(int mpi_rank, int mpi_size) {
  memset(&g_sim, 0, sizeof(g_sim));
  g_sim.mpi_rank = mpi_rank; g_sim.mpi_size = mpi_size;
  setup_process_domain(mpi_rank, mpi_size);
}

static void free_simulation(void) {
  if (g_gfields) {
    for (int g=0; g<g_decomp->n_groups; g++) group_free_field(&g_gfields[g]);
    free(g_gfields); g_gfields = NULL;
  }
  mythread_decomp_free((mythread_decomp*)g_decomp); g_decomp = NULL;
  memset(&g_sim, 0, sizeof(g_sim));
}

/* ── Dirichlet 边界 ── */
static void enforce_dirichlet_boundaries(GroupField *gf, int gid) {
  const mythread_tile *tile = &g_decomp->group_tiles[gid];
  int s = gf->stride, h = gf->ny_padded;
  if (g_sim.local_x_begin == 0)
    for (int y=0; y<h; y++) gf->u_curr[GFIDX(gf,y,HALO)] = 0.0;
  if (g_sim.local_x_end == NX) {
    int xr = HALO + tile->nx - 1;
    for (int y=0; y<h; y++) gf->u_curr[GFIDX(gf,y,xr)] = 0.0;
  }
  if (g_sim.neighbor_left < 0 && mythread_decomp_is_domain_boundary(g_decomp,gid,MYTHREAD_NEIGHBOR_LEFT))
    for (int y=0; y<h; y++) gf->u_curr[GFIDX(gf,y,0)] = 0.0;
  if (g_sim.neighbor_right < 0 && mythread_decomp_is_domain_boundary(g_decomp,gid,MYTHREAD_NEIGHBOR_RIGHT)) {
    int xh = HALO + tile->nx;
    for (int y=0; y<h; y++) gf->u_curr[GFIDX(gf,y,xh)] = 0.0;
  }
  if (g_sim.neighbor_down < 0 && mythread_decomp_is_domain_boundary(g_decomp,gid,MYTHREAD_NEIGHBOR_DOWN))
    memset(&gf->u_curr[0], 0, (size_t)s * sizeof(double));
  if (g_sim.neighbor_up < 0 && mythread_decomp_is_domain_boundary(g_decomp,gid,MYTHREAD_NEIGHBOR_UP))
    memset(&gf->u_curr[(tile->ny+HALO)*s], 0, (size_t)s * sizeof(double));
}

static void zero_physical_y_boundaries(GroupField *gf, int gid) {
  if (g_sim.local_y_begin == 0 && mythread_decomp_is_domain_boundary(g_decomp,gid,MYTHREAD_NEIGHBOR_DOWN))
    memset(&gf->u_curr[GFIDX(gf,HALO,0)], 0, (size_t)gf->stride*sizeof(double));
  if (g_sim.local_y_end == NY && mythread_decomp_is_domain_boundary(g_decomp,gid,MYTHREAD_NEIGHBOR_UP))
    memset(&gf->u_curr[GFIDX(gf,g_decomp->group_tiles[gid].ny,0)], 0, (size_t)gf->stride*sizeof(double));
}

static void apply_dirichlet_all_groups(void) {
  for (int g=0; g<g_decomp->n_groups; g++) {
    enforce_dirichlet_boundaries(&g_gfields[g], g);
    zero_physical_y_boundaries(&g_gfields[g], g);
  }
}
static void copy_curr_to_prev_all_groups(void) {
  for (int g=0; g<g_decomp->n_groups; g++)
    memcpy(g_gfields[g].u_prev, g_gfields[g].u_curr, g_gfields[g].plane_bytes);
}

/* ── 任务池回调（GroupField 版本）── */
static void task_init_field(void *ctx) {
  RowTaskCtx *c = (RowTaskCtx*)ctx; GroupField *gf = c->gf; int gid = c->gid;
  for (int y=c->y_begin; y<c->y_end; y++) {
    int gy = global_y_from_local(gid, y);
    for (int x=c->x_begin; x<c->x_end; x++) {
      int gx = global_x_from_local(gid, x);
      size_t p = GFIDX(gf, y, x);
      double v = 0.0;
      if (gx!=0 && gx!=NX-1 && gy!=0 && gy!=NY-1) v = initial_condition_value(gy, gx);
      gf->u_curr[p] = v; gf->u_prev[p] = v; gf->u_next[p] = 0.0;
    }
  }
}

static void task_compute_interior(void *ctx) {
  RowTaskCtx *c = (RowTaskCtx*)ctx; GroupField *gf = c->gf; int gid = c->gid;
  const mythread_tile *tile = &g_decomp->group_tiles[gid];
  int yb=c->y_begin, ye=c->y_end, xb=c->x_begin, xe=c->x_end;

  if (g_sim.neighbor_down>=0 && mythread_decomp_is_domain_boundary(g_decomp,gid,MYTHREAD_NEIGHBOR_DOWN) && yb<HALO+1) yb=HALO+1;
  if (g_sim.neighbor_up>=0 && mythread_decomp_is_domain_boundary(g_decomp,gid,MYTHREAD_NEIGHBOR_UP) && ye>tile->ny) ye=tile->ny;
  if (g_sim.neighbor_left>=0 && mythread_decomp_is_domain_boundary(g_decomp,gid,MYTHREAD_NEIGHBOR_LEFT) && xb<HALO+1) xb=HALO+1;
  if (g_sim.neighbor_right>=0 && mythread_decomp_is_domain_boundary(g_decomp,gid,MYTHREAD_NEIGHBOR_RIGHT) && xe>tile->nx) xe=tile->nx;
  if (yb>=ye || xb>=xe) return;

  for (int y=yb; y<ye; y++) {
    int gy = global_y_from_local(gid, y);
    if (gy==0 || gy==NY-1) continue;
    for (int x=xb; x<xe; x++) {
      int gx = global_x_from_local(gid, x);
      if (gx==0 || gx==NX-1) continue;
      double u_ij = gf->u_curr[GFIDX(gf,y,x)];
      double d2x = (gf->u_curr[GFIDX(gf,y,x-1)] - 2.0*u_ij + gf->u_curr[GFIDX(gf,y,x+1)])/(DX*DX);
      double d2y = (gf->u_curr[GFIDX(gf,y-1,x)] - 2.0*u_ij + gf->u_curr[GFIDX(gf,y+1,x)])/(DY*DY);
      gf->u_next[GFIDX(gf,y,x)] = 2.0*u_ij - gf->u_prev[GFIDX(gf,y,x)] + C0*C0*DT2*(d2x+d2y);
    }
  }
}

static void task_compute_boundary(void *ctx) {
  RowTaskCtx *c = (RowTaskCtx*)ctx; GroupField *gf = c->gf; int gid = c->gid;
  const mythread_tile *tile = &g_decomp->group_tiles[gid];
  int lr=HALO, ur=tile->ny, lc=HALO, rc=tile->nx;
  int xb=c->x_begin, xe=c->x_end, yb=c->y_begin, ye=c->y_end;
  if (xb<HALO) xb=HALO; if (xe>HALO+tile->nx) xe=HALO+tile->nx;
  if (yb<HALO) yb=HALO; if (ye>HALO+tile->ny) ye=HALO+tile->ny;

  int nd = g_sim.neighbor_down >= 0 && mythread_decomp_is_domain_boundary(g_decomp,gid,MYTHREAD_NEIGHBOR_DOWN);
  int nu = g_sim.neighbor_up >= 0 && mythread_decomp_is_domain_boundary(g_decomp,gid,MYTHREAD_NEIGHBOR_UP);
  int nl = g_sim.neighbor_left >= 0 && mythread_decomp_is_domain_boundary(g_decomp,gid,MYTHREAD_NEIGHBOR_LEFT);
  int nr = g_sim.neighbor_right >= 0 && mythread_decomp_is_domain_boundary(g_decomp,gid,MYTHREAD_NEIGHBOR_RIGHT);

  if (nd && yb<=lr && lr<ye) {
    for (int x=xb; x<xe; x++) {
      int gx=global_x_from_local(gid,x); if (gx==0||gx==NX-1) continue;
      double u_ij=gf->u_curr[GFIDX(gf,lr,x)];
      gf->u_next[GFIDX(gf,lr,x)] = 2.0*u_ij - gf->u_prev[GFIDX(gf,lr,x)]
        + C0*C0*DT2*((gf->u_curr[GFIDX(gf,lr,x-1)]-2.0*u_ij+gf->u_curr[GFIDX(gf,lr,x+1)])/(DX*DX)
        +(gf->u_curr[GFIDX(gf,lr-1,x)]-2.0*u_ij+gf->u_curr[GFIDX(gf,lr+1,x)])/(DY*DY));
    }
  }
  if (nu && ur!=lr && yb<=ur && ur<ye) {
    for (int x=xb; x<xe; x++) {
      int gx=global_x_from_local(gid,x); if (gx==0||gx==NX-1) continue;
      double u_ij=gf->u_curr[GFIDX(gf,ur,x)];
      gf->u_next[GFIDX(gf,ur,x)] = 2.0*u_ij - gf->u_prev[GFIDX(gf,ur,x)]
        + C0*C0*DT2*((gf->u_curr[GFIDX(gf,ur,x-1)]-2.0*u_ij+gf->u_curr[GFIDX(gf,ur,x+1)])/(DX*DX)
        +(gf->u_curr[GFIDX(gf,ur-1,x)]-2.0*u_ij+gf->u_curr[GFIDX(gf,ur+1,x)])/(DY*DY));
    }
  }
  if (nl && xb<=lc && lc<xe) {
    for (int y=yb; y<ye; y++) {
      int gy=global_y_from_local(gid,y); if (gy==0||gy==NY-1) continue;
      double u_ij=gf->u_curr[GFIDX(gf,y,lc)];
      gf->u_next[GFIDX(gf,y,lc)] = 2.0*u_ij - gf->u_prev[GFIDX(gf,y,lc)]
        + C0*C0*DT2*((gf->u_curr[GFIDX(gf,y,lc-1)]-2.0*u_ij+gf->u_curr[GFIDX(gf,y,lc+1)])/(DX*DX)
        +(gf->u_curr[GFIDX(gf,y-1,lc)]-2.0*u_ij+gf->u_curr[GFIDX(gf,y+1,lc)])/(DY*DY));
    }
  }
  if (nr && rc!=lc && xb<=rc && rc<xe) {
    for (int y=yb; y<ye; y++) {
      int gy=global_y_from_local(gid,y); if (gy==0||gy==NY-1) continue;
      double u_ij=gf->u_curr[GFIDX(gf,y,rc)];
      gf->u_next[GFIDX(gf,y,rc)] = 2.0*u_ij - gf->u_prev[GFIDX(gf,y,rc)]
        + C0*C0*DT2*((gf->u_curr[GFIDX(gf,y,rc-1)]-2.0*u_ij+gf->u_curr[GFIDX(gf,y,rc+1)])/(DX*DX)
        +(gf->u_curr[GFIDX(gf,y-1,rc)]-2.0*u_ij+gf->u_curr[GFIDX(gf,y+1,rc)])/(DY*DY));
    }
  }
}

static void task_compute_energy(void *ctx) {
  RowTaskCtx *c = (RowTaskCtx*)ctx; GroupField *gf = c->gf; int gid = c->gid;
  const mythread_tile *tile = &g_decomp->group_tiles[gid];
  if (c->x_begin>=c->x_end || c->y_begin>=c->y_end) return;
  double ke=0.0, px=0.0, py=0.0, ca=DX*DY, l2=0.0, ma=0.0;
  struct { int x,y; } maxp = {0,0};
  int xeb=(c->x_begin==HALO)?(HALO-1):c->x_begin, xee=c->x_end;
  if (xee>HALO+tile->nx) xee=HALO+tile->nx;
  for (int y=c->y_begin; y<c->y_end; y++) {
    int gy=global_y_from_local(gid,y);
    for (int x=c->x_begin; x<c->x_end; x++) {
      double u = gf->u_curr[GFIDX(gf,y,x)];
      l2 += u*u;
      double au = fabs(u);
      if (au > ma) { ma = au; maxp.x = global_x_from_local(gid,x); maxp.y = gy; }
    }
    if (gy>0 && gy<NY-1)
      for (int x=c->x_begin; x<c->x_end; x++) {
        double ut=(gf->u_curr[GFIDX(gf,y,x)]-gf->u_prev[GFIDX(gf,y,x)])/DT; ke+=ut*ut;
      }
    for (int x=xeb; x<xee; x++) {
      double dc=(gf->u_curr[GFIDX(gf,y,x+1)]-gf->u_curr[GFIDX(gf,y,x)])/DX;
      double dp=(gf->u_prev[GFIDX(gf,y,x+1)]-gf->u_prev[GFIDX(gf,y,x)])/DX; px+=dc*dp;
    }
    if (gy<NY-1)
      for (int x=c->x_begin; x<c->x_end; x++) {
        double dc=(gf->u_curr[GFIDX(gf,y+1,x)]-gf->u_curr[GFIDX(gf,y,x)])/DY;
        double dp=(gf->u_prev[GFIDX(gf,y+1,x)]-gf->u_prev[GFIDX(gf,y,x)])/DY; py+=dc*dp;
      }
  }
  double partial = 0.5*(ke+C0*C0*(px+py))*ca;
  union { double d; uint64_t i; } old, nw;
  do { old.d=*c->energy_acc; nw.d=old.d+partial; } while (!__sync_bool_compare_and_swap((volatile uint64_t*)c->energy_acc,old.i,nw.i));
  do { old.d=*c->l2_acc; nw.d=old.d+l2; } while (!__sync_bool_compare_and_swap((volatile uint64_t*)c->l2_acc,old.i,nw.i));

  if (c->max_amp) {
    for (;;) {
      old.d = *c->max_amp;
      if (ma <= old.d) break;
      nw.d = ma;
      if (__sync_bool_compare_and_swap((volatile uint64_t*)c->max_amp, old.i, nw.i)) {
        *c->max_amp_gx = maxp.x; *c->max_amp_gy = maxp.y;
        break;
      }
    }
  }
}

/* ── 度量汇总 ── */
static double accumulate_worker_energy(void) {
  double e=0.0;
  if (uses_group_threads()) {
    for (int g=0; g<g_decomp->n_groups; g++) e += g_gfields[g].group_energy;
  } else {
    for (int t=0; t<md.Nthreads-1; t++) e += ((ThreadTask*)md.threads[t].td)->partial_energy;
  }
  return e;
}
static double accumulate_l2(void) {
  double s=0.0;
  if (uses_group_threads()) { for (int g=0; g<g_decomp->n_groups; g++) s+=g_group_l2[g]; }
  else s=g_l2_acc;
  return sqrt(s * DX * DY);
}
static void reduce_max_amp(double *amp, int *gx, int *gy) {
  *amp=0.0; *gx=*gy=0;
  if (uses_group_threads()) {
    for (int g=0; g<g_decomp->n_groups; g++)
      if (g_group_max[g] > *amp) { *amp=g_group_max[g]; *gx=g_group_max_x[g]; *gy=g_group_max_y[g]; }
  } else { *amp=g_max_amp; *gx=g_max_amp_gx; *gy=g_max_amp_gy; }
}
static void reset_metrics(void) {
  g_energy_acc=g_l2_acc=g_max_amp=0.0; g_max_amp_gx=g_max_amp_gy=0;
  if (uses_group_threads()) for (int g=0; g<g_decomp->n_groups; g++)
    g_gfields[g].group_energy=g_group_l2[g]=g_group_max[g]=g_group_max_x[g]=g_group_max_y[g]=0;
}

/* ── 任务分配 ── */
static void setup_group_thread_tasks(void) {
  g_group_plans = (GroupTaskPlan*)calloc((size_t)g_decomp->n_groups, sizeof(*g_group_plans));
  if (!g_group_plans) { fprintf(stderr,"alloc plans failed\n"); MPI_Abort(MPI_COMM_WORLD,1); }

  for (int g=0; g<g_decomp->n_groups; g++) {
    const mythread_tile *gtile = &g_decomp->group_tiles[g];
    GroupField *gf = &g_gfields[g];
    int n_rows = gtile->ny, n_cols = gtile->nx;
    size_t ny_chunks = ((size_t)n_rows+CHUNK_ROWS-1)/CHUNK_ROWS;
    size_t nx_chunks = ((size_t)n_cols+CHUNK_COLS-1)/CHUNK_COLS;
    GroupTaskPlan *plan = &g_group_plans[g];
    plan->n_tasks = ny_chunks * nx_chunks;
    plan->tasks = (RowTaskCtx*)calloc(plan->n_tasks, sizeof(*plan->tasks));
    if (!plan->tasks) { fprintf(stderr,"alloc tasks failed\n"); MPI_Abort(MPI_COMM_WORLD,1); }

    for (size_t iy=0; iy<ny_chunks; iy++) {
      int cy = HALO + (int)iy*CHUNK_ROWS, cye = cy+CHUNK_ROWS;
      if (cye > HALO+n_rows) cye = HALO+n_rows;
      for (size_t ix=0; ix<nx_chunks; ix++) {
        int cx = HALO + (int)ix*CHUNK_COLS, cxe = cx+CHUNK_COLS;
        if (cxe > HALO+n_cols) cxe = HALO+n_cols;
        size_t idx = iy*nx_chunks + ix;
        plan->tasks[idx].y_begin=cy; plan->tasks[idx].y_end=cye;
        plan->tasks[idx].x_begin=cx; plan->tasks[idx].x_end=cxe;
        plan->tasks[idx].gid = g;
        plan->tasks[idx].gf = gf;
        plan->tasks[idx].energy_acc = &gf->group_energy;
        plan->tasks[idx].l2_acc = &g_group_l2[g];
        plan->tasks[idx].max_amp = &g_group_max[g];
        plan->tasks[idx].max_amp_gx = &g_group_max_x[g];
        plan->tasks[idx].max_amp_gy = &g_group_max_y[g];
      }
    }

    threadGroup *pg = md.grps[g];
    for (int t=0; t<pg->Nthreads; t++) {
      ThreadTask *task = (ThreadTask*)pg->threads[t].td;
      task->gid=g; task->tid=t; task->gf=gf;
      task->x_begin=HALO; task->x_end=HALO+gtile->nx;
      task->y_begin=HALO; task->y_end=HALO+gtile->ny;
      reset_task_timers(task);
    }
  }
}

static void setup_single_group_tasks(void) {
  int n_rows = g_decomp->group_tiles[0].ny, n_cols = g_decomp->group_tiles[0].nx;
  int nw = g_decomp->n_workers_per_group;
  size_t ny_chunks = ((size_t)n_rows+CHUNK_ROWS-1)/CHUNK_ROWS;
  size_t nx_chunks = ((size_t)n_cols+CHUNK_COLS-1)/CHUNK_COLS;
  g_n_flat_tasks = (int)(ny_chunks * nx_chunks);
  g_flat_tasks = (RowTaskCtx*)calloc((size_t)g_n_flat_tasks, sizeof(*g_flat_tasks));
  if (!g_flat_tasks) { fprintf(stderr,"alloc flat tasks failed\n"); MPI_Abort(MPI_COMM_WORLD,1); }
  GroupField *gf = &g_gfields[0];
  for (size_t iy=0; iy<ny_chunks; iy++) {
    int cy=HALO+(int)iy*CHUNK_ROWS, cye=cy+CHUNK_ROWS;
    if (cye>HALO+n_rows) cye=HALO+n_rows;
    for (size_t ix=0; ix<nx_chunks; ix++) {
      int cx=HALO+(int)ix*CHUNK_COLS, cxe=cx+CHUNK_COLS;
      if (cxe>HALO+n_cols) cxe=HALO+n_cols;
      size_t idx = iy*nx_chunks + ix;
      g_flat_tasks[idx].y_begin=cy; g_flat_tasks[idx].y_end=cye;
      g_flat_tasks[idx].x_begin=cx; g_flat_tasks[idx].x_end=cxe;
      g_flat_tasks[idx].gid=0;
      g_flat_tasks[idx].gf=gf; g_flat_tasks[idx].energy_acc=&g_energy_acc;
      g_flat_tasks[idx].l2_acc=&g_l2_acc;
      g_flat_tasks[idx].max_amp=&g_max_amp;
      g_flat_tasks[idx].max_amp_gx=&g_max_amp_gx;
      g_flat_tasks[idx].max_amp_gy=&g_max_amp_gy;
    }
  }
  for (int t=0; t<nw; t++) {
    ThreadTask *task=(ThreadTask*)md.threads[t].td; task->gid=0; task->tid=t; task->gf=gf;
    reset_task_timers(task);
  }
}

static void setup_thread_tasks(void) {
  if (uses_group_threads()) setup_group_thread_tasks();
  else setup_single_group_tasks();
}

static void free_task_plans(void) {
  if (g_group_plans) {
    for (int g=0; g<g_decomp->n_groups; g++) free(g_group_plans[g].tasks);
    free(g_group_plans); g_group_plans=NULL;
  }
  free(g_flat_tasks); g_flat_tasks=NULL; g_n_flat_tasks=0;
}

/* ═══════════════════════════════════════════════════════════════
 *  线程入口
 * ═══════════════════════════════════════════════════════════════ */

static void submit_energy_tasks(int slot, GroupTaskPlan *plan) {
  for (size_t i=0; i<plan->n_tasks; i++)
    mt_taskpool_submit(slot, task_compute_energy, &plan->tasks[i]);
}

static void pool_dispatch_flat(int slot, mt_task_fn fn) {
  mt_taskpool_begin(slot);
  for (int i=0; i<g_n_flat_tasks; i++) mt_taskpool_submit(slot, fn, &g_flat_tasks[i]);
  mt_taskpool_close(slot); mt_taskpool_wait(slot);
}

static void worker_thread(void) {
  int typ = uses_group_threads() ? 1 : 2;
  while (!GetLocV(typ, SLOT_TP, NULL)) ntdelay(1);
  mt_taskpool_worker_loop(SLOT_TP);
}

static void group_main_thread(void) {
  ThreadTask *task = (ThreadTask*)ti->td; int gid = ti->igrp, slot=SLOT_TP;
  GroupTaskPlan *plan = &g_group_plans[gid]; double t0;

  if (group_alloc_field(&g_gfields[gid], gid, g_decomp) != 0) {
    fprintf(stderr,"[Error] MPI=%d GMT gid=%d: alloc failed\n",mpi_id,gid);
    MPI_Abort(MPI_COMM_WORLD,1);
  }
  task->gf = &g_gfields[gid];
  mt_taskpool_attach(slot, TP_CAP, MT_TASKPOOL_BLOCK);

  t0=wall_time(); gWaitMain(init_fields_state()); task->t_wait_init+=wall_time()-t0;
  t0=wall_time();
  mt_taskpool_begin(slot);
  for (size_t i=0; i<plan->n_tasks; i++) mt_taskpool_submit(slot, task_init_field, &plan->tasks[i]);
  mt_taskpool_close(slot); mt_taskpool_wait(slot);
  task->t_work_init+=wall_time()-t0; gSetMain(init_fields_state());

  t0=wall_time(); gWaitMain(initial_energy_state()); task->t_wait_energy0+=wall_time()-t0;
  t0=wall_time();
  g_gfields[gid].group_energy = 0.0;
  mt_taskpool_begin(slot); submit_energy_tasks(slot, plan);
  mt_taskpool_close(slot); mt_taskpool_wait(slot);
  task->t_work_energy0+=wall_time()-t0; gSetMain(initial_energy_state());

  for (int step=0; step<NT; step++) {
    int cs=compute_phase_state(step), bs=boundary_phase_state(step), es=energy_phase_state(step);
    t0=wall_time(); gWaitMain(cs); task->t_wait_compute+=wall_time()-t0;
    t0=wall_time();
    mt_taskpool_begin(slot);
    for (size_t i=0; i<plan->n_tasks; i++) mt_taskpool_submit(slot, task_compute_interior, &plan->tasks[i]);
    mt_taskpool_close(slot); mt_taskpool_wait(slot);
    task->t_work_compute+=wall_time()-t0; gSetMain(cs);

    t0=wall_time(); gWaitMain(bs); task->t_wait_boundary+=wall_time()-t0;
    t0=wall_time();
    mt_taskpool_begin(slot);
    for (size_t i=0; i<plan->n_tasks; i++) mt_taskpool_submit(slot, task_compute_boundary, &plan->tasks[i]);
    mt_taskpool_close(slot); mt_taskpool_wait(slot);
    task->t_work_boundary+=wall_time()-t0; gSetMain(bs);

    if (should_measure_energy_step(step)) {
      t0=wall_time(); gWaitMain(es); task->t_wait_energy+=wall_time()-t0;
      t0=wall_time();
      g_gfields[gid].group_energy = 0.0;
      mt_taskpool_begin(slot); submit_energy_tasks(slot, plan);
      mt_taskpool_close(slot); mt_taskpool_wait(slot);
      task->t_work_energy+=wall_time()-t0; task->energy_steps+=1; gSetMain(es);
    }
  }
  mt_taskpool_shutdown(slot); mt_taskpool_detach(slot);
  printf("Group-Summary%d-%d: %.3f %.3f\n",mpi_id,gid,task->t_wait_compute,task->t_work_compute);
}

static void ungrouped_main(void) {
  ThreadTask *task = (ThreadTask*)ti->td;
  if (group_alloc_field(&g_gfields[0],0,g_decomp)!=0) {
    fprintf(stderr,"[Error] MPI=%d: alloc failed\n",mpi_id); MPI_Abort(MPI_COMM_WORLD,1);
  }
  for (int t=0; t<g_decomp->n_workers_per_group; t++) ((ThreadTask*)md.threads[t].td)->gf=&g_gfields[0];
  reset_task_timers(task); int slot=SLOT_TP; double t0;
  mt_taskpool_attach(slot, TP_CAP, MT_TASKPOOL_BLOCK);

  t0=wall_time(); pool_dispatch_flat(slot, task_init_field); task->t_wait_init+=wall_time()-t0;

  apply_dirichlet_all_groups();
  t0=wall_time(); mythread_halo_exchange_intra(g_gfields,g_decomp);
  mythread_halo_exchange_mpi(g_gfields,g_decomp,&g_mpi_ctx,(uintptr_t)MPI_COMM_WORLD);
  apply_dirichlet_all_groups(); task->t_comm+=wall_time()-t0;
  copy_curr_to_prev_all_groups();

  t0=wall_time(); g_energy_acc=g_l2_acc=g_max_amp=0.0; pool_dispatch_flat(slot, task_compute_energy);
  double local_e=g_energy_acc, global_e; task->t_wait_energy0+=wall_time()-t0;
  double l2=accumulate_l2(), ma; int max,may; reduce_max_amp(&ma,&max,&may);
  t0=wall_time(); MPI_Allreduce(&local_e,&global_e,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  task->t_allreduce+=wall_time()-t0; g_sim.initial_energy=global_e;
  if (mpi_id==0) printf("[Main] Initial: E=%.6f L2=%.6f max|u|=%.6f@(%d,%d)\n",
                         g_sim.initial_energy,l2,ma,max,may);

  double start_time=MPI_Wtime(), prev_time=start_time;
  for (int step=0; step<NT; step++) {
    MPI_Barrier(MPI_COMM_WORLD);
    int need_e=should_measure_energy_step(step);

    apply_dirichlet_all_groups();
    t0=wall_time(); mythread_halo_exchange_intra(g_gfields,g_decomp);
    mythread_halo_exchange_mpi(g_gfields,g_decomp,&g_mpi_ctx,(uintptr_t)MPI_COMM_WORLD);
    apply_dirichlet_all_groups(); task->t_comm+=wall_time()-t0;

    t0=wall_time(); pool_dispatch_flat(slot, task_compute_interior); task->t_wait_compute+=wall_time()-t0;
    t0=wall_time(); pool_dispatch_flat(slot, task_compute_boundary); task->t_wait_boundary+=wall_time()-t0;

    for (int g=0; g<g_decomp->n_groups; g++) group_field_swap(&g_gfields[g]);

    if (need_e) {
      apply_dirichlet_all_groups();
      t0=wall_time(); mythread_halo_exchange_intra(g_gfields,g_decomp);
      mythread_halo_exchange_mpi(g_gfields,g_decomp,&g_mpi_ctx,(uintptr_t)MPI_COMM_WORLD);
      apply_dirichlet_all_groups(); task->t_comm+=wall_time()-t0;

      t0=wall_time(); g_energy_acc=g_l2_acc=g_max_amp=0.0;
      pool_dispatch_flat(slot, task_compute_energy);
      local_e=g_energy_acc; task->t_wait_energy+=wall_time()-t0;
      t0=wall_time(); MPI_Allreduce(&local_e,&global_e,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      task->t_allreduce+=wall_time()-t0; task->energy_steps+=1;
      if (mpi_id==0) {
        double ct=MPI_Wtime()-prev_time; prev_time=MPI_Wtime();
        printf("[Main] Step %4d/%d, time %.3f,  E=%.6f L2=%.6f max|u|=%.6f@(%d,%d)\\n",step+1,NT,ct,global_e,accumulate_l2(),g_max_amp,g_max_amp_gx,g_max_amp_gy);
      }
    }
  }
  mt_taskpool_shutdown(slot); mt_taskpool_detach(slot);
  double elapsed=MPI_Wtime()-start_time;
  if (mpi_id==0) {
    printf("[Main] Simulation completed in %.3f seconds\n",elapsed);
    printf("[Main] comm:%.3f compute:%.3f boundary:%.3f\n",task->t_comm,task->t_wait_compute,task->t_wait_boundary);
    printf("[Main] Throughput: %.2f Mpoint-updates/s\n",(double)NX*NY*NT/elapsed/1.0e6);
  }
}

static void main_thread(void) {
  ThreadTask *task = (ThreadTask*)ti->td;
  if (!uses_group_threads()) { ungrouped_main(); return; }
  reset_task_timers(task);
  double t0, t0_val, local_e, global_e;

  t0=wall_time(); start_phase_from_main(init_fields_state()); wait_phase_from_main(init_fields_state());
  task->t_wait_init+=wall_time()-t0; group_field_link_buffers(g_gfields,g_decomp);

  apply_dirichlet_all_groups();
  t0=wall_time(); mythread_halo_exchange_intra(g_gfields,g_decomp);
  mythread_halo_exchange_mpi(g_gfields,g_decomp,&g_mpi_ctx,(uintptr_t)MPI_COMM_WORLD);
  apply_dirichlet_all_groups(); task->t_comm+=wall_time()-t0; copy_curr_to_prev_all_groups();

  t0=wall_time(); start_phase_from_main(initial_energy_state()); wait_phase_from_main(initial_energy_state());
  task->t_wait_energy0+=wall_time()-t0;
  t0_val=wall_time(); local_e=accumulate_worker_energy(); task->t_work_energy0+=wall_time()-t0_val;
  double l2_=accumulate_l2(), ma_; int max_, may_; reduce_max_amp(&ma_,&max_,&may_);

  t0_val=wall_time(); MPI_Allreduce(&local_e,&global_e,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  task->t_allreduce+=wall_time()-t0_val; g_sim.initial_energy=global_e;
  if (mpi_id==0) printf("[Main] Initial: E=%.6f L2=%.6f max|u|=%.6f@(%d,%d)\n",
                         g_sim.initial_energy,l2_,ma_,max_,may_);

  double start_time=MPI_Wtime(), prev_time=start_time;
  for (int step=0; step<NT; step++) {
    MPI_Barrier(MPI_COMM_WORLD);
    int cs=compute_phase_state(step), bs=boundary_phase_state(step), es=energy_phase_state(step);
    int need_e=should_measure_energy_step(step);

    apply_dirichlet_all_groups();
    t0=wall_time(); mythread_halo_exchange_intra(g_gfields,g_decomp);
    mythread_halo_exchange_mpi(g_gfields,g_decomp,&g_mpi_ctx,(uintptr_t)MPI_COMM_WORLD);
    apply_dirichlet_all_groups(); task->t_comm+=wall_time()-t0;

    t0=wall_time(); start_phase_from_main(cs); wait_phase_from_main(cs); task->t_wait_compute+=wall_time()-t0;
    t0=wall_time(); start_phase_from_main(bs); wait_phase_from_main(bs); task->t_wait_boundary+=wall_time()-t0;

    for (int g=0; g<g_decomp->n_groups; g++) group_field_swap(&g_gfields[g]);

    if (need_e) {
      apply_dirichlet_all_groups();
      t0=wall_time(); mythread_halo_exchange_intra(g_gfields,g_decomp);
      mythread_halo_exchange_mpi(g_gfields,g_decomp,&g_mpi_ctx,(uintptr_t)MPI_COMM_WORLD);
      apply_dirichlet_all_groups(); task->t_comm+=wall_time()-t0;

      t0=wall_time(); start_phase_from_main(es); wait_phase_from_main(es); task->t_wait_energy+=wall_time()-t0;
      t0_val=wall_time(); local_e=accumulate_worker_energy(); task->t_work_energy+=wall_time()-t0_val;
      double l2=accumulate_l2(), ma; int max,may; reduce_max_amp(&ma,&max,&may);
      t0_val=wall_time(); MPI_Allreduce(&local_e,&global_e,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      task->t_allreduce+=wall_time()-t0_val; task->energy_steps+=1;
      if (mpi_id==0) {
        double ct=MPI_Wtime()-prev_time; prev_time=MPI_Wtime();
        printf("[Main] Step %4d/%d, time %.3f,  E=%.6f L2=%.6f max|u|=%.6f@(%d,%d)\\n",step+1,NT,ct,global_e,l2,ma,max,may);
      }
    }
  }
  double elapsed=MPI_Wtime()-start_time;
  if (mpi_id==0) {
    printf("[Main] Simulation completed in %.3f seconds\n",elapsed);
    printf("[Main] comm:%.3f compute:%.3f boundary:%.3f\n",task->t_comm,task->t_wait_compute,task->t_wait_boundary);
    printf("[Main] Throughput: %.2f Mpoint-updates/s\n",(double)NX*NY*NT/elapsed/1.0e6);
  }
}

void thread_run(void) {
  ThreadTask *task = (ThreadTask*)ti->td;
  if (task) { task->gid=ti->igrp; task->tid=ti->ind; task->cpu_id=getcpuid(); }
  if (ti->igrp < 0) { main_thread(); }
  else if (!uses_group_threads()) { worker_thread(); }
  else if (ti->ind == 0) { group_main_thread(); }
  else { worker_thread(); }
}

/* ── 计时报告 ── */
static const char *task_role_name(const ThreadTask *task) {
  if (!task) return "null"; if (task->gid<0) return "main";
  if (uses_group_threads() && task->tid==0) return "group-main"; return "worker";
}
static double task_total_time(const ThreadTask *task) {
  if (!task) return 0.0;
  return task->t_wait_init+task->t_work_init+task->t_wait_energy0+task->t_work_energy0
    +task->t_wait_compute+task->t_work_compute+task->t_wait_boundary+task->t_work_boundary
    +task->t_wait_energy+task->t_work_energy+task->t_comm+task->t_allreduce;
}
static void print_one_task_timing(int mr,int ms,int ns,const ThreadTask *task) {
  (void)mr;(void)ms;(void)ns;(void)task;
#ifdef DEBUG
  int gyb=-1,gye=-1; double tot=task_total_time(task);
  if (task && task->y_end>task->y_begin && task->y_begin>=HALO) {
    gyb = global_y_from_local(task->gid,task->y_begin); gye = global_y_from_local(task->gid,task->y_end);
  }
  printf("[Timing] rank %d/%d gid=%d tid=%d role=%s cpu=%d y=[%d,%d) gy=[%d,%d) "
    "w(i=%.4f e0=%.4f c=%.4f b=%.4f e=%.4f) wk(i=%.4f e0=%.4f c=%.4f b=%.4f e=%.4f) "
    "comm=%.4f allreduce=%.4f tot=%.4f\n", mr,ms, task?task->gid:-1,task?task->tid:-1,
    task_role_name(task),task?task->cpu_id:-1, task?task->y_begin:0,task?task->y_end:0,gyb,gye,
    task?task->t_wait_init:0.,task?task->t_wait_energy0:0.,task?task->t_wait_compute:0.,
    task?task->t_wait_boundary:0.,task?task->t_wait_energy:0.,task?task->t_work_init:0.,
    task?task->t_work_energy0:0.,task?task->t_work_compute:0.,task?task->t_work_boundary:0.,
    task?task->t_work_energy:0.,task?task->t_comm:0.,task?task->t_allreduce:0.,tot);
#endif
}
static void print_timing_report(int mr,int ms,int ns) {
  for (int r=0; r<ms; r++) { MPI_Barrier(MPI_COMM_WORLD);
    if (mr==r) {
      print_one_task_timing(mr,ms,ns,(ThreadTask*)md.tm.td);
      if (uses_group_threads()) {
        for (int g=0; g<g_decomp->n_groups; g++)
          for (int t=0; t<md.grps[g]->Nthreads; t++)
            print_one_task_timing(mr,ms,ns,(ThreadTask*)md.grps[g]->threads[t].td);
      } else for (int t=0; t<md.Nthreads-1; t++) print_one_task_timing(mr,ms,ns,(ThreadTask*)md.threads[t].td);
      fflush(stdout);
    }
  } MPI_Barrier(MPI_COMM_WORLD);
}

/* ═══════════════════════════════════════════════════════════════
 *  main
 * ═══════════════════════════════════════════════════════════════ */
int main(int argc,char **argv) {
  int mpi_rank,mpi_size,NCorePClu,NCluPNode,NCorePGrp,NThPGrp,NGrpPProc,NProcPNode=1,ManageCoreId;
  int local_ready=1,global_ready=1,err;
  MPI_Init(&argc,&argv);

  { const char *cp=mythread_env_get("WAVE_CASE_CFG","config/case.cfg");
    const char *hp=mythread_env_get("WAVE_HARDWARE_CFG","config/hardware.cfg");
    cfg_load_from_files(cp,hp); }
  NCorePClu=cfg_NCorePClu; NCluPNode=cfg_NCluPNode; NCorePGrp=cfg_NCorePGrp;
  NThPGrp=cfg_THREADS_PER_GROUP; NGrpPProc=cfg_N_GROUPS; ManageCoreId=cfg_ManageCoreId;

  MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank); MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
  MPI_Comm nc; MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,0,MPI_INFO_NULL,&nc);
  int ns; MPI_Comm_size(nc,&ns); NProcPNode=ns;

  if (NX<3||NY<3) { if(mpi_rank==0) fprintf(stderr,"[Error] NX/NY >=3\n"); MPI_Finalize(); return 1; }
  if ((long long)mpi_size>(long long)NX*NY) { if(mpi_rank==0) fprintf(stderr,"[Error] too many procs\n"); MPI_Finalize(); return 1; }
  if (CFL_SUM2>1.0) { if(mpi_rank==0) fprintf(stderr,"[Error] CFL>1\n"); MPI_Finalize(); return 1; }

  init_simulation(mpi_rank,mpi_size);
  if (g_sim.local_nx<=0||g_sim.local_ny<=0) { fprintf(stderr,"[Error] invalid decomp\n"); local_ready=0; }
  MPI_Allreduce(&local_ready,&global_ready,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  if (!global_ready) { free_simulation(); MPI_Finalize(); return 1; }

  g_decomp = mythread_decomp_create(g_sim.local_nx,g_sim.local_ny,HALO,NGrpPProc,
    uses_group_threads()?(NThPGrp-1):(NThPGrp-1),cfg_GROUP_DECOMP);
  if (!g_decomp) { fprintf(stderr,"[Error] decomp_create failed\n"); MPI_Finalize(); return 1; }
  if (!validate_per_group_memory(mpi_rank, ns, g_decomp)) {
    mythread_decomp_free((mythread_decomp*)g_decomp); g_decomp=NULL; MPI_Finalize(); return 1;
  }
  g_gfields = (GroupField*)xcalloc((size_t)g_decomp->n_groups,sizeof(GroupField));
  g_group_l2   = (double*)xcalloc((size_t)g_decomp->n_groups,sizeof(double));
  g_group_max  = (double*)xcalloc((size_t)g_decomp->n_groups,sizeof(double));
  g_group_max_x = (int*)xcalloc((size_t)g_decomp->n_groups,sizeof(int));
  g_group_max_y = (int*)xcalloc((size_t)g_decomp->n_groups,sizeof(int));

  if (mpi_rank==0) {
    printf("============================================\n");
    printf("  Wave Equation (NUMA-aware, taskpool)\n");
    printf("============================================\n");
    printf("Global grid       : %d x %d\n",NX,NY); printf("Time steps        : %d\n",NT);
    printf("MPI processes     : %d (%d x %d)\n",mpi_size,g_sim.proc_px,g_sim.proc_py);
    printf("Thread groups     : %d\n",NGrpPProc); printf("Threads/group     : %d\n",cfg_THREADS_PER_GROUP);
    printf("Group decomp      : %s\n",g_decomp->policy==MYTHREAD_DECOMP_Y_ONLY?"Y_ONLY":"XY_2D");
    printf("Group grid        : %d x %d\n",g_decomp->gx,g_decomp->gy);
    printf("CFL               : x=%.4f y=%.4f\n",CFL_X,CFL_Y);
    printf("============================================\n");
  }
  MPI_Barrier(MPI_COMM_WORLD);
  printf("[Domain] rank %d/%d ns=%d proc=(%d,%d)/(%d,%d) local=(nx=%d,ny=%d) nbr(L=%d R=%d D=%d U=%d)\n",
    mpi_rank,mpi_size,ns,g_sim.proc_x,g_sim.proc_y,g_sim.proc_px,g_sim.proc_py,
    g_sim.local_nx,g_sim.local_ny,g_sim.neighbor_left,g_sim.neighbor_right,
    g_sim.neighbor_down,g_sim.neighbor_up);
  MPI_Barrier(MPI_COMM_WORLD);

  err=InitThreads(mpi_rank,NCorePClu,NCluPNode,NCorePGrp,NThPGrp,NGrpPProc,NProcPNode,&ManageCoreId);
  if (err!=0) { fprintf(stderr,"[Error] InitThreads failed: %d\n",err); free_simulation(); MPI_Finalize(); return 1; }

  setup_thread_tasks();
  StartThreads(thread_run); thread_run(); EndThreads();
  print_timing_report(mpi_rank,mpi_size,ns);
  free_task_plans(); free_simulation(); MPI_Finalize(); return 0;
}
