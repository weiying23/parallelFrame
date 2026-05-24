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
 * NUMA 感知版二维波动方程。
 *
 * 与旧版的关键区别:
 * - 每个线程组分配独立 GroupField，物理内存绑在组所在 NUMA 节点
 * - 组间 halo 通过 mythread_halo_exchange_intra 在 MMT 中集中完成
 * - MPI halo 通过 mythread_halo_exchange_mpi 自动处理边界组
 * - 域分解使用 mythread_decomp (Y_ONLY / XY_2D 可配置)
 * - Worker 计算 tile 由 mythread_decomp 自动分配
 */

/* ── 运行时参数（编译默认值 → config 文件 → 环境变量）── */
static int    cfg_NX=14000, cfg_NY=14000, cfg_NT=480;
static double cfg_A=10.0, cfg_DT=0.001, cfg_C0=0.1;
static int    cfg_USE_FIXED_DOMAIN=0;
static double cfg_DX, cfg_DY, cfg_LX, cfg_LY, cfg_DT2, cfg_CFL_X, cfg_CFL_Y, cfg_CFL_SUM2;
static int    cfg_HALO=1, cfg_N_GROUPS=2, cfg_N_WORKERS=3, cfg_THREADS_PER_GROUP=4;
static int    cfg_ENERGY_REPORT_INTERVAL=60, cfg_GROUP_DECOMP=0;
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

/* 便捷宏（引用运行时变量） */
#define NX    cfg_NX
#define NY    cfg_NY
#define NT    cfg_NT
#define A     cfg_A
#define DT    cfg_DT
#define C0    cfg_C0
#define USE_FIXED_DOMAIN cfg_USE_FIXED_DOMAIN
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


/* ── 应用数据结构 ── */

typedef struct {
  int local_x_begin, local_x_end, local_nx;
  int local_y_begin, local_y_end, local_ny;
  int mpi_rank, mpi_size;
  int proc_x, proc_y, proc_px, proc_py;
  int neighbor_left, neighbor_right;
  int neighbor_up, neighbor_down;
  double initial_energy;
} SimulationData;

typedef struct {
  int gid, tid;
  int cpu_id;
  int x_begin, x_end;          /* worker tile 在 GroupField 内的局部坐标 */
  int y_begin, y_end;
  GroupField *gf;               /* 本组 GroupField 指针 */
  double partial_energy;
  double partial_l2;
  double partial_max; int partial_max_gx, partial_max_gy;
  double t_wait_init,   t_work_init;
  double t_wait_energy0, t_work_energy0;
  double t_wait_compute, t_work_compute;
  double t_wait_boundary,t_work_boundary;
  double t_wait_energy,  t_work_energy;
  double t_comm, t_recv, t_waitr, t_allreduce;
  int energy_steps;
} ThreadTask;

/* ── 全局变量 ── */
static SimulationData g_sim = {0};
static GroupField *g_gfields = NULL;
static const mythread_decomp *g_decomp = NULL;
static mythread_mpi_ctx g_mpi_ctx = {0};
static double *g_group_l2   = NULL;
static double *g_group_max  = NULL;
static int    *g_group_max_x = NULL, *g_group_max_y = NULL;

int _gettdsize_() { return (int)sizeof(ThreadTask); }
int _getgdsize_() { return 0; }

/* ── 坐标工具 ── */
static inline int global_x_from_local(int gid, int local_x) {
  return g_sim.local_x_begin +
         (g_decomp->group_tiles[gid].x_begin + local_x - HALO);
}

static inline int global_y_from_local(int gid, int local_y) {
  return g_sim.local_y_begin +
         (g_decomp->group_tiles[gid].y_begin + local_y - HALO);
}

static inline double wall_time(void) { return MPI_Wtime(); }

static void reset_task_timers(ThreadTask *task) {
  if (!task) return;
  task->cpu_id = -1;
  task->partial_energy = 0.0;
  task->partial_l2 = 0.0;
  task->partial_max = 0.0; task->partial_max_gx = 0; task->partial_max_gy = 0;
  task->t_wait_init = task->t_work_init = 0.0;
  task->t_wait_energy0 = task->t_work_energy0 = 0.0;
  task->t_wait_compute = task->t_work_compute = 0.0;
  task->t_wait_boundary = task->t_work_boundary = 0.0;
  task->t_wait_energy = task->t_work_energy = 0.0;
  task->t_comm = task->t_recv = task->t_waitr = task->t_allreduce = 0.0;
  task->energy_steps = 0;
}

/* ── 同步状态编码 ── */
static int init_fields_state(void)     { return 1; }
static int initial_energy_state(void)  { return 2; }
static int compute_phase_state(int s)  { return 3 * s + 3; }
static int boundary_phase_state(int s) { return 3 * s + 4; }
static int energy_phase_state(int s)   { return 3 * s + 5; }

static int should_measure_energy_step(int step) {
  if (step == 0 || step == NT - 1) return 1;
  if (ENERGY_REPORT_INTERVAL > 0 && ((step + 1) % ENERGY_REPORT_INTERVAL) == 0)
    return 1;
  return 0;
}

static int uses_group_threads(void) { return ThreadG != 0; }

/* ── 同步助手 ── */
static void start_phase_from_main(int state) {
  if (uses_group_threads()) mSetGrps(state);
  else mSetSubs(state);
}
static void wait_phase_from_main(int state) {
  if (uses_group_threads()) mWaitGrps(state);
  else mWaitSubs(state);
}
static void wait_phase_from_worker(int state) {
  if (uses_group_threads()) sWaitGrp(state);
  else sWaitState(state);
}
static void finish_phase_from_worker(int state) {
  if (uses_group_threads()) sSetGrp(state);
  else sSetState(state);
}

/* ── 内存分配 ── */
static void *xcalloc(size_t count, size_t size) {
  void *p = calloc(count, size);
  if (!p) {
    double gib = ((double)count * (double)size) / (1024.0 * 1024.0 * 1024.0);
    fprintf(stderr, "[Error] MPI=%d: allocation failed (%zu x %zu, %.3f GiB)\n",
            g_sim.mpi_rank, count, size, gib);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  return p;
}

static int add_size_checked(size_t a, size_t b, size_t *out) {
  if (a > SIZE_MAX - b) return 0;
  *out = a + b; return 1;
}
static int multiply_size_checked(size_t a, size_t b, size_t *out) {
  if (a != 0 && b > SIZE_MAX / a) return 0;
  *out = a * b; return 1;
}

static int validate_per_group_memory(int mpi_rank, int node_size,
                                      const mythread_decomp *dc) {
  size_t total = 0;
  for (int g = 0; g < dc->n_groups; g++) {
    const mythread_tile *t = &dc->group_tiles[g];
    size_t pw, ph, plane_b, planes, halo_y, halo_x, mpi_bufs, gf_total;
    if (!add_size_checked((size_t)t->nx, 2u*(size_t)dc->halo, &pw)) return 0;
    if (!add_size_checked((size_t)t->ny, 2u*(size_t)dc->halo, &ph)) return 0;
    if (!multiply_size_checked(ph, pw, &plane_b)) return 0;
    if (!multiply_size_checked(plane_b, sizeof(double)*3, &planes)) return 0;
    if (!multiply_size_checked((size_t)t->nx, sizeof(double)*2, &halo_y)) return 0;
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
  if (mpi_rank == 0)
    printf("Estimated rank memory: %.3f GiB (%d groups)\n", (double)total/(1024.*1024.*1024.), dc->n_groups);
  if (node_size < 1) node_size = 1;
  size_t nb;
  if (!multiply_size_checked(total, (size_t)node_size, &nb)) return 0;
  if (mpi_rank == 0)
    printf("Estimated node memory: %.3f GiB (%d ranks/node)\n", (double)nb/(1024.*1024.*1024.), node_size);
#if defined(__linux__)
  {
    struct sysinfo info;
    if (sysinfo(&info) == 0) {
      unsigned long long vis =
        (unsigned long long)info.totalram * (unsigned long long)info.mem_unit;
      if (mpi_rank == 0)
        printf("Visible node memory : %.3f GiB\n",
               (double)vis / (1024.0 * 1024.0 * 1024.0));
      if ((unsigned long long)nb > vis) {
        if (mpi_rank == 0)
          fprintf(stderr, "[Error] estimated memory exceeds visible memory\n");
        return 0;
      }
    }
  }
#endif
  return 1;
}

/* ── MPI 拓扑 ── */
static void setup_process_domain(int mpi_rank, int mpi_size) {
  int px, py;
  /* 复用 mythread_decomp 的 2D 网格选择: 创建临时对象仅用于 choose_2d_grid */
  {
    mythread_decomp *tmp = mythread_decomp_create(NX, NY, HALO, mpi_size, 0,
                                                   MYTHREAD_DECOMP_XY_2D);
    px = tmp->gx; py = tmp->gy;
    mythread_decomp_free(tmp);
  }

  g_sim.proc_px = px; g_sim.proc_py = py;
  g_sim.proc_x  = mpi_rank % px;
  g_sim.proc_y  = mpi_rank / px;

  /* 切分进程域 */
  {
    int t_nx = NX / px, r_nx = NX % px;
    int t_ny = NY / py, r_ny = NY % py;
    int off_x = g_sim.proc_x * t_nx + (g_sim.proc_x < r_nx ? g_sim.proc_x : r_nx);
    int off_y = g_sim.proc_y * t_ny + (g_sim.proc_y < r_ny ? g_sim.proc_y : r_ny);
    g_sim.local_x_begin = off_x;
    g_sim.local_x_end   = off_x + t_nx + (g_sim.proc_x < r_nx ? 1 : 0);
    g_sim.local_y_begin = off_y;
    g_sim.local_y_end   = off_y + t_ny + (g_sim.proc_y < r_ny ? 1 : 0);
  }

  g_sim.local_nx = g_sim.local_x_end - g_sim.local_x_begin;
  g_sim.local_ny = g_sim.local_y_end - g_sim.local_y_begin;

  g_sim.neighbor_left  = (g_sim.proc_x > 0)      ? (mpi_rank - 1)   : -1;
  g_sim.neighbor_right = (g_sim.proc_x + 1 < px)  ? (mpi_rank + 1)   : -1;
  g_sim.neighbor_down  = (g_sim.proc_y > 0)       ? (mpi_rank - px)  : -1;
  g_sim.neighbor_up    = (g_sim.proc_y + 1 < py)  ? (mpi_rank + px)  : -1;

  /* 填充 MPI 上下文 */
  g_mpi_ctx.mpirank_up    = g_sim.neighbor_up;
  g_mpi_ctx.mpirank_down  = g_sim.neighbor_down;
  g_mpi_ctx.mpirank_left  = g_sim.neighbor_left;
  g_mpi_ctx.mpirank_right = g_sim.neighbor_right;
  g_mpi_ctx.mpi_tag_base  = 100;
}

static double initial_condition_value(int global_y, int global_x) {
  double cx = 0.5 * LX, cy = 0.5 * LY;
  double sigma = 0.06 * ((LX < LY) ? LX : LY);
  double dx = global_x * DX - cx, dy = global_y * DY - cy;
  return A * exp(-(dx * dx + dy * dy) / (2.0 * sigma * sigma));
}

/* ── 初始化 ── */
static void init_simulation(int mpi_rank, int mpi_size) {
  memset(&g_sim, 0, sizeof(g_sim));
  g_sim.mpi_rank = mpi_rank;
  g_sim.mpi_size = mpi_size;
  setup_process_domain(mpi_rank, mpi_size);
  /* 波场平面由各组 GMT 独立分配，此处不分配 */
}

static void free_simulation(void) {
  if (g_gfields) {
    for (int g = 0; g < g_decomp->n_groups; g++)
      group_free_field(&g_gfields[g]);
    free(g_gfields);
    g_gfields = NULL;
  }
  mythread_decomp_free((mythread_decomp*)g_decomp);
  g_decomp = NULL;
  memset(&g_sim, 0, sizeof(g_sim));
}

/* ── Dirichlet 边界 ── */
static void enforce_dirichlet_boundaries(GroupField *gf, int gid) {
  const mythread_tile *tile = &g_decomp->group_tiles[gid];
  int stride = gf->stride;
  int height = gf->ny_padded;

  /* 全局 X 边界 */
  if (g_sim.local_x_begin == 0) {
    int x_col = HALO;
    for (int y = 0; y < height; y++)
      gf->u_curr[GFIDX(gf, y, x_col)] = 0.0;
  }
  if (g_sim.local_x_end == NX) {
    int x_col = HALO + tile->nx - 1;
    for (int y = 0; y < height; y++)
      gf->u_curr[GFIDX(gf, y, x_col)] = 0.0;
  }

  /* MPI X halo — 仅当 MPI 无邻居 且 组在域边界时归零 */
  if (g_sim.neighbor_left < 0 &&
      mythread_decomp_is_domain_boundary(g_decomp, gid, MYTHREAD_NEIGHBOR_LEFT)) {
    for (int y = 0; y < height; y++)
      gf->u_curr[GFIDX(gf, y, 0)] = 0.0;
  }
  if (g_sim.neighbor_right < 0 &&
      mythread_decomp_is_domain_boundary(g_decomp, gid, MYTHREAD_NEIGHBOR_RIGHT)) {
    int x_halo = HALO + tile->nx;
    for (int y = 0; y < height; y++)
      gf->u_curr[GFIDX(gf, y, x_halo)] = 0.0;
  }

  /* MPI Y halo — 仅当 MPI 无邻居 且 组在域边界时归零 */
  if (g_sim.neighbor_down < 0 &&
      mythread_decomp_is_domain_boundary(g_decomp, gid, MYTHREAD_NEIGHBOR_DOWN)) {
    memset(&gf->u_curr[0], 0, (size_t)stride * sizeof(double));
  }
  if (g_sim.neighbor_up < 0 &&
      mythread_decomp_is_domain_boundary(g_decomp, gid, MYTHREAD_NEIGHBOR_UP)) {
    memset(&gf->u_curr[(tile->ny + HALO) * stride], 0,
           (size_t)stride * sizeof(double));
  }
}

static void zero_physical_y_boundaries(GroupField *gf, int gid) {
  const mythread_tile *tile = &g_decomp->group_tiles[gid];
  if (g_sim.local_y_begin == 0 &&
      mythread_decomp_is_domain_boundary(g_decomp, gid, MYTHREAD_NEIGHBOR_DOWN)) {
    memset(&gf->u_curr[GFIDX(gf, HALO, 0)], 0,
           (size_t)gf->stride * sizeof(double));
  }
  if (g_sim.local_y_end == NY &&
      mythread_decomp_is_domain_boundary(g_decomp, gid, MYTHREAD_NEIGHBOR_UP)) {
    memset(&gf->u_curr[GFIDX(gf, tile->ny, 0)], 0,
           (size_t)gf->stride * sizeof(double));
  }
}

/* ── 初始化波场 ── */
static void init_field_block(GroupField *gf, int gid,
                              int y_begin, int y_end,
                              int x_begin, int x_end) {
  for (int y = y_begin; y < y_end; y++) {
    int global_y = global_y_from_local(gid, y);
    for (int x = x_begin; x < x_end; x++) {
      int global_x = global_x_from_local(gid, x);
      size_t p = GFIDX(gf, y, x);
      double value = 0.0;
      if (global_x != 0 && global_x != NX - 1 &&
          global_y != 0 && global_y != NY - 1) {
        value = initial_condition_value(global_y, global_x);
      }
      gf->u_curr[p] = value;
      gf->u_prev[p] = value;
      gf->u_next[p] = 0.0;
    }
  }
}

/* ── 计算核心 ── */
static void compute_block_region(GroupField *gf, int gid,
                                  int y_begin, int y_end,
                                  int x_begin, int x_end) {
  for (int y = y_begin; y < y_end; y++) {
    int global_y = global_y_from_local(gid, y);
    if (global_y == 0 || global_y == NY - 1) continue;

    for (int x = x_begin; x < x_end; x++) {
      int global_x = global_x_from_local(gid, x);
      if (global_x == 0 || global_x == NX - 1) continue;

      double u_ij = gf->u_curr[GFIDX(gf, y, x)];
      double d2x = (gf->u_curr[GFIDX(gf, y, x - 1)] -
                    2.0 * u_ij +
                    gf->u_curr[GFIDX(gf, y, x + 1)]) / (DX * DX);
      double d2y = (gf->u_curr[GFIDX(gf, y - 1, x)] -
                    2.0 * u_ij +
                    gf->u_curr[GFIDX(gf, y + 1, x)]) / (DY * DY);
      gf->u_next[GFIDX(gf, y, x)] =
        2.0 * u_ij - gf->u_prev[GFIDX(gf, y, x)] + C0 * C0 * DT2 * (d2x + d2y);
    }
  }
}

static void compute_interior_block(GroupField *gf, int gid,
                                    int y_begin, int y_end,
                                    int x_begin, int x_end) {
  const mythread_tile *tile = &g_decomp->group_tiles[gid];
  int ny = tile->ny, nx = tile->nx;

  /* 有 MPI 邻居的方向：跳过依赖 halo 数据的行/列 */
  if (g_sim.neighbor_down >= 0 &&
      mythread_decomp_is_domain_boundary(g_decomp, gid, MYTHREAD_NEIGHBOR_DOWN) &&
      y_begin < HALO + 1)
    y_begin = HALO + 1;
  if (g_sim.neighbor_up >= 0 &&
      mythread_decomp_is_domain_boundary(g_decomp, gid, MYTHREAD_NEIGHBOR_UP) &&
      y_end > ny)
    y_end = ny;
  if (g_sim.neighbor_left >= 0 &&
      mythread_decomp_is_domain_boundary(g_decomp, gid, MYTHREAD_NEIGHBOR_LEFT) &&
      x_begin < HALO + 1)
    x_begin = HALO + 1;
  if (g_sim.neighbor_right >= 0 &&
      mythread_decomp_is_domain_boundary(g_decomp, gid, MYTHREAD_NEIGHBOR_RIGHT) &&
      x_end > nx)
    x_end = nx;

  if (y_begin < y_end && x_begin < x_end)
    compute_block_region(gf, gid, y_begin, y_end, x_begin, x_end);
}

static void compute_boundary_block(GroupField *gf, int gid,
                                    int y_begin, int y_end,
                                    int x_begin, int x_end) {
  const mythread_tile *tile = &g_decomp->group_tiles[gid];
  int ny = tile->ny, nx = tile->nx;
  int lower_row = HALO, upper_row = ny;
  int left_col  = HALO, right_col = nx;

  /* clamp 到有效范围 */
  if (x_begin < HALO) x_begin = HALO;
  if (x_end > HALO + nx) x_end = HALO + nx;
  if (y_begin < HALO) y_begin = HALO;
  if (y_end > HALO + ny) y_end = HALO + ny;

  /* 域下边界：有 MPI down 邻居的行 */
  if (g_sim.neighbor_down >= 0 &&
      mythread_decomp_is_domain_boundary(g_decomp, gid, MYTHREAD_NEIGHBOR_DOWN) &&
      y_begin <= lower_row && lower_row < y_end) {
    compute_block_region(gf, gid, lower_row, lower_row + 1, x_begin, x_end);
  }
  /* 域上边界：有 MPI up 邻居的行 */
  if (g_sim.neighbor_up >= 0 &&
      mythread_decomp_is_domain_boundary(g_decomp, gid, MYTHREAD_NEIGHBOR_UP) &&
      y_begin <= upper_row && upper_row < y_end) {
    compute_block_region(gf, gid, upper_row, upper_row + 1, x_begin, x_end);
  }
  /* X 方向域边界 */
  if (g_sim.neighbor_left >= 0 &&
      mythread_decomp_is_domain_boundary(g_decomp, gid, MYTHREAD_NEIGHBOR_LEFT) &&
      x_begin <= left_col && left_col < x_end) {
    compute_block_region(gf, gid, y_begin, y_end, left_col, left_col + 1);
  }
  if (g_sim.neighbor_right >= 0 &&
      mythread_decomp_is_domain_boundary(g_decomp, gid, MYTHREAD_NEIGHBOR_RIGHT) &&
      x_begin <= right_col && right_col < x_end) {
    compute_block_region(gf, gid, y_begin, y_end, right_col, right_col + 1);
  }
}

static double compute_energy_block(GroupField *gf, int gid,
                                    int y_begin, int y_end,
                                    int x_begin, int x_end,
                                    double *l2_out, double *max_out,
                                    int *max_x, int *max_y) {
  if (x_begin >= x_end || y_begin >= y_end) {
    if (l2_out) *l2_out=0.0; if (max_out) *max_out=0.0; return 0.0;
  }

  const mythread_tile *tile = &g_decomp->group_tiles[gid];
  double kinetic = 0.0, potential_x = 0.0, potential_y = 0.0, l2=0.0, ma=0.0;
  double cell_area = DX * DY;
  int mx=0, my=0;

  int x_edge_begin = (x_begin == HALO) ? (HALO - 1) : x_begin;
  int x_edge_end = x_end;
  if (x_edge_end > HALO + tile->nx) x_edge_end = HALO + tile->nx;

  for (int y = y_begin; y < y_end; y++) {
    int global_y = global_y_from_local(gid, y);
    for (int x = x_begin; x < x_end; x++) {
      double u = gf->u_curr[GFIDX(gf, y, x)]; l2 += u*u;
      double au = fabs(u); if (au > ma) { ma=au; mx=global_x_from_local(gid,x); my=global_y; }
    }
    if (global_y > 0 && global_y < NY - 1) {
      for (int x = x_begin; x < x_end; x++) {
        double ut = (gf->u_curr[GFIDX(gf, y, x)] -
                     gf->u_prev[GFIDX(gf, y, x)]) / DT;
        kinetic += ut * ut;
      }
    }
    for (int x = x_edge_begin; x < x_edge_end; x++) {
      double du_curr = (gf->u_curr[GFIDX(gf, y, x + 1)] -
                        gf->u_curr[GFIDX(gf, y, x)]) / DX;
      double du_prev = (gf->u_prev[GFIDX(gf, y, x + 1)] -
                        gf->u_prev[GFIDX(gf, y, x)]) / DX;
      potential_x += du_curr * du_prev;
    }
    if (global_y < NY - 1) {
      for (int x = x_begin; x < x_end; x++) {
        double du_curr = (gf->u_curr[GFIDX(gf, y + 1, x)] -
                          gf->u_curr[GFIDX(gf, y, x)]) / DY;
        double du_prev = (gf->u_prev[GFIDX(gf, y + 1, x)] -
                          gf->u_prev[GFIDX(gf, y, x)]) / DY;
        potential_y += du_curr * du_prev;
      }
    }
  }
  if (l2_out) *l2_out=l2; if (max_out) *max_out=ma; if (max_x) *max_x=mx; if (max_y) *max_y=my;
  return 0.5 * (kinetic + C0 * C0 * (potential_x + potential_y)) * cell_area;
}

/* ── 能量汇总 ── */
static double accumulate_worker_energy(void) {
  double e = 0.0;
  if (uses_group_threads()) {
    for (int g = 0; g < g_decomp->n_groups; g++)
      e += g_gfields[g].group_energy;
  } else {
    for (int t = 0; t < md.Nthreads - 1; t++) {
      ThreadTask *task = (ThreadTask*)md.threads[t].td;
      e += task->partial_energy;
    }
  }
  return e;
}

static double reduce_group_worker_energy(int gid) {
  double e = 0.0;
  threadGroup *pg = md.grps[gid];
  for (int t = 1; t < pg->Nthreads; t++) {
    ThreadTask *task = (ThreadTask*)pg->threads[t].td;
    e += task->partial_energy;
  }
  return e;
}

static double reduce_group_worker_l2(int gid) {
  double s = 0.0;
  threadGroup *pg = md.grps[gid];
  for (int t = 1; t < pg->Nthreads; t++) {
    ThreadTask *task = (ThreadTask*)pg->threads[t].td;
    s += task->partial_l2;
  }
  return s;
}

static void reduce_group_worker_max(int gid, double *max_out, int *gx, int *gy) {
  double ma = 0.0; int mx = 0, my = 0;
  threadGroup *pg = md.grps[gid];
  for (int t = 1; t < pg->Nthreads; t++) {
    ThreadTask *task = (ThreadTask*)pg->threads[t].td;
    if (task->partial_max > ma) { ma = task->partial_max; mx = task->partial_max_gx; my = task->partial_max_gy; }
  }
  *max_out = ma; *gx = mx; *gy = my;
}

static double local_l2_sum(void) {
  double s = 0.0;
  if (uses_group_threads()) {
    for (int g = 0; g < g_decomp->n_groups; g++) s += g_group_l2[g];
  } else {
    for (int t = 0; t < md.Nthreads - 1; t++)
      s += ((ThreadTask*)md.threads[t].td)->partial_l2;
  }
  return s;
}
static void reduce_max_amp(double *amp, int *gx, int *gy) {
  *amp=0.0; *gx=*gy=0;
  if (uses_group_threads()) {
    for (int g=0; g<g_decomp->n_groups; g++)
      if (g_group_max[g] > *amp) { *amp=g_group_max[g]; *gx=g_group_max_x[g]; *gy=g_group_max_y[g]; }
  } else {
    for (int t = 0; t < md.Nthreads - 1; t++) {
      ThreadTask *task = (ThreadTask*)md.threads[t].td;
      if (task->partial_max > *amp) { *amp = task->partial_max; *gx = task->partial_max_gx; *gy = task->partial_max_gy; }
    }
  }
}
static void reset_group_metrics(int gid) {
  g_gfields[gid].group_energy=0.0; g_group_l2[gid]=0.0;
  g_group_max[gid]=0.0; g_group_max_x[gid]=g_group_max_y[gid]=0;
}

/* ── 配置加载 ── */
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
  cfg_GROUP_DECOMP = hw ? mythread_cfg_get_int(hw, "", "GROUP_DECOMP", cfg_GROUP_DECOMP) : cfg_GROUP_DECOMP;
  cfg_NCorePClu  = hw ? mythread_cfg_get_int(hw, "", "NCorePClu",  cfg_NCorePClu)  : cfg_NCorePClu;
  cfg_NCluPNode  = hw ? mythread_cfg_get_int(hw, "", "NCluPNode",  cfg_NCluPNode)  : cfg_NCluPNode;
  cfg_NCorePGrp  = hw ? mythread_cfg_get_int(hw, "", "NCorePGrp",  cfg_NCorePGrp)  : cfg_NCorePGrp;
  cfg_ManageCoreId = hw ? mythread_cfg_get_int(hw, "", "ManageCoreId", cfg_ManageCoreId) : cfg_ManageCoreId;

  mythread_cfg_free(cs);
  mythread_cfg_free(hw);

  /* 环境变量最高优先级 */

  cfg_compute_derived();
}

/* ── 线程任务分配 ── */
/* ── 负载不均衡配置（环境变量）── */
static int g_imbalance_pct = 0;
static int g_imbalance_worker = -1;

static void parse_imbalance_config(void) {
  /* 已由 cfg_load_from_files 通过 mythread_env_get_int 处理，
     此处保留兼容性（直接读环境变量作为 fallback） */
  g_imbalance_pct    = mythread_env_get_int("WAVE_IMBALANCE_PCT",    g_imbalance_pct);
  g_imbalance_worker = mythread_env_get_int("WAVE_IMBALANCE_WORKER", g_imbalance_worker);
}

/*
 * 对一组 worker 施加 Y 方向的不均衡。
 * 选定的 heavy worker 多承担 WAVE_IMBALANCE_PCT% 的行数，
 * 其余 worker 均分剩余行。
 */
static void apply_imbalance_to_group(ThreadTask **workers, int nw) {
  if (g_imbalance_pct <= 0 || nw < 2) return;

  int big_id = (g_imbalance_worker >= 0 && g_imbalance_worker < nw)
                 ? g_imbalance_worker : (nw - 1);

  int total_rows = 0;
  for (int w = 0; w < nw; w++)
    total_rows += workers[w]->y_end - workers[w]->y_begin;

  int avg = total_rows / nw;
  int big_rows = avg + (avg * g_imbalance_pct) / 100;
  if (big_rows >= total_rows) big_rows = total_rows - (nw - 1);
  if (big_rows < 1) big_rows = 1;

  int remaining = total_rows - big_rows;
  int small_rows = remaining / (nw - 1);
  int small_rem  = remaining % (nw - 1);

  int y_cur = workers[0]->y_begin;
  for (int w = 0; w < nw; w++) {
    int rows;
    if (w == big_id) {
      rows = big_rows;
    } else {
      rows = small_rows + (small_rem > 0 ? 1 : 0);
      if (small_rem > 0) small_rem--;
    }
    workers[w]->y_begin = y_cur;
    workers[w]->y_end   = y_cur + rows;
    y_cur += rows;
  }
  workers[nw - 1]->y_end = workers[0]->y_begin + total_rows;

  if (mpi_id == 0) {
    printf("[Imbalance] pct=%d heavy=%d total=%d heavy_rows=%d other_rows=%d\n",
           g_imbalance_pct, big_id, total_rows, big_rows, small_rows);
    for (int w = 0; w < nw; w++)
      printf("  worker %d: y=[%d,%d) rows=%d\n",
             w, workers[w]->y_begin, workers[w]->y_end,
             workers[w]->y_end - workers[w]->y_begin);
  }
}

static void setup_group_thread_tasks(void) {
  parse_imbalance_config();

  for (int g = 0; g < g_decomp->n_groups; g++) {
    threadGroup *pg = md.grps[g];
    const mythread_tile *gtile = &g_decomp->group_tiles[g];
    int nw = pg->Nthreads - 1;
    ThreadTask **wtasks = (nw > 0)
      ? (ThreadTask**)alloca((size_t)nw * sizeof(ThreadTask*)) : NULL;

    for (int t = 0; t < pg->Nthreads; t++) {
      ThreadTask *task = (ThreadTask*)pg->threads[t].td;
      task->gid = g;
      task->tid = t;
      task->gf  = &g_gfields[g];

      if (t == 0 || nw <= 0) {
        task->x_begin = gtile->x_begin;
        task->x_end   = gtile->x_begin;
        task->y_begin = gtile->y_begin;
        task->y_end   = gtile->y_begin;
      } else {
        wtasks[t - 1] = task;
      }
      reset_task_timers(task);
    }

    /*
     * 每组用实际 worker 数重新做 Y 方向切分。
     * 不用 mythread_decomp 的 worker_tiles（它假设统一 worker 数），
     * 因为 group 0 可能少一个 worker（MMT 占槽位）。
     */
    if (nw > 0) {
      int total = gtile->ny;
      int base  = total / nw;
      int rem   = total % nw;
      /* GroupField 内每组 interior 始终从 y=HALO, x=HALO 开始 */
      int y_cur = HALO;
      for (int w = 0; w < nw; w++) {
        int rows = base + (w < rem ? 1 : 0);
        wtasks[w]->x_begin = HALO;
        wtasks[w]->x_end   = HALO + gtile->nx;
        wtasks[w]->y_begin = y_cur;
        wtasks[w]->y_end   = y_cur + rows;
        y_cur += rows;
      }
      apply_imbalance_to_group(wtasks, nw);
    }
  }
}

static void setup_single_group_tasks(void) {
  parse_imbalance_config();

  int nw = g_decomp->n_workers_per_group;
  if (nw <= 0) return;
  const mythread_tile *gtile = &g_decomp->group_tiles[0];

  ThreadTask **wtasks = (ThreadTask**)alloca((size_t)nw * sizeof(ThreadTask*));

  int total = gtile->ny;
  int base  = total / nw;
  int rem   = total % nw;
  int y_cur = HALO;

  for (int t = 0; t < nw; t++) {
    int rows = base + (t < rem ? 1 : 0);
    ThreadTask *task = (ThreadTask*)md.threads[t].td;
    task->gid = 0;
    task->tid = t;
    task->gf  = &g_gfields[0];
    task->x_begin = HALO;
    task->x_end   = HALO + gtile->nx;
    task->y_begin = y_cur;
    task->y_end   = y_cur + rows;
    y_cur += rows;
    wtasks[t] = task;
    reset_task_timers(task);
  }
  apply_imbalance_to_group(wtasks, nw);
}

static void setup_thread_tasks(void) {
  if (uses_group_threads())
    setup_group_thread_tasks();
  else
    setup_single_group_tasks();
}

/* ═══════════════════════════════════════════════════════════════
 *  线程入口
 * ═══════════════════════════════════════════════════════════════ */

static void worker_thread(void) {
  ThreadTask *task = (ThreadTask*)ti->td;
  GroupField *gf = task->gf;
  int gid = task->gid;
  double t0;

  /* init fields */
  t0 = wall_time();
  wait_phase_from_worker(init_fields_state());
  task->t_wait_init += wall_time() - t0;
  t0 = wall_time();
  init_field_block(gf, gid, task->y_begin, task->y_end,
                   task->x_begin, task->x_end);
  task->t_work_init += wall_time() - t0;
  finish_phase_from_worker(init_fields_state());

  /* initial energy */
  t0 = wall_time();
  wait_phase_from_worker(initial_energy_state());
  task->t_wait_energy0 += wall_time() - t0;
  t0 = wall_time();
  double l2_,ma_; int mx_,my_;
  task->partial_energy = compute_energy_block(gf, gid,
      task->y_begin, task->y_end, task->x_begin, task->x_end,
      &l2_, &ma_, &mx_, &my_);
  task->partial_l2 = l2_;
  task->partial_max = ma_; task->partial_max_gx = mx_; task->partial_max_gy = my_;
  task->t_work_energy0 += wall_time() - t0;
  finish_phase_from_worker(initial_energy_state());

  for (int step = 0; step < NT; step++) {
    int cs = compute_phase_state(step);
    int bs = boundary_phase_state(step);
    int es = energy_phase_state(step);

    t0 = wall_time();
    wait_phase_from_worker(cs);
    task->t_wait_compute += wall_time() - t0;
    t0 = wall_time();
    compute_interior_block(gf, gid,
        task->y_begin, task->y_end, task->x_begin, task->x_end);
    task->t_work_compute += wall_time() - t0;
    finish_phase_from_worker(cs);

    t0 = wall_time();
    wait_phase_from_worker(bs);
    task->t_wait_boundary += wall_time() - t0;
    t0 = wall_time();
    compute_boundary_block(gf, gid,
        task->y_begin, task->y_end, task->x_begin, task->x_end);
    task->t_work_boundary += wall_time() - t0;
    finish_phase_from_worker(bs);

    if (should_measure_energy_step(step)) {
      t0 = wall_time();
      wait_phase_from_worker(es);
      task->t_wait_energy += wall_time() - t0;
      t0 = wall_time();
      double l2__,ma__; int mx__,my__;
      task->partial_energy = compute_energy_block(gf, gid,
          task->y_begin, task->y_end, task->x_begin, task->x_end,
          &l2__, &ma__, &mx__, &my__);
      task->partial_l2 = l2__;
      task->partial_max = ma__; task->partial_max_gx = mx__; task->partial_max_gy = my__;
      task->t_work_energy += wall_time() - t0;
      task->energy_steps += 1;
      finish_phase_from_worker(es);
    }
  }
}

static void group_main_thread(void) {
  ThreadTask *task = (ThreadTask*)ti->td;
  int gid = ti->igrp;
  double t0;

  /* ── 分配本组 GroupField（此时已 bindcpu，NUMA 节点正确）── */
  if (group_alloc_field(&g_gfields[gid], gid, g_decomp) != 0) {
    fprintf(stderr, "[Error] MPI=%d GMT gid=%d: group_alloc_field failed\n",
            mpi_id, gid);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  task->gf = &g_gfields[gid];

  /* init fields */
  t0 = wall_time();
  gWaitMain(init_fields_state());
  task->t_wait_init += wall_time() - t0;
  gSetSubs(init_fields_state());
  t0 = wall_time();
  gWaitSubs(init_fields_state());
  task->t_wait_init += wall_time() - t0;
  gSetMain(init_fields_state());

  /* initial energy */
  t0 = wall_time();
  gWaitMain(initial_energy_state());
  task->t_wait_energy0 += wall_time() - t0;
  gSetSubs(initial_energy_state());
  t0 = wall_time();
  gWaitSubs(initial_energy_state());
  task->t_wait_energy0 += wall_time() - t0;
  t0 = wall_time();
  g_gfields[gid].group_energy = reduce_group_worker_energy(gid);
  g_group_l2[gid] = reduce_group_worker_l2(gid);
  reduce_group_worker_max(gid, &g_group_max[gid], &g_group_max_x[gid], &g_group_max_y[gid]);
  task->t_work_energy0 += wall_time() - t0;
  gSetMain(initial_energy_state());

  for (int step = 0; step < NT; step++) {
    int cs = compute_phase_state(step);
    int bs = boundary_phase_state(step);
    int es = energy_phase_state(step);

    t0 = wall_time();
    gWaitMain(cs);
    double t1 = wall_time();
    task->t_wait_compute += t1 - t0;
    gSetSubs(cs);
    t0 = wall_time();
    gWaitSubs(cs);
    task->t_wait_compute += wall_time() - t0;
    gSetMain(cs);

    t0 = wall_time();
    gWaitMain(bs);
    task->t_wait_boundary += wall_time() - t0;
    gSetSubs(bs);
    t0 = wall_time();
    gWaitSubs(bs);
    task->t_wait_boundary += wall_time() - t0;
    gSetMain(bs);

    if (should_measure_energy_step(step)) {
      t0 = wall_time();
      gWaitMain(es);
      task->t_wait_energy += wall_time() - t0;
      gSetSubs(es);
      t0 = wall_time();
      gWaitSubs(es);
      task->t_wait_energy += wall_time() - t0;
      t0 = wall_time();
      g_gfields[gid].group_energy = reduce_group_worker_energy(gid);
      g_group_l2[gid] = reduce_group_worker_l2(gid);
      reduce_group_worker_max(gid, &g_group_max[gid], &g_group_max_x[gid], &g_group_max_y[gid]);
      task->t_work_energy += wall_time() - t0;
      task->energy_steps += 1;
      gSetMain(es);
    }
  }
  printf("Group-Summary%d-%d: %.3f %.3f\n",
         mpi_id, gid, task->t_wait_compute,
         wall_time() - task->t_wait_compute); /* rough work est */
}

/* Dirichlet 边界应用到所有组 */
static void apply_dirichlet_all_groups(void) {
  for (int g = 0; g < g_decomp->n_groups; g++) {
    enforce_dirichlet_boundaries(&g_gfields[g], g);
    zero_physical_y_boundaries(&g_gfields[g], g);
  }
}

/* 复制 u_curr 到 u_prev（初始化后同步 halo） */
static void copy_curr_to_prev_all_groups(void) {
  for (int g = 0; g < g_decomp->n_groups; g++)
    memcpy(g_gfields[g].u_prev, g_gfields[g].u_curr, g_gfields[g].plane_bytes);
}

static void main_thread(void) {
  ThreadTask *task = (ThreadTask*)ti->td;
  double start_time, local_energy, global_energy;
  double t0;

  reset_task_timers(task);

  /* Phase 0: alloc GroupField */
  if (!uses_group_threads()) {
    if (group_alloc_field(&g_gfields[0], 0, g_decomp) != 0) {
      fprintf(stderr, "[Error] MPI=%d: group_alloc_field failed\n", mpi_id);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    for (int t = 0; t < g_decomp->n_workers_per_group; t++) {
      ThreadTask *wt = (ThreadTask*)md.threads[t].td;
      wt->gf = &g_gfields[0];
    }
  }

  /* Wait for GMTs to alloc + workers to init fields */
  t0 = wall_time();
  start_phase_from_main(init_fields_state());
  wait_phase_from_main(init_fields_state());
  task->t_wait_init += wall_time() - t0;

  group_field_link_buffers(g_gfields, g_decomp);

  /* Dirichlet + halo exchange for u_curr */
  apply_dirichlet_all_groups();
  t0 = wall_time();
  mythread_halo_exchange_intra(g_gfields, g_decomp);
  mythread_halo_exchange_mpi(g_gfields, g_decomp, &g_mpi_ctx,
                              (uintptr_t)MPI_COMM_WORLD);
  apply_dirichlet_all_groups();
  task->t_comm += wall_time() - t0;

  /* sync u_prev halos */
  copy_curr_to_prev_all_groups();

  /* initial energy */
  t0 = wall_time();
  start_phase_from_main(initial_energy_state());
  wait_phase_from_main(initial_energy_state());
  task->t_wait_energy0 += wall_time() - t0;
  t0 = wall_time();
  local_energy = accumulate_worker_energy();
  task->t_work_energy0 += wall_time() - t0;
  t0 = wall_time();
  MPI_Allreduce(&local_energy, &global_energy, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  task->t_allreduce += wall_time() - t0;
  g_sim.initial_energy = global_energy;
  { /* L2: MPI-reduce the sum-of-squares, then compute sqrt */
    double local_l2 = local_l2_sum();
    double global_l2_sum;
    MPI_Allreduce(&local_l2, &global_l2_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double global_l2 = sqrt(global_l2_sum * DX * DY);
    /* max|u|: MPI-reduce position-aware.  Use separate MPI_MAX on double
       (avoids MPI_DOUBLE_INT struct-padding portability issues), then
       broadcast coordinates from the winning rank.  All ranks must
       participate in MPI_Bcast regardless of who holds the max. */
    double local_ma; int lmx, lmy; reduce_max_amp(&local_ma, &lmx, &lmy);
    double global_ma;
    MPI_Allreduce(&local_ma, &global_ma, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    int has_max = (local_ma >= global_ma) ? mpi_id : -1;
    int winner;
    MPI_Allreduce(&has_max, &winner, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    int gx = lmx, gy = lmy;
    MPI_Bcast(&gx, 1, MPI_INT, winner, MPI_COMM_WORLD);
    MPI_Bcast(&gy, 1, MPI_INT, winner, MPI_COMM_WORLD);
    if (mpi_id == 0)
      printf("[Main] Initial: E=%.6f L2=%.6f max|u|=%.6f@(%d,%d)\n",
             g_sim.initial_energy, global_l2, global_ma, gx, gy);
  }

  start_time = MPI_Wtime();
  double prev_time = start_time;

  for (int step = 0; step < NT; step++) {
    MPI_Barrier(MPI_COMM_WORLD);
    int cs = compute_phase_state(step);
    int bs = boundary_phase_state(step);
    int es = energy_phase_state(step);
    int need_energy = should_measure_energy_step(step);

    /* Dirichlet + halo */
    apply_dirichlet_all_groups();
    t0 = wall_time();
    mythread_halo_exchange_intra(g_gfields, g_decomp);
    mythread_halo_exchange_mpi(g_gfields, g_decomp, &g_mpi_ctx,
                                (uintptr_t)MPI_COMM_WORLD);
    apply_dirichlet_all_groups();
    task->t_comm += wall_time() - t0;

    /* ── compute interior ── */
    t0 = wall_time();
    start_phase_from_main(cs);
    wait_phase_from_main(cs);
    task->t_wait_compute += wall_time() - t0;

    /* ── compute boundary ── */
    t0 = wall_time();
    start_phase_from_main(bs);
    wait_phase_from_main(bs);
    task->t_wait_boundary += wall_time() - t0;

    /* ── swap fields（每组独立）── */
    for (int g = 0; g < g_decomp->n_groups; g++)
      group_field_swap(&g_gfields[g]);

    if (need_energy) {
      /* energy needs fresh halo */
      apply_dirichlet_all_groups();
      t0 = wall_time();
      mythread_halo_exchange_intra(g_gfields, g_decomp);
      mythread_halo_exchange_mpi(g_gfields, g_decomp, &g_mpi_ctx,
                                  (uintptr_t)MPI_COMM_WORLD);
      apply_dirichlet_all_groups();
      task->t_comm += wall_time() - t0;

      t0 = wall_time();
      start_phase_from_main(es);
      wait_phase_from_main(es);
      task->t_wait_energy += wall_time() - t0;
      t0 = wall_time();
      local_energy = accumulate_worker_energy();
      task->t_work_energy += wall_time() - t0;
      t0 = wall_time();
      MPI_Allreduce(&local_energy, &global_energy, 1, MPI_DOUBLE, MPI_SUM,
                    MPI_COMM_WORLD);
      task->t_allreduce += wall_time() - t0;
      task->energy_steps += 1;

      if (mpi_id == 0) {
        double cur_time = MPI_Wtime();
        double compute_time = cur_time - prev_time;
        prev_time = cur_time;
        /* L2: MPI-reduce */
        double local_l2 = local_l2_sum();
        double global_l2_sum;
        MPI_Allreduce(&local_l2, &global_l2_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        double global_l2 = sqrt(global_l2_sum * DX * DY);
        /* max|u|: MPI-reduce position-aware (same safe pattern as initial) */
        double local_ma; int lmx, lmy; reduce_max_amp(&local_ma, &lmx, &lmy);
        double global_ma;
        MPI_Allreduce(&local_ma, &global_ma, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        int has_max2 = (local_ma >= global_ma) ? mpi_id : -1;
        int winner2;
        MPI_Allreduce(&has_max2, &winner2, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        int gx = lmx, gy = lmy;
        MPI_Bcast(&gx, 1, MPI_INT, winner2, MPI_COMM_WORLD);
        MPI_Bcast(&gy, 1, MPI_INT, winner2, MPI_COMM_WORLD);
        printf("[Main] Step %4d/%d, time %.3f,  E=%.6f L2=%.6f max|u|=%.6f@(%d,%d)\n",
               step + 1, NT, compute_time, global_energy, global_l2, global_ma, gx, gy);
      }
    }
  }

  double elapsed = MPI_Wtime() - start_time;
  if (mpi_id == 0) {
    double points = (double)NX * (double)NY * (double)NT;
    printf("[Main] Simulation completed in %.3f seconds\n", elapsed);
    printf("[Main] comm time: %.3f, compute wait: %.3f, boundary wait: %.3f\n",
           task->t_comm, task->t_wait_compute, task->t_wait_boundary);
    printf("[Main] Throughput: %.2f Mpoint-updates/s\n",
           points / elapsed / 1.0e6);
  }
  printf("Main-summary-%d: %d %d %.3f %.3f %.3f %.3f\n",
         mpi_id, g_sim.local_nx, g_sim.local_ny + HALO,
         elapsed, task->t_comm, task->t_wait_compute, task->t_wait_boundary);
}

void thread_run(void) {
  ThreadTask *task = (ThreadTask*)ti->td;
#ifdef DEBUG
  printf("[info] mpi %d gid %d, pid %d, c %d, getcore %d\n",
         mpi_id, ti->igrp, ti->ind, ti->indg, getcpuid());
#endif
  if (task) {
    task->gid = ti->igrp;
    task->tid = ti->ind;
    task->cpu_id = getcpuid();
  }

  if (ti->igrp < 0) {
    main_thread();
  } else if (!uses_group_threads()) {
    worker_thread();
  } else if (ti->ind == 0) {
    group_main_thread();
  } else {
    worker_thread();
  }
}

/* ── 计时报告 ── */
static const char *task_role_name(const ThreadTask *task) {
  if (!task) return "null";
  if (task->gid < 0) return "main";
  if (uses_group_threads() && task->tid == 0) return "group-main";
  return "worker";
}

static double task_total_time(const ThreadTask *task) {
  if (!task) return 0.0;
  return task->t_wait_init + task->t_work_init +
         task->t_wait_energy0 + task->t_work_energy0 +
         task->t_wait_compute + task->t_work_compute +
         task->t_wait_boundary + task->t_work_boundary +
         task->t_wait_energy + task->t_work_energy +
         task->t_comm + task->t_allreduce;
}

static void print_one_task_timing(int mpi_rank, int mpi_size, int node_size,
                                  const ThreadTask *task) {
  double total = task_total_time(task);
  (void)mpi_rank; (void)mpi_size; (void)node_size; (void)total;
#ifdef DEBUG
  int gyb = -1, gye = -1;
  if (task && task->y_end > task->y_begin && task->y_begin >= HALO) {
    gyb = global_y_from_local(task->gid, task->y_begin);
    gye = global_y_from_local(task->gid, task->y_end);
  }
  printf("[Timing] rank %d/%d gid=%d tid=%d role=%s cpu=%d "
         "x=[%d,%d) y=[%d,%d) global-y=[%d,%d) "
         "wait(i=%.4f e0=%.4f c=%.4f b=%.4f e=%.4f) "
         "work(i=%.4f e0=%.4f c=%.4f b=%.4f e=%.4f) "
         "comm=%.4f allreduce=%.4f steps=%d total=%.4f\n",
         mpi_rank, mpi_size,
         task ? task->gid : -1, task ? task->tid : -1,
         task_role_name(task), task ? task->cpu_id : -1,
         task ? task->x_begin : 0, task ? task->x_end : 0,
         task ? task->y_begin : 0, task ? task->y_end : 0,
         gyb, gye,
         task ? task->t_wait_init : 0.0,
         task ? task->t_wait_energy0 : 0.0,
         task ? task->t_wait_compute : 0.0,
         task ? task->t_wait_boundary : 0.0,
         task ? task->t_wait_energy : 0.0,
         task ? task->t_work_init : 0.0,
         task ? task->t_work_energy0 : 0.0,
         task ? task->t_work_compute : 0.0,
         task ? task->t_work_boundary : 0.0,
         task ? task->t_work_energy : 0.0,
         task ? task->t_comm : 0.0,
         task ? task->t_allreduce : 0.0,
         task ? task->energy_steps : 0, total);
#endif
}

static void print_timing_report(int mpi_rank, int mpi_size, int node_size) {
  for (int r = 0; r < mpi_size; r++) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (mpi_rank == r) {
      print_one_task_timing(mpi_rank, mpi_size, node_size,
                            (ThreadTask*)md.tm.td);
      if (uses_group_threads()) {
        for (int g = 0; g < g_decomp->n_groups; g++) {
          threadGroup *pg = md.grps[g];
          for (int t = 0; t < pg->Nthreads; t++)
            print_one_task_timing(mpi_rank, mpi_size, node_size,
                                  (ThreadTask*)pg->threads[t].td);
        }
      } else {
        for (int t = 0; t < md.Nthreads - 1; t++)
          print_one_task_timing(mpi_rank, mpi_size, node_size,
                                (ThreadTask*)md.threads[t].td);
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
  int mpi_rank, mpi_size;
  int NCorePClu, NCluPNode, NCorePGrp;
  int NThPGrp, NGrpPProc, NProcPNode = 1;
  int ManageCoreId;
  int local_ready = 1, global_ready = 1;
  int err;

  MPI_Init(&argc, &argv);

  /* 加载配置文件（config 目录相对于工作目录，可用环境变量覆盖路径） */
  {
    const char *cp = mythread_env_get("WAVE_CASE_CFG",    "config/case.cfg");
    const char *hp = mythread_env_get("WAVE_HARDWARE_CFG","config/hardware.cfg");
    cfg_load_from_files(cp, hp);
  }

  /* 使用配置值 */
  NCorePClu  = cfg_NCorePClu;
  NCluPNode  = cfg_NCluPNode;
  NCorePGrp  = cfg_NCorePGrp;
  NThPGrp    = cfg_THREADS_PER_GROUP;
  NGrpPProc  = cfg_N_GROUPS;
  ManageCoreId = cfg_ManageCoreId;
  CoreOffset  = cfg_CoreOffset;
  ClustOffset = cfg_ClustOffset;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  MPI_Comm node_comm;
  MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,
                      MPI_INFO_NULL, &node_comm);
  int node_size;
  MPI_Comm_size(node_comm, &node_size);
  NProcPNode = node_size;

  if (NX < 3 || NY < 3) {
    if (mpi_rank == 0)
      fprintf(stderr, "[Error] NX and NY must both be >= 3\n");
    MPI_Finalize(); return 1;
  }
  if ((long long)mpi_size > (long long)NX * (long long)NY) {
    if (mpi_rank == 0)
      fprintf(stderr, "[Error] mpi_size too large\n");
    MPI_Finalize(); return 1;
  }
  if (CFL_SUM2 > 1.0) {
    if (mpi_rank == 0)
      fprintf(stderr, "[Error] unstable CFL=%.6f > 1\n", CFL_SUM2);
    MPI_Finalize(); return 1;
  }

  /* ── MPI 拓扑初始化（不分配波场）── */
  init_simulation(mpi_rank, mpi_size);
  if (g_sim.local_nx <= 0 || g_sim.local_ny <= 0) {
    fprintf(stderr, "[Error] MPI=%d: invalid decomposition\n", mpi_rank);
    local_ready = 0;
  }
  MPI_Allreduce(&local_ready, &global_ready, 1, MPI_INT, MPI_MIN,
                MPI_COMM_WORLD);
  if (!global_ready) { free_simulation(); MPI_Finalize(); return 1; }

  int n_workers = uses_group_threads() ? (NThPGrp - 1) : (NThPGrp - 1);
  g_decomp = mythread_decomp_create(g_sim.local_nx, g_sim.local_ny, HALO,
                                     NGrpPProc, n_workers, cfg_GROUP_DECOMP);
  if (!g_decomp) { fprintf(stderr, "[Error] decomp_create failed\n"); MPI_Finalize(); return 1; }
  if (!validate_per_group_memory(mpi_rank, node_size, g_decomp)) {
    mythread_decomp_free((mythread_decomp*)g_decomp); g_decomp=NULL; MPI_Finalize(); return 1;
  }

  /* ── 分配 GroupField 注册表（MMT 中填充）── */
  g_gfields = (GroupField*)xcalloc((size_t)g_decomp->n_groups,
                                    sizeof(GroupField));
  g_group_l2   = (double*)xcalloc((size_t)g_decomp->n_groups, sizeof(double));
  g_group_max  = (double*)xcalloc((size_t)g_decomp->n_groups, sizeof(double));
  g_group_max_x = (int*)xcalloc((size_t)g_decomp->n_groups, sizeof(int));
  g_group_max_y = (int*)xcalloc((size_t)g_decomp->n_groups, sizeof(int));

  if (mpi_rank == 0) {
    printf("============================================\n");
    printf("  Wave Equation (NUMA-aware, per-group fields)\n");
    printf("============================================\n");
    printf("Global grid       : %d x %d\n", NX, NY);
    printf("Time steps        : %d\n", NT);
    printf("MPI processes     : %d (%d x %d)\n",
           mpi_size, g_sim.proc_px, g_sim.proc_py);
    printf("Thread groups     : %d\n", NGrpPProc);
    printf("Threads/group     : %d\n", cfg_THREADS_PER_GROUP);
    printf("Group decomp      : %s\n",
           g_decomp->policy == MYTHREAD_DECOMP_Y_ONLY ? "Y_ONLY" : "XY_2D");
    printf("Group grid        : %d x %d\n", g_decomp->gx, g_decomp->gy);
    printf("DT=%.4f DX=%.4f CFL: x=%.4f y=%.4f sum2=%.6f\n",
           cfg_DT, DX, CFL_X, CFL_Y, CFL_SUM2);
    printf("============================================\n");
  }

  MPI_Barrier(MPI_COMM_WORLD);
  printf("[Domain] rank %d/%d node_size=%d proc=(%d,%d)/(%d,%d) "
         "local=(nx=%d,ny=%d) neighbors(L=%d R=%d D=%d U=%d)\n",
         mpi_rank, mpi_size, node_size,
         g_sim.proc_x, g_sim.proc_y, g_sim.proc_px, g_sim.proc_py,
         g_sim.local_nx, g_sim.local_ny,
         g_sim.neighbor_left, g_sim.neighbor_right,
         g_sim.neighbor_down, g_sim.neighbor_up);
  MPI_Barrier(MPI_COMM_WORLD);

  /* ── 初始化 mythread ── */
  err = InitThreads(mpi_rank, NCorePClu, NCluPNode, NCorePGrp,
                    NThPGrp, NGrpPProc, NProcPNode, &ManageCoreId);
  if (err != 0) {
    fprintf(stderr, "[Error] MPI=%d: InitThreads failed: %d\n", mpi_rank, err);
    free_simulation();
    MPI_Finalize();
    return 1;
  }

  setup_thread_tasks();

  StartThreads(thread_run);
  thread_run();
  EndThreads();

  print_timing_report(mpi_rank, mpi_size, node_size);

  free_simulation();
  MPI_Finalize();
  return 0;
}
