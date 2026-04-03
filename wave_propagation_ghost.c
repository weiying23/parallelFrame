#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mythread/mythread.h"

/*
 * 重新设计的二维波动方程示例。
 *
 * 并行布局：
 * - MPI 沿 Y 方向切分全局网格，每步交换一层 halo 行；
 * - 每个 MPI 进程内部再由 mythread 线程组切分本地 Y 行；
 * - 每个组内部的线程继续切分 X 方向内部列。
 *
 * 与旧版示例相比，这里每个 MPI 进程只保留一份共享局部网格，
 * 线程只负责计算不重叠的 tile，不再做线程级 ghost 交换。
 */

#define NX 4096
#define NY 4096
#define NT 400

#define DT 0.001
#define C0 0.1

/*
 * 增大计算量有两种常见方式：
 * 1. USE_FIXED_DOMAIN=1：固定物理区域，增大 NX/NY 表示网格加密、分辨率提高。
 * 2. USE_FIXED_DOMAIN=0：固定 DX/DY，增大 NX/NY 表示物理区域扩大、总网格点增加。
 */
#define USE_FIXED_DOMAIN 1

#if USE_FIXED_DOMAIN
#define LX 1.0
#define LY 1.0
#define DX (LX / (double)(NX - 1))
#define DY (LY / (double)(NY - 1))
#else
#define DX 0.01
#define DY 0.01
#define LX ((NX - 1) * DX)
#define LY ((NY - 1) * DY)
#endif

#define HALO 1
#define N_GROUPS 2
#define N_WORKERS 8
#define THREADS_PER_GROUP (N_WORKERS + 1)

#define DT2 (DT * DT)
#define CFL_X (C0 * DT / DX)
#define CFL_Y (C0 * DT / DY)
#define CFL_SUM2 (CFL_X * CFL_X + CFL_Y * CFL_Y)

typedef struct {
  double *u_prev;
  double *u_curr;
  double *u_next;
  int local_y_begin;
  int local_y_end;
  int local_ny;
  int mpi_rank;
  int mpi_size;
  int neighbor_up;
  int neighbor_down;
  double initial_energy;
} SimulationData;

typedef struct {
  int gid;
  int tid;
  int x_begin;
  int x_end;
  int y_begin;
  int y_end;
} ThreadTask;

static SimulationData g_sim = {0};
static double *g_send_up = NULL;
static double *g_recv_up = NULL;
static double *g_send_down = NULL;
static double *g_recv_down = NULL;

int _gettdsize_() {
  return (int)sizeof(ThreadTask);
}

int _getgdsize_() {
  return 0;
}

static inline int idx(int y, int x) {
  return y * NX + x;
}

static inline int global_y_from_local(int local_y) {
  return g_sim.local_y_begin + (local_y - HALO);
}

static void split_range(int begin,int end,int parts,int index,int *sub_begin,int *sub_end) {
  int total = end - begin;
  int base = total / parts;
  int rem = total % parts;
  int extra = (index < rem) ? 1 : 0;
  int offset = index * base + ((index < rem) ? index : rem);

  *sub_begin = begin + offset;
  *sub_end = *sub_begin + base + extra;
}

static void *xcalloc(size_t count,size_t size) {
  void *p = calloc(count,size);
  if (!p) {
    fprintf(stderr,"[Error] MPI=%d: allocation failed (%zu x %zu)\n", g_sim.mpi_rank, count, size);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  return p;
}

static void setup_process_domain(int mpi_rank,int mpi_size) {
  split_range(0, NY, mpi_size, mpi_rank, &g_sim.local_y_begin, &g_sim.local_y_end);
  g_sim.local_ny = g_sim.local_y_end - g_sim.local_y_begin;
  g_sim.neighbor_down = (g_sim.local_y_begin > 0) ? mpi_rank - 1 : -1;
  g_sim.neighbor_up = (g_sim.local_y_end < NY) ? mpi_rank + 1 : -1;
}

static void seed_initial_condition(void) {
  double cx = 0.5 * LX;
  double cy = 0.5 * LY;
  double sigma = 0.06 * ((LX < LY) ? LX : LY);

  for (int ly = 0; ly < g_sim.local_ny; ly++) {
    int gy = g_sim.local_y_begin + ly;
    double y = gy * DY;

    if (gy == 0 || gy == NY - 1) {
      continue;
    }

    for (int x = 1; x < NX - 1; x++) {
      double dx = x * DX - cx;
      double dy = y - cy;
      double value = exp(-(dx * dx + dy * dy) / (2.0 * sigma * sigma));
      int p = idx(ly + HALO, x);

      g_sim.u_curr[p] = value;
      g_sim.u_prev[p] = value;
    }
  }
}

static void init_simulation(int mpi_rank,int mpi_size) {
  size_t plane_size;

  memset(&g_sim, 0, sizeof(g_sim));
  g_sim.mpi_rank = mpi_rank;
  g_sim.mpi_size = mpi_size;

  setup_process_domain(mpi_rank, mpi_size);
  plane_size = (size_t)(g_sim.local_ny + 2 * HALO) * (size_t)NX;

  g_sim.u_prev = (double*)xcalloc(plane_size, sizeof(double));
  g_sim.u_curr = (double*)xcalloc(plane_size, sizeof(double));
  g_sim.u_next = (double*)xcalloc(plane_size, sizeof(double));

  g_send_up = (double*)xcalloc((size_t)NX, sizeof(double));
  g_recv_up = (double*)xcalloc((size_t)NX, sizeof(double));
  g_send_down = (double*)xcalloc((size_t)NX, sizeof(double));
  g_recv_down = (double*)xcalloc((size_t)NX, sizeof(double));

  seed_initial_condition();
}

static void free_simulation(void) {
  free(g_sim.u_prev);
  free(g_sim.u_curr);
  free(g_sim.u_next);
  free(g_send_up);
  free(g_recv_up);
  free(g_send_down);
  free(g_recv_down);

  memset(&g_sim, 0, sizeof(g_sim));
  g_send_up = NULL;
  g_recv_up = NULL;
  g_send_down = NULL;
  g_recv_down = NULL;
}

static void enforce_dirichlet_boundaries(double *field) {
  for (int y = 0; y < g_sim.local_ny + 2 * HALO; y++) {
    field[idx(y, 0)] = 0.0;
    field[idx(y, NX - 1)] = 0.0;
  }

  if (g_sim.neighbor_down < 0) {
    memset(&field[idx(0, 0)], 0, (size_t)NX * sizeof(double));
  }
  if (g_sim.neighbor_up < 0) {
    memset(&field[idx(g_sim.local_ny + HALO, 0)], 0, (size_t)NX * sizeof(double));
  }
}

static void zero_physical_y_boundaries(double *field) {
  if (g_sim.local_y_begin == 0) {
    memset(&field[idx(HALO, 0)], 0, (size_t)NX * sizeof(double));
  }
  if (g_sim.local_y_end == NY) {
    memset(&field[idx(g_sim.local_ny, 0)], 0, (size_t)NX * sizeof(double));
  }
}

static void exchange_y_halos_for(double *field) {
  MPI_Status status;

  enforce_dirichlet_boundaries(field);
  zero_physical_y_boundaries(field);

  if (g_sim.neighbor_up >= 0) {
    memcpy(g_send_up, &field[idx(g_sim.local_ny, 0)], (size_t)NX * sizeof(double));
    MPI_Sendrecv(
      g_send_up, NX, MPI_DOUBLE, g_sim.neighbor_up, 100,
      g_recv_up, NX, MPI_DOUBLE, g_sim.neighbor_up, 101,
      MPI_COMM_WORLD, &status
    );
    memcpy(&field[idx(g_sim.local_ny + HALO, 0)], g_recv_up, (size_t)NX * sizeof(double));
  } else {
    memset(&field[idx(g_sim.local_ny + HALO, 0)], 0, (size_t)NX * sizeof(double));
  }

  if (g_sim.neighbor_down >= 0) {
    memcpy(g_send_down, &field[idx(HALO, 0)], (size_t)NX * sizeof(double));
    MPI_Sendrecv(
      g_send_down, NX, MPI_DOUBLE, g_sim.neighbor_down, 101,
      g_recv_down, NX, MPI_DOUBLE, g_sim.neighbor_down, 100,
      MPI_COMM_WORLD, &status
    );
    memcpy(&field[idx(0, 0)], g_recv_down, (size_t)NX * sizeof(double));
  } else {
    memset(&field[idx(0, 0)], 0, (size_t)NX * sizeof(double));
  }

  enforce_dirichlet_boundaries(field);
  zero_physical_y_boundaries(field);
}

static void exchange_y_halos(void) {
  exchange_y_halos_for(g_sim.u_curr);
}

static void compute_block(ThreadTask *task) {
  for (int y = task->y_begin; y < task->y_end; y++) {
    int global_y = global_y_from_local(y);

    if (global_y == 0 || global_y == NY - 1) {
      continue;
    }
    for (int x = task->x_begin; x < task->x_end; x++) {
      double u_ij = g_sim.u_curr[idx(y, x)];
      double d2x = (g_sim.u_curr[idx(y, x - 1)] - 2.0 * u_ij + g_sim.u_curr[idx(y, x + 1)]) / (DX * DX);
      double d2y = (g_sim.u_curr[idx(y - 1, x)] - 2.0 * u_ij + g_sim.u_curr[idx(y + 1, x)]) / (DY * DY);
      double next = 2.0 * u_ij - g_sim.u_prev[idx(y, x)] + C0 * C0 * DT2 * (d2x + d2y);

      g_sim.u_next[idx(y, x)] = next;
    }
  }
}

static double compute_total_energy_local(void) {
  double kinetic = 0.0;
  double potential_x = 0.0;
  double potential_y = 0.0;
  double cell_area = DX * DY;

  for (int y = HALO; y < g_sim.local_ny + HALO; y++) {
    int global_y = global_y_from_local(y);

    for (int x = 1; x < NX - 1; x++) {
      double ut;

      if (global_y == 0 || global_y == NY - 1) {
        continue;
      }

      ut = (g_sim.u_curr[idx(y, x)] - g_sim.u_prev[idx(y, x)]) / DT;
      kinetic += ut * ut;
    }

    for (int x = 0; x < NX - 1; x++) {
      double du_curr = (g_sim.u_curr[idx(y, x + 1)] - g_sim.u_curr[idx(y, x)]) / DX;
      double du_prev = (g_sim.u_prev[idx(y, x + 1)] - g_sim.u_prev[idx(y, x)]) / DX;
      potential_x += du_curr * du_prev;
    }

    if (global_y < NY - 1) {
      for (int x = 1; x < NX - 1; x++) {
        double du_curr = (g_sim.u_curr[idx(y + 1, x)] - g_sim.u_curr[idx(y, x)]) / DY;
        double du_prev = (g_sim.u_prev[idx(y + 1, x)] - g_sim.u_prev[idx(y, x)]) / DY;
        potential_y += du_curr * du_prev;
      }
    }
  }

  return 0.5 * (kinetic + C0 * C0 * (potential_x + potential_y)) * cell_area;
}

static void swap_fields(void) {
  double *tmp = g_sim.u_prev;
  g_sim.u_prev = g_sim.u_curr;
  g_sim.u_curr = g_sim.u_next;
  g_sim.u_next = tmp;
}

static void setup_thread_tasks(void) {
  for (int g = 0; g < md.ngrp; g++) {
    threadGroup *pg = md.grps[g];
    int y_begin, y_end;

    split_range(HALO, g_sim.local_ny + HALO, md.ngrp, g, &y_begin, &y_end);

    for (int t = 0; t < pg->Nthreads; t++) {
      THREADINFO *pti = &pg->threads[t];
      ThreadTask *task = (ThreadTask*)pti->td;
      int x_begin, x_end;

      split_range(1, NX - 1, pg->Nthreads, t, &x_begin, &x_end);

      task->gid = g;
      task->tid = t;
      task->x_begin = x_begin;
      task->x_end = x_end;
      task->y_begin = y_begin;
      task->y_end = y_end;

      printf(
        "[Init] MPI=%d Group=%d Thread=%d global-y=[%d,%d) x=[%d,%d)\n",
        mpi_id,
        g,
        t,
        g_sim.local_y_begin + (y_begin - HALO),
        g_sim.local_y_begin + (y_end - HALO),
        x_begin,
        x_end
      );
    }
  }
}

static void worker_thread(void) {
  ThreadTask *task = (ThreadTask*)ti->td;

  for (int step = 0; step < NT; step++) {
    int state = step + 1;

    sWaitGrp(state);
    compute_block(task);
    sSetGrp(state);
  }
}

static void group_main_thread(void) {
  ThreadTask *task = (ThreadTask*)ti->td;

  for (int step = 0; step < NT; step++) {
    int state = step + 1;

    gWaitMain(state);
    gSetSubs(state);
    compute_block(task);
    gWaitSubs(state);
    gSetMain(state);
  }
}

static void main_thread(void) {
  double start_time;
  double local_energy = 0.0;
  double global_energy = 0.0;

  exchange_y_halos_for(g_sim.u_curr);
  exchange_y_halos_for(g_sim.u_prev);
  local_energy = compute_total_energy_local();
  MPI_Allreduce(&local_energy, &global_energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  g_sim.initial_energy = global_energy;

  if (mpi_id == 0) {
    printf("[Main] Initial total energy: %.6f\n", g_sim.initial_energy);
  }

  start_time = MPI_Wtime();

  for (int step = 0; step < NT; step++) {
    int state = step + 1;
    size_t plane_size = (size_t)(g_sim.local_ny + 2 * HALO) * (size_t)NX;

    memset(g_sim.u_next, 0, plane_size * sizeof(double));

    mSetGrps(state);
    mWaitGrps(state);

    swap_fields();
    exchange_y_halos();
    local_energy = compute_total_energy_local();
    MPI_Allreduce(&local_energy, &global_energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if (mpi_id == 0 && (step == 0 || (step + 1) % 100 == 0 || step == NT - 1)) {
      double rel_diff = (g_sim.initial_energy > 0.0)
        ? fabs(global_energy - g_sim.initial_energy) / g_sim.initial_energy
        : 0.0;

      printf(
        "[Main] Step %4d/%d  Energy=%.6f  RelDiff=%.3e\n",
        step + 1,
        NT,
        global_energy,
        rel_diff
      );
    }
  }

  if (mpi_id == 0) {
    double elapsed = MPI_Wtime() - start_time;
    double points = (double)NX * (double)NY * (double)NT;
    printf("[Main] Simulation completed in %.3f seconds\n", elapsed);
    printf("[Main] Throughput: %.2f Mpoint-updates/s\n", points / elapsed / 1.0e6);
  }
}

void thread_run(void) {
  if (ti->igrp < 0) {
    main_thread();
  } else if (ti->ind == 0) {
    group_main_thread();
  } else {
    worker_thread();
  }
}

int main(int argc,char **argv) {
  int mpi_rank, mpi_size;
  int NCorePClu = 38;
  int NCluPNode = 16;
  int NCorePGrp = 8;
  int NThPGrp = THREADS_PER_GROUP;
  int NGrpPProc = N_GROUPS;
  int NProcPNode = 4;
  int ManageCoreId = -2;
  int err;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  if (NX < 3 || NY < 3) {
    if (mpi_rank == 0) {
      fprintf(stderr, "[Error] NX and NY must both be >= 3\n");
    }
    MPI_Finalize();
    return 1;
  }
  if (mpi_size > NY) {
    if (mpi_rank == 0) {
      fprintf(stderr, "[Error] mpi_size=%d is larger than NY=%d\n", mpi_size, NY);
    }
    MPI_Finalize();
    return 1;
  }
  if (CFL_SUM2 > 1.0) {
    if (mpi_rank == 0) {
      fprintf(stderr, "[Error] unstable time step: CFLx^2 + CFLy^2 = %.6f > 1\n", CFL_SUM2);
    }
    MPI_Finalize();
    return 1;
  }

  init_simulation(mpi_rank, mpi_size);

  if (mpi_rank == 0) {
    printf("============================================\n");
    printf("  Redesigned 2D Wave Equation Example\n");
    printf("============================================\n");
    printf("Global grid          : %d x %d\n", NX, NY);
    printf("Physical size        : %.6f x %.6f\n", LX, LY);
    printf("Grid spacing         : dx=%.6e  dy=%.6e\n", DX, DY);
    printf("Time steps           : %d\n", NT);
    printf("MPI processes        : %d\n", mpi_size);
    printf("Thread groups/rank   : %d\n", N_GROUPS);
    printf("Threads/group target : %d\n", THREADS_PER_GROUP);
    printf("CFL numbers          : c*dt/dx=%.4f  c*dt/dy=%.4f\n", CFL_X, CFL_Y);
    printf("Resolution mode      : %s\n", USE_FIXED_DOMAIN ? "fixed domain" : "fixed spacing");
    printf("Local decomposition  : MPI(Y) + Group(Y) + Thread(X)\n");
    printf("============================================\n");
  }

  NProcPNode = mpi_size;
  err = InitThreads(
    mpi_rank,
    NCorePClu,
    NCluPNode,
    NCorePGrp,
    NThPGrp,
    NGrpPProc,
    NProcPNode,
    &ManageCoreId
  );

  if (err != 0) {
    fprintf(stderr, "[Error] MPI=%d: InitThreads failed with code %d\n", mpi_rank, err);
    free_simulation();
    MPI_Finalize();
    return 1;
  }

  setup_thread_tasks();

  StartThreads(thread_run);
  thread_run();
  EndThreads();

  free_simulation();
  MPI_Finalize();
  return 0;
}
