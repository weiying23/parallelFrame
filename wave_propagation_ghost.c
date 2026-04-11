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

#define NX 32000
#define NY 8000
#define NT 480

#define DT 0.001
#define C0 0.1

/*
 * 增大计算量有两种常见方式：
 * 1. USE_FIXED_DOMAIN=1：固定物理区域，增大 NX/NY 表示网格加密、分辨率提高。
 * 2. USE_FIXED_DOMAIN=0：固定 DX/DY，增大 NX/NY 表示物理区域扩大、总网格点增加。
 */
#define USE_FIXED_DOMAIN 0

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
#define N_GROUPS 1
#define N_WORKERS 35
#define THREADS_PER_GROUP (N_WORKERS + 1)
#define ENERGY_REPORT_INTERVAL 60

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
  int cpu_id;
  int x_begin;
  int x_end;
  int y_begin;
  int y_end;
  double partial_energy;
  double t_wait_init;
  double t_work_init;
  double t_wait_energy0;
  double t_work_energy0;
  double t_wait_compute;
  double t_work_compute;
  double t_wait_boundary;
  double t_work_boundary;
  double t_wait_energy;
  double t_work_energy;
  double t_comm;
  double t_allreduce;
  int energy_steps;
} ThreadTask;

static SimulationData g_sim = {0};
static double *g_send_up = NULL;
static double *g_recv_up = NULL;
static double *g_send_down = NULL;
static double *g_recv_down = NULL;
static double *g_group_energy = NULL;

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

static inline double wall_time(void) {
  return MPI_Wtime();
}

static void reset_task_timers(ThreadTask *task) {
  if (!task) {
    return;
  }
  task->cpu_id = -1;
  task->partial_energy = 0.0;
  task->t_wait_init = 0.0;
  task->t_work_init = 0.0;
  task->t_wait_energy0 = 0.0;
  task->t_work_energy0 = 0.0;
  task->t_wait_compute = 0.0;
  task->t_work_compute = 0.0;
  task->t_wait_boundary = 0.0;
  task->t_work_boundary = 0.0;
  task->t_wait_energy = 0.0;
  task->t_work_energy = 0.0;
  task->t_comm = 0.0;
  task->t_allreduce = 0.0;
  task->energy_steps = 0;
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

static int initial_energy_state(void) {
  return 2;
}

static int compute_phase_state(int step) {
  return 3 * step + 3;
}

static int boundary_phase_state(int step) {
  return 3 * step + 4;
}

static int energy_phase_state(int step) {
  return 3 * step + 5;
}

static int should_measure_energy_step(int step) {
  if (step == 0 || step == NT - 1) {
    return 1;
  }
  if (ENERGY_REPORT_INTERVAL > 0 && ((step + 1) % ENERGY_REPORT_INTERVAL) == 0) {
    return 1;
  }
  return 0;
}

/* mythread treats NGrpPProc <= 1 as plain non-group mode rather than a one-group layout. */
static int uses_group_threads(void) {
  return ThreadG != 0;
}

static int init_fields_state(void) {
  return 1;
}

static void start_phase_from_main(int state) {
  if (uses_group_threads()) {
    mSetGrps(state);
  } else {
    mSetSubs(state);
  }
}

static void wait_phase_from_main(int state) {
  if (uses_group_threads()) {
    mWaitGrps(state);
  } else {
    mWaitSubs(state);
  }
}

static void wait_phase_from_worker(int state) {
  if (uses_group_threads()) {
    sWaitGrp(state);
  } else {
    sWaitState(state);
  }
}

static void finish_phase_from_worker(int state) {
  if (uses_group_threads()) {
    sSetGrp(state);
  } else {
    sSetState(state);
  }
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

static double initial_condition_value(int global_y,int x) {
  double cx = 0.5 * LX;
  double cy = 0.5 * LY;
  double sigma = 0.06 * ((LX < LY) ? LX : LY);
  double dx = x * DX - cx;
  double dy = global_y * DY - cy;

  return exp(-(dx * dx + dy * dy) / (2.0 * sigma * sigma));
}

static void initialize_field_block(const ThreadTask *task) {
  for (int y = task->y_begin; y < task->y_end; y++) {
    int global_y = global_y_from_local(y);

    for (int x = task->x_begin; x < task->x_end; x++) {
      int p = idx(y, x);
      double value = 0.0;

      if (global_y != 0 && global_y != NY - 1) {
        value = initial_condition_value(global_y, x);
      }
      g_sim.u_curr[p] = value;
      g_sim.u_prev[p] = value;
      g_sim.u_next[p] = 0.0;
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

static void free_group_energy(void) {
  free(g_group_energy);
  g_group_energy = NULL;
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

static void begin_y_halo_exchange_for(double *field,MPI_Request requests[4],int *request_count) {
  int count = 0;

  enforce_dirichlet_boundaries(field);
  zero_physical_y_boundaries(field);

  if (g_sim.neighbor_up >= 0) {
    memcpy(g_send_up, &field[idx(g_sim.local_ny, 0)], (size_t)NX * sizeof(double));
    MPI_Irecv(g_recv_up, NX, MPI_DOUBLE, g_sim.neighbor_up, 101, MPI_COMM_WORLD, &requests[count++]);
    MPI_Isend(g_send_up, NX, MPI_DOUBLE, g_sim.neighbor_up, 100, MPI_COMM_WORLD, &requests[count++]);
  } else {
    memset(&field[idx(g_sim.local_ny + HALO, 0)], 0, (size_t)NX * sizeof(double));
  }

  if (g_sim.neighbor_down >= 0) {
    memcpy(g_send_down, &field[idx(HALO, 0)], (size_t)NX * sizeof(double));
    MPI_Irecv(g_recv_down, NX, MPI_DOUBLE, g_sim.neighbor_down, 100, MPI_COMM_WORLD, &requests[count++]);
    MPI_Isend(g_send_down, NX, MPI_DOUBLE, g_sim.neighbor_down, 101, MPI_COMM_WORLD, &requests[count++]);
  } else {
    memset(&field[idx(0, 0)], 0, (size_t)NX * sizeof(double));
  }

  *request_count = count;
}

static void end_y_halo_exchange_for(double *field,MPI_Request requests[4],int request_count) {
  if (request_count > 0) {
    MPI_Waitall(request_count, requests, MPI_STATUSES_IGNORE);
  }

  if (g_sim.neighbor_up >= 0) {
    memcpy(&field[idx(g_sim.local_ny + HALO, 0)], g_recv_up, (size_t)NX * sizeof(double));
  } else {
    memset(&field[idx(g_sim.local_ny + HALO, 0)], 0, (size_t)NX * sizeof(double));
  }

  if (g_sim.neighbor_down >= 0) {
    memcpy(&field[idx(0, 0)], g_recv_down, (size_t)NX * sizeof(double));
  } else {
    memset(&field[idx(0, 0)], 0, (size_t)NX * sizeof(double));
  }

  enforce_dirichlet_boundaries(field);
  zero_physical_y_boundaries(field);
}

static void compute_block_rows(const ThreadTask *task,int y_begin,int y_end) {
  for (int y = y_begin; y < y_end; y++) {
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

static void compute_interior_block(const ThreadTask *task) {
  int y_begin = task->y_begin;
  int y_end = task->y_end;

  if (g_sim.neighbor_down >= 0 && y_begin < HALO + 1) {
    y_begin = HALO + 1;
  }
  if (g_sim.neighbor_up >= 0 && y_end > g_sim.local_ny) {
    y_end = g_sim.local_ny;
  }
  if (y_begin < y_end) {
    compute_block_rows(task, y_begin, y_end);
  }
}

static void compute_boundary_block(const ThreadTask *task) {
  int lower_row = HALO;
  int upper_row = g_sim.local_ny;

  if (g_sim.neighbor_down >= 0 && task->y_begin <= lower_row && lower_row < task->y_end) {
    compute_block_rows(task, lower_row, lower_row + 1);
  }
  if (g_sim.neighbor_up >= 0 &&
      upper_row != lower_row &&
      task->y_begin <= upper_row &&
      upper_row < task->y_end) {
    compute_block_rows(task, upper_row, upper_row + 1);
  }
}

static double compute_energy_block(const ThreadTask *task) {
  double kinetic = 0.0;
  double potential_x = 0.0;
  double potential_y = 0.0;
  double cell_area = DX * DY;
  int x_edge_begin;
  int x_edge_end;

  if (task->x_begin >= task->x_end || task->y_begin >= task->y_end) {
    return 0.0;
  }

  x_edge_begin = (task->x_begin == 1) ? 0 : task->x_begin;
  x_edge_end = task->x_end;

  for (int y = task->y_begin; y < task->y_end; y++) {
    int global_y = global_y_from_local(y);

    if (global_y > 0 && global_y < NY - 1) {
      for (int x = task->x_begin; x < task->x_end; x++) {
        double ut = (g_sim.u_curr[idx(y, x)] - g_sim.u_prev[idx(y, x)]) / DT;
        kinetic += ut * ut;
      }
    }

    for (int x = x_edge_begin; x < x_edge_end; x++) {
      double du_curr = (g_sim.u_curr[idx(y, x + 1)] - g_sim.u_curr[idx(y, x)]) / DX;
      double du_prev = (g_sim.u_prev[idx(y, x + 1)] - g_sim.u_prev[idx(y, x)]) / DX;
      potential_x += du_curr * du_prev;
    }

    if (global_y < NY - 1) {
      for (int x = task->x_begin; x < task->x_end; x++) {
        double du_curr = (g_sim.u_curr[idx(y + 1, x)] - g_sim.u_curr[idx(y, x)]) / DY;
        double du_prev = (g_sim.u_prev[idx(y + 1, x)] - g_sim.u_prev[idx(y, x)]) / DY;
        potential_y += du_curr * du_prev;
      }
    }
  }

  return 0.5 * (kinetic + C0 * C0 * (potential_x + potential_y)) * cell_area;
}

static double accumulate_worker_energy(void) {
  double local_energy = 0.0;

  if (uses_group_threads()) {
    for (int g = 0; g < md.ngrp; g++) {
      local_energy += g_group_energy[g];
    }
  } else {
    for (int t = 0; t < md.Nthreads - 1; t++) {
      ThreadTask *task = (ThreadTask*)md.threads[t].td;
      local_energy += task->partial_energy;
    }
  }

  return local_energy;
}

static double reduce_group_worker_energy(int group_id) {
  double energy = 0.0;
  threadGroup *pg = md.grps[group_id];

  for (int t = 1; t < pg->Nthreads; t++) {
    ThreadTask *task = (ThreadTask*)pg->threads[t].td;
    energy += task->partial_energy;
  }

  return energy;
}

static void swap_fields(void) {
  double *tmp = g_sim.u_prev;
  g_sim.u_prev = g_sim.u_curr;
  g_sim.u_curr = g_sim.u_next;
  g_sim.u_next = tmp;
}

static void setup_group_thread_tasks(void) {
  for (int g = 0; g < md.ngrp; g++) {
    threadGroup *pg = md.grps[g];
    int worker_count = pg->Nthreads - 1;
    int y_begin, y_end;

    split_range(HALO, g_sim.local_ny + HALO, md.ngrp, g, &y_begin, &y_end);

    for (int t = 0; t < pg->Nthreads; t++) {
      THREADINFO *pti = &pg->threads[t];
      ThreadTask *task = (ThreadTask*)pti->td;
      int x_begin, x_end;

      if (t == 0 || worker_count <= 0) {
        x_begin = 1;
        x_end = 1;
      } else {
        split_range(1, NX - 1, worker_count, t - 1, &x_begin, &x_end);
      }

      task->gid = g;
      task->tid = t;
      task->x_begin = x_begin;
      task->x_end = x_end;
      task->y_begin = y_begin;
      task->y_end = y_end;
      reset_task_timers(task);

#ifdef DEBUG
      printf(
        "[Init] MPI=%d Group=%d Thread=%d role=%s global-y=[%d,%d) x=[%d,%d)\n",
        mpi_id,
        g,
        t,
        (t == 0) ? "group-main" : "worker",
        g_sim.local_y_begin + (y_begin - HALO),
        g_sim.local_y_begin + (y_end - HALO),
        x_begin,
        x_end
      );
#endif // DEBUG
    }
  }
}

static void setup_single_group_tasks(void) {
  int worker_count = md.Nthreads - 1;
  int y_begin = HALO;
  int y_end = g_sim.local_ny + HALO;

  for (int t = 0; t < worker_count; t++) {
    THREADINFO *pti = &md.threads[t];
    ThreadTask *task = (ThreadTask*)pti->td;
    int x_begin, x_end;

    split_range(1, NX - 1, worker_count, t, &x_begin, &x_end);

    task->gid = 0;
    task->tid = t;
    task->x_begin = x_begin;
    task->x_end = x_end;
    task->y_begin = y_begin;
    task->y_end = y_end;
    reset_task_timers(task);

#ifdef DEBUG
    printf(
      "[Init] MPI=%d Group=%d Thread=%d role=worker global-y=[%d,%d) x=[%d,%d)\n",
      mpi_id,
      0,
      t,
      g_sim.local_y_begin,
      g_sim.local_y_end,
      x_begin,
      x_end
    );
#endif // DEBUG
  }
}

static void setup_thread_tasks(void) {
  if (uses_group_threads()) {
    setup_group_thread_tasks();
  } else {
    setup_single_group_tasks();
  }
}

static void worker_thread(void) {
  ThreadTask *task = (ThreadTask*)ti->td;
  double t0;

  t0 = wall_time();
  wait_phase_from_worker(init_fields_state());
  task->t_wait_init += wall_time() - t0;
  t0 = wall_time();
  initialize_field_block(task);
  task->t_work_init += wall_time() - t0;
  finish_phase_from_worker(init_fields_state());

  t0 = wall_time();
  wait_phase_from_worker(initial_energy_state());
  task->t_wait_energy0 += wall_time() - t0;
  t0 = wall_time();
  task->partial_energy = compute_energy_block(task);
  task->t_work_energy0 += wall_time() - t0;
  finish_phase_from_worker(initial_energy_state());

  for (int step = 0; step < NT; step++) {
    int compute_state = compute_phase_state(step);
    int boundary_state = boundary_phase_state(step);
    int energy_state = energy_phase_state(step);

    t0 = wall_time();
    wait_phase_from_worker(compute_state);
    task->t_wait_compute += wall_time() - t0;
    t0 = wall_time();
    compute_interior_block(task);
    task->t_work_compute += wall_time() - t0;
    finish_phase_from_worker(compute_state);

    t0 = wall_time();
    wait_phase_from_worker(boundary_state);
    task->t_wait_boundary += wall_time() - t0;
    t0 = wall_time();
    compute_boundary_block(task);
    task->t_work_boundary += wall_time() - t0;
    finish_phase_from_worker(boundary_state);

    if (should_measure_energy_step(step)) {
      t0 = wall_time();
      wait_phase_from_worker(energy_state);
      task->t_wait_energy += wall_time() - t0;
      t0 = wall_time();
      task->partial_energy = compute_energy_block(task);
      task->t_work_energy += wall_time() - t0;
      task->energy_steps += 1;
      finish_phase_from_worker(energy_state);
    }
  }
}


static void group_main_thread(void) {
  ThreadTask *task = (ThreadTask*)ti->td;
  double t0;

  t0 = wall_time();
  gWaitMain(init_fields_state());
  task->t_wait_init += wall_time() - t0;
  gSetSubs(init_fields_state());
  t0 = wall_time();
  gWaitSubs(init_fields_state());
  task->t_wait_init += wall_time() - t0;
  gSetMain(init_fields_state());

  t0 = wall_time();
  gWaitMain(initial_energy_state());
  task->t_wait_energy0 += wall_time() - t0;
  gSetSubs(initial_energy_state());
  t0 = wall_time();
  gWaitSubs(initial_energy_state());
  task->t_wait_energy0 += wall_time() - t0;
  t0 = wall_time();
  g_group_energy[ti->igrp] = reduce_group_worker_energy(ti->igrp);
  task->t_work_energy0 += wall_time() - t0;
  gSetMain(initial_energy_state());

  for (int step = 0; step < NT; step++) {
    int compute_state = compute_phase_state(step);
    int boundary_state = boundary_phase_state(step);
    int energy_state = energy_phase_state(step);

    t0 = wall_time();
    gWaitMain(compute_state);
    task->t_wait_compute += wall_time() - t0;
    gSetSubs(compute_state);
    t0 = wall_time();
    gWaitSubs(compute_state);
    task->t_wait_compute += wall_time() - t0;
    gSetMain(compute_state);

    t0 = wall_time();
    gWaitMain(boundary_state);
    task->t_wait_boundary += wall_time() - t0;
    gSetSubs(boundary_state);
    t0 = wall_time();
    gWaitSubs(boundary_state);
    task->t_wait_boundary += wall_time() - t0;
    gSetMain(boundary_state);

    if (should_measure_energy_step(step)) {
      t0 = wall_time();
      gWaitMain(energy_state);
      task->t_wait_energy += wall_time() - t0;
      gSetSubs(energy_state);
      t0 = wall_time();
      gWaitSubs(energy_state);
      task->t_wait_energy += wall_time() - t0;
      t0 = wall_time();
      g_group_energy[ti->igrp] = reduce_group_worker_energy(ti->igrp);
      task->t_work_energy += wall_time() - t0;
      task->energy_steps += 1;
      gSetMain(energy_state);
    }
  }
}

static void main_thread(void) {
  ThreadTask *task = (ThreadTask*)ti->td;
  double start_time;
  double local_energy = 0.0;
  double global_energy = 0.0;
  int current_halo_ready = 1;
  double t0;

  reset_task_timers(task);

  t0 = wall_time();
  start_phase_from_main(init_fields_state());
  wait_phase_from_main(init_fields_state());
  task->t_wait_init += wall_time() - t0;

  t0 = wall_time();
  exchange_y_halos_for(g_sim.u_curr);
  exchange_y_halos_for(g_sim.u_prev);
  task->t_comm += wall_time() - t0;

  t0 = wall_time();
  start_phase_from_main(initial_energy_state());
  wait_phase_from_main(initial_energy_state());
  task->t_wait_energy0 += wall_time() - t0;
  t0 = wall_time();
  local_energy = accumulate_worker_energy();
  task->t_work_energy0 += wall_time() - t0;
  t0 = wall_time();
  MPI_Allreduce(&local_energy, &global_energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  task->t_allreduce += wall_time() - t0;
  g_sim.initial_energy = global_energy;

  if (mpi_id == 0) {
    printf("[Main] Initial total energy: %.6f\n", g_sim.initial_energy);
  }

  start_time = MPI_Wtime();
  double prev_time = start_time;

  for (int step = 0; step < NT; step++) {
    int compute_state = compute_phase_state(step);
    int boundary_state = boundary_phase_state(step);
    int energy_state = energy_phase_state(step);
    int need_energy = should_measure_energy_step(step);
    MPI_Request requests[4];
    int request_count = 0;

    if (!current_halo_ready) {
      t0 = wall_time();
      begin_y_halo_exchange_for(g_sim.u_curr, requests, &request_count);
      task->t_comm += wall_time() - t0;
    }

    t0 = wall_time();
    start_phase_from_main(compute_state);
    wait_phase_from_main(compute_state);
    task->t_wait_compute += wall_time() - t0;

    if (!current_halo_ready) {
      t0 = wall_time();
      end_y_halo_exchange_for(g_sim.u_curr, requests, request_count);
      task->t_comm += wall_time() - t0;
      current_halo_ready = 1;
    }

    t0 = wall_time();
    start_phase_from_main(boundary_state);
    wait_phase_from_main(boundary_state);
    task->t_wait_boundary += wall_time() - t0;

    swap_fields();
    current_halo_ready = 0;

    if (need_energy) {
      t0 = wall_time();
      begin_y_halo_exchange_for(g_sim.u_curr, requests, &request_count);
      end_y_halo_exchange_for(g_sim.u_curr, requests, request_count);
      task->t_comm += wall_time() - t0;
      current_halo_ready = 1;

      t0 = wall_time();
      start_phase_from_main(energy_state);
      wait_phase_from_main(energy_state);
      task->t_wait_energy += wall_time() - t0;
      t0 = wall_time();
      local_energy = accumulate_worker_energy();
      task->t_work_energy += wall_time() - t0;
      t0 = wall_time();
      MPI_Allreduce(&local_energy, &global_energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      task->t_allreduce += wall_time() - t0;
      task->energy_steps += 1;

      if (mpi_id == 0) {
        double rel_diff = (g_sim.initial_energy > 0.0)
          ? fabs(global_energy - g_sim.initial_energy) / g_sim.initial_energy
          : 0.0;
        double cur_time = MPI_Wtime();
        double compute_time = cur_time - prev_time;
        prev_time = cur_time;

        printf(
          "[Main] Step %4d/%d, time %.3f,  Energy=%.6f  RelDiff=%.3e\n",
          step + 1,
          NT,
          compute_time,
          global_energy,
          rel_diff
        );
      }
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
  ThreadTask *task = (ThreadTask*)ti->td;
#ifdef DEBUG
  printf("[info] mpi %d gid %d, pid %d, c %d, getcore %d\n",
    mpi_id, ti->igrp, ti->ind, ti->indg, getcpuid());
#endif // DEBUG

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

static const char *task_role_name(const ThreadTask *task) {
  if (!task) {
    return "null";
  }
  if (task->gid < 0) {
    return "main";
  }
  if (uses_group_threads() && task->tid == 0) {
    return "group-main";
  }
  return "worker";
}

static double task_total_time(const ThreadTask *task) {
  if (!task) {
    return 0.0;
  }
  return task->t_wait_init + task->t_work_init +
         task->t_wait_energy0 + task->t_work_energy0 +
         task->t_wait_compute + task->t_work_compute +
         task->t_wait_boundary + task->t_work_boundary +
         task->t_wait_energy + task->t_work_energy +
         task->t_comm + task->t_allreduce;
}

static void print_one_task_timing(int mpi_rank,int mpi_size,int node_size,const ThreadTask *task) {
  int global_y_begin = -1;
  int global_y_end = -1;
  double total = task_total_time(task);

  if (task && task->y_end > task->y_begin && task->y_begin >= HALO) {
    global_y_begin = g_sim.local_y_begin + (task->y_begin - HALO);
    global_y_end = g_sim.local_y_begin + (task->y_end - HALO);
  }

  // printf(
    // "[Timing] rank %d/%d node_size=%d gid=%d tid=%d role=%s cpu=%d "
    // "x=[%d,%d) y=[%d,%d) global-y=[%d,%d) "
    // "wait(init=%.6f energy0=%.6f compute=%.6f boundary=%.6f energy=%.6f) "
    // "work(init=%.6f energy0=%.6f compute=%.6f boundary=%.6f energy=%.6f) "
    // "comm=%.6f allreduce=%.6f energy_steps=%d total=%.6f\n",
    // mpi_rank,
    // mpi_size,
    // node_size,
    // task ? task->gid : -999,
    // task ? task->tid : -999,
    // task_role_name(task),
    // task ? task->cpu_id : -1,
    // task ? task->x_begin : 0,
    // task ? task->x_end : 0,
    // task ? task->y_begin : 0,
    // task ? task->y_end : 0,
    // global_y_begin,
    // global_y_end,
    // task ? task->t_wait_init : 0.0,
    // task ? task->t_wait_energy0 : 0.0,
    // task ? task->t_wait_compute : 0.0,
    // task ? task->t_wait_boundary : 0.0,
    // task ? task->t_wait_energy : 0.0,
    // task ? task->t_work_init : 0.0,
    // task ? task->t_work_energy0 : 0.0,
    // task ? task->t_work_compute : 0.0,
    // task ? task->t_work_boundary : 0.0,
    // task ? task->t_work_energy : 0.0,
    // task ? task->t_comm : 0.0,
    // task ? task->t_allreduce : 0.0,
    // task ? task->energy_steps : 0,
    // total
  // );
}

static void print_timing_report(int mpi_rank,int mpi_size,int node_size) {
  for (int r = 0; r < mpi_size; r++) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (mpi_rank == r) {
      printf("========== Timing Report rank %d/%d ==========\n", mpi_rank, mpi_size);
      print_one_task_timing(mpi_rank, mpi_size, node_size, (ThreadTask*)md.tm.td);
      if (uses_group_threads()) {
        for (int g = 0; g < md.ngrp; g++) {
          threadGroup *pg = md.grps[g];
          for (int t = 0; t < pg->Nthreads; t++) {
            print_one_task_timing(mpi_rank, mpi_size, node_size, (ThreadTask*)pg->threads[t].td);
          }
        }
      } else {
        int worker_count = md.Nthreads - 1;
        for (int t = 0; t < worker_count; t++) {
          print_one_task_timing(mpi_rank, mpi_size, node_size, (ThreadTask*)md.threads[t].td);
        }
      }
      fflush(stdout);
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

int main(int argc,char **argv) {
  int mpi_rank, mpi_size;
  int NCorePClu = 38;
  int NCluPNode = 16;
  int NCorePGrp = 37;
  int NThPGrp = THREADS_PER_GROUP;
  int NGrpPProc = N_GROUPS;
  int NProcPNode = 16;
  int ManageCoreId = 36;
  int err;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  
  MPI_Comm node_comm;
  MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &node_comm);
  int node_size;
  MPI_Comm_size(node_comm, &node_size);

  NProcPNode = node_size;

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
    printf("Energy diagnostics   : initial + every %d steps + final\n", ENERGY_REPORT_INTERVAL);
    printf("============================================\n");
  }
  fflush(stdout);

  // for (int r = 0; r < mpi_size; r++) {
    MPI_Barrier(MPI_COMM_WORLD);
    // if (mpi_rank == r) {
      printf(
        "[Domain] rank %d/%d node_size=%d global-y=[%d,%d) local_ny=%d halo=%d neighbors(down=%d up=%d) x-update=[1,%d)\n",
        mpi_rank,
        mpi_size,
        node_size,
        g_sim.local_y_begin,
        g_sim.local_y_end,
        g_sim.local_ny,
        HALO,
        g_sim.neighbor_down,
        g_sim.neighbor_up,
        NX - 1
      );
      // fflush(stdout);
    // }
  // }
  MPI_Barrier(MPI_COMM_WORLD);

  // NProcPNode = mpi_size;
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
  if (uses_group_threads()) {
    g_group_energy = (double*)xcalloc((size_t)md.ngrp, sizeof(double));
  }

  StartThreads(thread_run);
  thread_run();
  EndThreads();

  print_timing_report(mpi_rank, mpi_size, node_size);

  free_group_energy();
  free_simulation();
  MPI_Finalize();
  return 0;
}
