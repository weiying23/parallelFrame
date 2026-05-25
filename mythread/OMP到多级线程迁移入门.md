# OMP 到 mythread 多级线程迁移入门

本文面向第一次把 OpenMP 应用迁移到 `mythread` 的开发者。目标不是解释框架内部实现，而是回答三个实际问题：

1. 应用代码需要改哪些地方？
2. 原来的 `#pragma omp parallel/for/barrier/single/reduction` 分别怎么替换？
3. 怎样从一个最小两级线程版本，逐步扩展到多级线程版本？

可以对照阅读：

- `wave_propagation_ghost_mpi_omp.c`：MPI + OpenMP 版本。
- `wave_propagation_ghost_mpi_mythread.c`：基于上面文件改出的两级 `mythread_sync` 版本。
- `wave_propagation_ghost_fix.c`：更完整的多线程组、NUMA、组间 halo 版本。

---

## 1. 先建立心智模型

OpenMP 的模型是“进入一个并行区域，编译器帮你创建线程、分配循环、做 barrier”：

```c
#pragma omp parallel
{
  #pragma omp single
  do_serial_work();

  #pragma omp for
  for (...) compute(...);

  #pragma omp barrier
}
```

`mythread` 的模型是“你显式启动线程，每个线程进入同一个 `thread_run()`，再根据角色执行不同逻辑”：

```c
StartThreads(thread_run);
thread_run();      /* 当前主线程也进入 */
EndThreads();
```

线程角色分两种常用层级。

### 两级线程

适合从 OMP 快速迁移，结构简单：

```text
MMT 主线程
  ├── Worker 0
  ├── Worker 1
  └── Worker N
```

主线程负责串行工作、MPI 通信、发阶段信号、等待 worker 完成。worker 只负责自己那一块计算。

对应同步接口：

```text
主线程 → worker：mSetSubs(state)
主线程等 worker：mWaitSubs(state)
worker 等主线程：sWaitState(state)
worker → 主线程：sSetState(state)
```

### 多级线程

适合 NUMA、多线程组、每组独立数据分配：

```text
MMT 主线程
  ├── GMT 组主 0
  │   ├── Worker 0.1
  │   └── Worker 0.2
  └── GMT 组主 1
      ├── Worker 1.1
      └── Worker 1.2
```

对应同步接口：

```text
MMT → GMT：mSetGrps(state)
MMT 等 GMT：mWaitGrps(state)
GMT 等 MMT：gWaitMain(state)
GMT → MMT：gSetMain(state)

GMT → Worker：gSetSubs(state)
GMT 等 Worker：gWaitSubs(state)
Worker 等 GMT：sWaitGrp(state)
Worker → GMT：sSetGrp(state)
```

新手建议先实现两级线程，验证数值正确后，再升级成多级线程。

---

## 2. OMP 语法到 mythread 写法速查

| OpenMP 写法 | mythread 应用写法 |
|---|---|
| `#include <omp.h>` | `#include "mythread/mythread.h"` |
| `#pragma omp parallel` | `StartThreads(thread_run); thread_run(); EndThreads();` |
| `omp_get_thread_num()` | `ti->ind`，当前线程编号 |
| `omp_get_max_threads()` | `md.Nthreads - 1` 或配置里的 worker 数 |
| `#pragma omp single` | 放到主线程函数里执行 |
| `#pragma omp for schedule(static)` | 启动前手动给每个 worker 分块 |
| `#pragma omp barrier` | 用 `mSetSubs/mWaitSubs` 或 `GSS/GWS` 等状态同步 |
| `reduction(+:x)` | 每个线程写自己的 `partial_x`，主线程汇总 |
| `critical` | 尽量避免；优先让每个线程只写私有区域 |

---

## 3. 应用需要做的核心修改

### 修改 1：引入 mythread 头文件

从：

```c
#include <omp.h>
#include "mythread/mythread_config.h"
#include "mythread/mythread_decomp.h"
```

改成：

```c
#include "mythread/mythread.h"
```

`mythread.h` 已经包含线程、同步、配置、域分解等常用接口。

### 修改 2：定义每个线程自己的任务数据

OpenMP 里线程私有变量通常写在并行区域内部。`mythread` 里建议显式定义一个任务结构：

```c
typedef struct {
  int y_begin, y_end;
  int x_begin, x_end;
  double partial_energy;
} ThreadTask;
```

然后告诉框架每个线程要分配多少私有数据：

```c
int _gettdsize_(void) { return (int)sizeof(ThreadTask); }
int _getgdsize_(void) { return 0; }
```

在线程函数中这样取：

```c
ThreadTask *task = (ThreadTask*)ti->td;
```

### 修改 3：把并行循环改成 block 函数

OpenMP 版本通常直接并行整个循环：

```c
#pragma omp parallel for schedule(static)
for (int y = yb; y < ye; y++) {
  for (int x = xb; x < xe; x++) {
    compute_point(y, x);
  }
}
```

迁移时先把循环包装成“只算一个块”的函数：

```c
static void compute_block(Field *f,
                          int y_begin, int y_end,
                          int x_begin, int x_end) {
  for (int y = y_begin; y < y_end; y++) {
    for (int x = x_begin; x < x_end; x++) {
      compute_point(f, y, x);
    }
  }
}
```

这样 worker 只需要调用：

```c
compute_block(f, task->y_begin, task->y_end,
                 task->x_begin, task->x_end);
```

### 修改 4：启动前手动切分任务

OpenMP 的 `schedule(static)` 由编译器分配。`mythread` 需要应用自己切分，例如按 Y 方向均分：

```c
static void setup_thread_tasks(const mythread_tile *tile) {
  int nworkers = md.Nthreads - 1;
  int base = tile->ny / nworkers;
  int rem  = tile->ny % nworkers;
  int y = HALO;

  for (int i = 0; i < nworkers; i++) {
    ThreadTask *task = (ThreadTask*)md.threads[i].td;
    int rows = base + (i < rem ? 1 : 0);

    task->x_begin = HALO;
    task->x_end   = HALO + tile->nx;
    task->y_begin = y;
    task->y_end   = y + rows;
    task->partial_energy = 0.0;

    y += rows;
  }
}
```

注意：`setup_thread_tasks()` 要在 `InitThreads()` 之后调用，因为此时 `md.threads[i].td` 才已经分配好。

### 修改 5：用状态编号替代 barrier

`mythread_sync` 是基于状态值的同步。每个阶段使用一个单调递增的状态编号：

```c
static int init_state(void) { return 1; }
static int compute_state(int step) { return 2 * step + 2; }
static int boundary_state(int step) { return 2 * step + 3; }
static int energy_state(void) { return 2 * NT + 4; }
```

不要在不同阶段复用同一个状态值。状态值单调递增，能避免 worker 收到上一阶段残留信号。

### 修改 6：把 OpenMP 并行区域拆成两个函数

两级线程的最小结构是：

```c
static void worker_thread(void) {
  ThreadTask *task = (ThreadTask*)ti->td;

  sWaitState(init_state());
  init_block(task);
  sSetState(init_state());

  for (int step = 0; step < NT; step++) {
    int cs = compute_state(step);
    int bs = boundary_state(step);

    sWaitState(cs);
    compute_interior_block(task);
    sSetState(cs);

    sWaitState(bs);
    compute_boundary_block(task);
    sSetState(bs);
  }

  sWaitState(energy_state());
  task->partial_energy = compute_energy_block(task);
  sSetState(energy_state());
}
```

主线程负责串行部分和发信号：

```c
static void main_thread_run(void) {
  mSetSubs(init_state());
  mWaitSubs(init_state());

  exchange_halo();

  for (int step = 0; step < NT; step++) {
    apply_boundary();
    exchange_halo();

    mSetSubs(compute_state(step));
    mWaitSubs(compute_state(step));

    mSetSubs(boundary_state(step));
    mWaitSubs(boundary_state(step));

    swap_fields();
  }

  mSetSubs(energy_state());
  mWaitSubs(energy_state());
  reduce_worker_results();
}
```

线程入口只做角色分发：

```c
void thread_run(void) {
  if (ti->igrp < 0)
    main_thread_run();
  else
    worker_thread();
}
```

在 `NGrpPProc=1` 的两级模式下，`StartThreads()` 创建的是 sub threads，当前调用 `thread_run()` 的线程是 MMT，满足上面的判断。

### 修改 7：替换 main 里的线程生命周期

OpenMP 版本通常只需要设置线程数：

```c
omp_set_num_threads(nthreads);
run_simulation();
```

`mythread` 版本需要初始化、启动、结束：

```c
int manage_core = cfg_ManageCoreId;

CoreOffset  = cfg_CoreOffset;
ClustOffset = cfg_ClustOffset;

int err = InitThreads(mpi_rank,
                      cfg_NCorePClu,
                      cfg_NCluPNode,
                      cfg_NCorePGrp,
                      cfg_N_WORKERS + 1,  /* 总线程数 = worker + 主控槽 */
                      1,                  /* 两级模式：每进程 1 组 */
                      node_size,
                      &manage_core);
if (err != 0) {
  /* 报错并退出 */
}

setup_thread_tasks(tile);

StartThreads(thread_run);
thread_run();
EndThreads();
```

---

## 4. 从两级线程升级到多级线程

两级线程先跑通后，如果要使用多个线程组，需要做四类扩展。

### 扩展 1：线程角色变成三类

```c
void thread_run(void) {
  if (ti->igrp < 0) {
    main_thread();          /* MMT */
  } else if (ti->ind == 0) {
    group_main_thread();    /* GMT */
  } else {
    worker_thread();        /* Worker */
  }
}
```

### 扩展 2：同步从一跳变成两跳

两级线程：

```text
MMT -- mSetSubs/mWaitSubs -- Worker
```

多级线程：

```text
MMT -- mSetGrps/mWaitGrps -- GMT -- gSetSubs/gWaitSubs -- Worker
```

GMT 的典型写法：

```c
static void group_main_thread(void) {
  for (int step = 0; step < NT; step++) {
    int cs = compute_state(step);

    gWaitMain(cs);  /* 等 MMT */
    gSetSubs(cs);   /* 放行本组 worker */
    gWaitSubs(cs);  /* 等本组 worker */
    gSetMain(cs);   /* 通知 MMT */
  }
}
```

MMT 的典型写法：

```c
mSetGrps(cs);
mWaitGrps(cs);
```

Worker 的典型写法：

```c
sWaitGrp(cs);
compute_block(...);
sSetGrp(cs);
```

### 扩展 3：数据从一份 Field 变成每组一份 GroupField

两级线程通常所有 worker 共享一个 rank 内 `Field`，只切分 Y 区间。

多级线程为了 NUMA 和组间 halo，推荐每个 group 分配自己的 `GroupField`：

```c
GroupField *g_gfields;

/* GMT 中分配，保证 first touch / NUMA 亲和 */
group_alloc_field(&g_gfields[gid], gid, g_decomp);
```

所有组分配完成后，MMT 建立组间 halo buffer 连接：

```c
group_field_link_buffers(g_gfields, g_decomp);
```

### 扩展 4：halo 交换拆成组内和 MPI 两层

两级线程中一般只有 MPI halo：

```c
halo_x(&g_f);
```

多级线程中，rank 内不同 group 之间也要交换 halo：

```c
mythread_halo_exchange_intra(g_gfields, g_decomp);
mythread_halo_exchange_mpi(g_gfields, g_decomp, &g_mpi_ctx,
                           (uintptr_t)MPI_COMM_WORLD);
```

这就是为什么完整版本会比两级版本多一些 `GroupField`、`mythread_halo` 和 `mythread_decomp` 的代码。

---

## 5. 迁移时最容易出错的地方

### 状态编号不能复用

错误示例：

```c
for (int step = 0; step < NT; step++) {
  mSetSubs(1);
  mWaitSubs(1);
}
```

worker 可能把上一轮的状态当成下一轮信号。请使用随 step 递增的状态：

```c
int cs = compute_state(step);
mSetSubs(cs);
mWaitSubs(cs);
```

### `setup_thread_tasks()` 必须在 `InitThreads()` 后

`ti->td` 和 `md.threads[i].td` 是框架在线程系统初始化时分配的。太早访问会得到空指针。

正确顺序：

```text
InitThreads
setup_thread_tasks
StartThreads
thread_run
EndThreads
```

### 主线程也要调用 `thread_run()`

`StartThreads(thread_run)` 只启动框架创建的线程；当前调用者也要进入入口函数，成为 MMT：

```c
StartThreads(thread_run);
thread_run();
EndThreads();
```

### MPI 通信建议放在 MMT

除非你明确使用 `MPI_Init_thread` 并确认 MPI 支持多线程通信，否则最稳妥的方式是只让 MMT 调 MPI：

```c
MPI_Barrier(...);
MPI_Isend(...);
MPI_Irecv(...);
MPI_Allreduce(...);
```

worker 只做纯计算。

### 避免多个 worker 写同一个单元

从 OpenMP 迁移出来后，不要依赖隐式循环调度。确保每个 worker 的 `(x,y)` 区间互不重叠。归约值写到 `task->partial_*`，最后由主线程汇总。

---

## 6. 最小迁移清单

把一个 OMP stencil 应用改成两级 `mythread_sync`，按这个清单做：

- 替换头文件为 `mythread/mythread.h`。
- 定义 `ThreadTask`，实现 `_gettdsize_()`。
- 把主要计算循环改成 `*_block(..., y_begin, y_end, x_begin, x_end)`。
- 新增 `setup_thread_tasks()`，手动给 worker 切分区间。
- 为每个阶段定义单调递增状态编号。
- 新增 `worker_thread()`：等待状态、计算、上报状态。
- 新增 `main_thread_run()`：串行逻辑、MPI 通信、发信号、等待 worker。
- 新增 `thread_run()` 做 MMT/Worker 分发。
- 在 `main()` 中增加 `InitThreads -> setup_thread_tasks -> StartThreads -> thread_run -> EndThreads`。
- 用小算例先跑 `mpirun -np 1`，再跑 `np=2`，最后跑 `np=4` 覆盖 2D 分解。

---

## 7. 建议的验证步骤

1. 先用很小算例，例如 `NX=64, NY=64, NT=5`。
2. `np=1` 验证线程同步和单 rank 数值。
3. `np=2` 验证上下 halo。
4. `np=4` 验证左右 halo 和 2D MPI 分解。
5. 对比 OMP 版和 mythread 版的 final energy。
6. 再放大网格做性能测试。

如果小算例 energy 和 OMP 版一致，说明迁移后的任务划分、halo、同步顺序基本正确。

