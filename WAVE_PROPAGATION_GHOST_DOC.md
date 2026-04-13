# 重新设计的二维波动方程示例说明

## 1. 总览

`wave_propagation_ghost.c` 是一个 `MPI + mythread` 的二维波动方程示例，用来计算：

```text
u_tt = c^2 (u_xx + u_yy)
```

当前版本的目标有三点：

1. 保持数值逻辑清晰，边界与能量诊断自洽。
2. 尽量把计算下放到从线程，减轻 `main_thread()` 和 `group_main_thread()` 的负担。
3. 减少多组线程下的远程内存访问，把大块数组尽量 first-touch 到对应线程组附近的本地内存。

## 2. 并行结构

整体分解是三层：

1. `MPI` 采用二维进程网格 `(Px, Py)`，同时切分全局 `X` 与 `Y` 区间（每个 rank 持有一个矩形子域）。
2. 每个 MPI rank 内部，线程组再次把本地子域切分为二维网格 `(Gx, Gy)`（每个组持有一个矩形 tile）。
3. 每个组内部，worker 线程沿 `X` 方向在该组 tile 内均分列区间；`Y` 区间在组内共享。

可以把它看成：

```text
全局网格
  -> MPI rank 拥有一个 (X,Y) 矩形子域
     -> 每个线程组拥有该 rank 的一个 (X,Y) 矩形 tile
        -> 每个 worker 线程拥有该 tile 中的一段 X（共享同一段 Y）
```

需要注意的是，`mythread` 中 `NGrpPProc <= 1` 并不表示“只有 1 个组”，而是直接退回非分组模式。因此当前示例同时兼容两种执行方式：

- `ThreadG != 0`：分组模式，使用 `mSetGrps / gWaitMain / sWaitGrp` 这套接口。
- `ThreadG == 0`：非分组模式，退回 `mSetSubs / sWaitState` 这套接口。

兼容入口在 [`wave_propagation_ghost.c:221`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c#L221)。

## 3. 主要数据结构

### 3.1 `SimulationData`

进程级共享状态：

```c
typedef struct {
  double *u_prev;
  double *u_curr;
  double *u_next;
  int local_x_begin;
  int local_x_end;
  int local_nx;
  int local_y_begin;
  int local_y_end;
  int local_ny;
  int mpi_rank;
  int mpi_size;
  int proc_x;
  int proc_y;
  int proc_px;
  int proc_py;
  int neighbor_left;
  int neighbor_right;
  int neighbor_up;
  int neighbor_down;
  double initial_energy;
} SimulationData;
```

含义：

- `u_prev / u_curr / u_next`：当前 MPI 进程共享的三个时间层。
- `local_x_begin / local_x_end`：本 rank 负责的全局 `X` 区间。
- `local_y_begin / local_y_end`：本 rank 负责的全局 `Y` 区间。
- `local_nx / local_ny`：本地真实物理列/行数。
- `proc_(x,y) / proc_(px,py)`：MPI 二维进程网格与本 rank 的网格坐标。
- `neighbor_left/right/up/down`：四个方向相邻 MPI rank（用于 halo 交换）。

局部数组四周都有 halo，因此实际分配尺寸是 `(local_ny + 2*HALO) x (local_nx + 2*HALO)`。

### 3.2 `ThreadTask`

线程级任务描述：

```c
typedef struct {
  int gid;
  int tid;
  int x_begin;
  int x_end;
  int y_begin;
  int y_end;
  double partial_energy;
} ThreadTask;
```

含义：

- `gid / tid`：所属组和组内线程号。
- `x_begin / x_end`：该线程负责的 `X` 区间。
- `y_begin / y_end`：该线程负责的本地 `Y` 区间。
- `partial_energy`：该线程负责 tile 的局部离散总能量。

## 4. 线程职责

### 4.1 主线程 `main_thread()`

主线程现在只负责：

1. 发起阶段同步。
2. 做 halo 通信。
3. 执行 `swap_fields()`。
4. 在需要时做 MPI 全局能量归约。

它不再负责大块差分更新，也不再扫描整个局部网格计算能量。

### 4.2 组主线程 `group_main_thread()`

组主线程现在只负责：

1. 等待主线程发来的阶段状态。
2. 唤醒本组从线程。
3. 等待本组从线程完成。
4. 在能量阶段做组内小归约。

它不再直接调用差分更新核。

### 4.3 从线程 `worker_thread()`

从线程负责全部主要数值工作：

1. 并行初始化本线程负责 tile 上的 `u_prev / u_curr / u_next`。
2. 计算内部行更新。
3. 计算边界行更新。
4. 在需要输出诊断时，计算本线程 tile 上的离散总能量。

## 5. 当前版本做过的关键优化

### 5.1 能量诊断按需触发

能量诊断不再每步都做，而是只在：

- 初始时刻
- 每 `ENERGY_REPORT_INTERVAL` 步
- 最后一步

默认参数在代码里是：

```c
#define ENERGY_REPORT_INTERVAL 60
```

这样可以显著减少每步都触发的能量阶段和 `MPI_Allreduce` 次数。

### 5.2 halo 通信与内部计算重叠

只有紧贴 MPI 边界的首末物理行依赖远端 halo，因此每步更新拆成两部分：

1. 内部行阶段：不依赖新的 MPI halo。
2. 边界行阶段：依赖新的 MPI halo。

主线程在 halo 通信进行时，先让从线程计算内部行；通信完成后，再补算边界行。这样可以把一部分 MPI 等待隐藏到从线程计算里。

非阻塞 halo 通信入口在 [`wave_propagation_ghost.c:325`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c#L325)。

### 5.3 通过并行 first-touch 缓解远程内存访问

这是当前针对多组线程性能问题新增的重点优化。

之前 `u_prev / u_curr` 的初始写入主要由主线程完成，容易导致：

- 页面 first-touch 落在主线程所在 cluster / NUMA 节点
- 其他线程组后续计算时大量读取远程内存
- 多组模式下不同 worker 的 `compute_interior_block()` 时间差显著放大

现在改成：

1. 先完成线程任务划分。
2. 主线程发起专门的“初始化阶段”。
3. 每个从线程在自己的 tile 上并行写入 `u_prev / u_curr / u_next`。

这样 bulk 页面会优先由负责该 tile 的 worker 触页，更容易落到该组本地内存。初始化核在 [`wave_propagation_ghost.c:327`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c#L327)。

## 6. 时间推进顺序

### 6.1 初始化阶段

进入时间循环前，顺序如下：

1. 主线程发起 `init_fields_state()`。
2. 所有 worker 并行初始化各自 tile 上的 `u_prev / u_curr / u_next`。
3. 主线程对 `u_curr` 和 `u_prev` 做初始 halo 交换。
4. 主线程发起初始能量阶段。
5. worker 计算局部能量，组主线程做组内归约，主线程做 MPI 全局归约。

### 6.2 每一步的执行顺序

对第 `step` 步：

1. 如有需要，主线程先发起新的 `u_curr` halo 通信。
2. 主线程发起内部行阶段。
3. worker 计算不依赖远端 halo 的内部行。
4. 主线程等待 halo 通信完成。
5. 主线程发起边界行阶段。
6. worker 计算依赖 halo 的边界行。
7. 主线程执行 `swap_fields()`。
8. 如果当前步需要输出能量：
   主线程先补齐新的 `u_curr` halo，再发起能量阶段。
9. worker 计算局部能量，组主线程做组内归约，主线程做 MPI 全局归约。

## 7. 同步状态

当前使用的状态值是：

```c
init_fields_state() = 1
initial_energy_state() = 2
compute_phase_state(step) = 3 * step + 3
boundary_phase_state(step) = 3 * step + 4
energy_phase_state(step) = 3 * step + 5
```

含义：

- `1`：并行初始化阶段
- `2`：初始能量阶段
- `3*step+3`：内部行阶段
- `3*step+4`：边界行阶段
- `3*step+5`：能量阶段

这些入口分别在：

- [`wave_propagation_ghost.c:225`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c#L225)
- [`wave_propagation_ghost.c:194`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c#L194)
- [`wave_propagation_ghost.c:198`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c#L198)
- [`wave_propagation_ghost.c:202`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c#L202)
- [`wave_propagation_ghost.c:206`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c#L206)

## 8. 数值更新与能量

更新核仍然使用二维显式二阶差分：

```text
u_next(y, x) =
  2*u_curr(y, x) - u_prev(y, x)
  + c^2 * dt^2 * (
      (u_curr(y, x-1) - 2*u_curr(y, x) + u_curr(y, x+1)) / dx^2
    + (u_curr(y-1, x) - 2*u_curr(y, x) + u_curr(y+1, x)) / dy^2
    )
```

边界规则：

- 左右边界固定为 0。
- 全局最上、最下物理边界固定为 0。
- MPI 之间交换 `X` 与 `Y` 两个方向的 halo。

离散总能量按 tile 分配给 worker，再做组内和全局归约。能量核在 [`wave_propagation_ghost.c:623`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c#L623)。

## 9. 分辨率与稳定性

代码提供两种分辨率模式：

```c
#define USE_FIXED_DOMAIN 0
```

- `USE_FIXED_DOMAIN=1`：固定物理区域，增大 `NX / NY` 表示网格加密。
- `USE_FIXED_DOMAIN=0`：固定 `DX / DY`，增大 `NX / NY` 表示物理区域扩大。

显式格式需要满足 CFL 条件：

```text
(c*dt/dx)^2 + (c*dt/dy)^2 <= 1
```

程序启动时会检查：

- `CFL_X = C0 * DT / DX`
- `CFL_Y = C0 * DT / DY`
- `CFL_SUM2 = CFL_X^2 + CFL_Y^2`

## 10. 关键代码位置

关键位置都在 [`wave_propagation_ghost.c`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c)：

- 并行初始化 first-touch：[`wave_propagation_ghost.c:327`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c#L327)
- 进程域初始化与内存分配：[`wave_propagation_ghost.c:296`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c#L296)
- 非阻塞 halo 交换：[`wave_propagation_ghost.c:446`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c#L446)
- 内部行更新：[`wave_propagation_ghost.c:554`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c#L554)
- 边界行更新：[`wave_propagation_ghost.c:579`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c#L579)
- 局部能量计算：[`wave_propagation_ghost.c:623`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c#L623)
- 任务划分：[`wave_propagation_ghost.c:705`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c#L705)
- worker 入口：[`wave_propagation_ghost.c:807`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c#L807)
- 组主线程入口：[`wave_propagation_ghost.c:862`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c#L862)
- 主线程入口：[`wave_propagation_ghost.c:927`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c#L927)

## 11. 编译与运行

### 11.1 编译

```bash
mpicc -Wall -O2 -o wave_propagation_ghost \
    wave_propagation_ghost.c \
    mythread/mythread.c \
    -lpthread -lm
```

如需做小规模快速验证或调参，可以在编译时用 `-D` 覆盖默认宏，例如：

```bash
mpicc -Wall -O2 -DNX=256 -DNY=128 -DNT=10 -DN_GROUPS=4 -DN_WORKERS=8 \
    -o wave_propagation_ghost_small wave_propagation_ghost.c mythread/mythread.c -lpthread -lm
```

### 11.2 运行

```bash
mpirun -np 4 ./wave_propagation_ghost
```

## 12. 本次针对远程内存访问的修改摘要

这次新增的直接修改点有：

- 去掉主线程里的 bulk 初值写入，改成 worker 并行初始化。
- 新增 `init_fields_state()` 阶段，用来统一调度并行 first-touch。
- 把 `u_prev / u_curr / u_next` 的 bulk 页面尽量 first-touch 到对应 worker 所在组附近。
- 线程组的任务划分从“只沿 `Y` 切分”更新为“组层二维切分本地 `(X,Y)` 子域”；组内 worker 在该组 tile 内沿 `X` 继续均分。
- 保留 `N_GROUPS=1` 时的兼容逻辑，单组退回非分组模式时仍然能走同样的并行初始化。
- 把 worker 内部的调试计时改成真正只包 `compute_interior_block()`，并在 `DEBUG` 下才输出，避免打印本身继续污染性能测量。

