# 重新设计的二维波动方程示例说明

## 1. 总览

`wave_propagation_ghost.c` 现在是一份更直接、更容易验证正确性的 `MPI + mythread` 混合并行示例，用来计算二维波动方程：

```text
u_tt = c^2 (u_xx + u_yy)
```

这次重构的核心目标不是继续扩展旧版“线程也带 ghost cell”的方案，而是先把下面这件事做扎实：

> 每个 MPI 进程只维护一份共享局部网格，MPI 负责进程间 halo 交换，线程只负责计算不重叠的 tile。

这样做之后，时间推进顺序、数据所有权和同步关系都会清晰很多。

## 2. 新的并行设计

### 2.1 分解方式

新版本把并行拆成三层：

1. `MPI` 沿 `Y` 方向切分全局网格。
2. 每个 MPI 进程内部，`mythread` 的线程组继续切分本地 `Y` 行。
3. 每个线程组内部，再由线程切分 `X` 方向内部列。

可以把它看成：

```text
全局网格
  -> MPI 进程拥有一个 Y 向子区域（含上下各一层 halo）
     -> 线程组拥有本进程的一段局部 Y 行
        -> 线程拥有该组内部的一块矩形 tile
```

### 2.2 和旧版 ghost 示例的根本区别

旧版设计里：

- 每个线程有自己的波场数组
- 进程间要交换边界
- 同一进程内组之间还要交换边界
- 组内线程之间还要交换边界

这会让“谁拥有哪份数据、什么时候交换、交换后谁再继续算”变得很难理清。

新版本改成：

- 每个 MPI 进程只有一份共享的 `u_prev / u_curr / u_next`
- 只保留进程间 `Y` 向 halo 交换
- 所有线程共享读取 `u_curr`
- 每个线程只写自己负责的 `u_next` tile

因此：

- 读依赖统一来自旧时间层
- 写入目标统一是新时间层
- 进程内不再需要线程级 ghost 交换
- 主线程可以统一做 MPI 通信和时间层轮换

## 3. 主要数据结构

### 3.1 `SimulationData`

进程级共享状态：

```c
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
```

含义是：

- `u_prev / u_curr / u_next`：当前 MPI 进程共享的三层时间数据
- `local_y_begin / local_y_end`：该进程负责的全局 `Y` 范围
- `local_ny`：本地内部行数
- `neighbor_up / neighbor_down`：上下相邻 MPI rank

本地数组实际分配的是 `local_ny + 2` 行，多出来的两行就是上下 halo。

### 3.2 `ThreadTask`

线程级计算块描述：

```c
typedef struct {
  int gid;
  int tid;
  int x_begin;
  int x_end;
  int y_begin;
  int y_end;
  double energy;
} ThreadTask;
```

每个线程只关心：

- 它属于哪个组 `gid`
- 它在组内的线程号 `tid`
- 它负责的列区间 `[x_begin, x_end)`
- 它负责的行区间 `[y_begin, y_end)`
- 它这一步计算出来的局部能量 `energy`

## 4. 每一步是怎么推进的

每个时间步都按下面顺序执行：

1. 主线程对 `u_curr` 做 MPI halo 交换。
2. 主线程把 `u_next` 清零。
3. 主线程通知所有组开始当前步计算。
4. 每个组主线程再通知自己的工作线程开始。
5. 所有线程在共享的 `u_curr` 上读数据，在各自的 tile 上写 `u_next`。
6. 组主线程等待本组工作线程完成。
7. 主线程等待所有组完成。
8. 主线程汇总本进程能量，并做 MPI 全局归约。
9. 主线程交换 `u_prev / u_curr / u_next` 指针。

这套顺序有个很重要的性质：

- 本步计算只读旧层 `u_curr`
- 本步结果只写新层 `u_next`

不会再出现“有些线程已经换到了下一层，有些线程还在读上一层”的混乱情况。

## 5. 同步关系

虽然数据布局简化了，但线程协作仍然沿用 `mythread` 的三层结构。

### 5.1 主线程 <-> 组主线程

- `mSetGrps(state)`：主线程发起当前步
- `mWaitGrps(state)`：主线程等待所有组完成
- `gWaitMain(state)`：组主线程等待主线程信号
- `gSetMain(state)`：组主线程通知主线程“本组完成”

### 5.2 组主线程 <-> 工作线程

- `gSetSubs(state)`：组主线程唤醒本组工作线程
- `gWaitSubs(state)`：组主线程等待工作线程完成
- `sWaitGrp(state)`：工作线程等待组主线程信号
- `sSetGrp(state)`：工作线程汇报“我这块算完了”

### 5.3 为什么状态值要递增

新版本使用：

```c
int state = step + 1;
```

也就是每一步都使用新的同步状态值，而不是反复复用同一个常量。这样可以避免旧版示例里那种“第一轮直接穿过、下一轮又可能卡死”的状态复用问题。

## 6. 数值更新公式

更新核使用标准的二维二阶显式差分：

```text
u_next(y, x) =
  2*u_curr(y, x) - u_prev(y, x)
  + (c*dt/dx)^2 * (
      u_curr(y-1, x) + u_curr(y+1, x)
    + u_curr(y, x-1) + u_curr(y, x+1)
    - 4*u_curr(y, x)
    )
```

这里：

- `Y` 向相邻点通过 MPI halo 行获得
- `X` 向相邻点直接从共享本地数组读取
- 左右物理边界固定为 0
- 如果某个 MPI rank 没有上/下邻居，对应 halo 行也固定为 0

## 7. 这次重构去掉了什么

和旧版相比，这次设计主动移除了下面这些复杂度：

- 每线程独立波场
- 组内共享边界缓冲区
- 线程级 `X` 向 ghost 交换
- `thread_done / boundary_ready` 之类的旁路状态
- 同一 MPI 进程内组之间的专门 ghost 复制逻辑

新版本只保留一条主线：

> 进程间边界由 MPI 处理，进程内线程只做共享网格上的分块计算。

## 8. 建议阅读顺序

如果你想快速理解新示例，推荐按这个顺序看：

1. 先读主线程 `main_thread()`，看全局时间步顺序。
2. 再读 `group_main_thread()` 和 `worker_thread()`，看分层同步。
3. 再读 `compute_block()`，看真正的差分更新。
4. 最后看 `exchange_y_halos()`，理解 MPI halo 交换。

## 9. 关键代码位置

[`wave_propagation_ghost.c`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c) 里最值得先看的位置有：

- 进程级初始化与内存分配：[`wave_propagation_ghost.c:127`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c#L127)
- MPI halo 交换：[`wave_propagation_ghost.c:179`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c#L179)
- 差分更新核：[`wave_propagation_ghost.c:211`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c#L211)
- 线程任务切分：[`wave_propagation_ghost.c:251`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c#L251)
- 工作线程入口：[`wave_propagation_ghost.c:287`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c#L287)
- 组主线程入口：[`wave_propagation_ghost.c:299`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c#L299)
- 主线程入口：[`wave_propagation_ghost.c:313`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c#L313)
- 程序入口：[`wave_propagation_ghost.c:383`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c#L383)

## 10. 编译与运行

### 10.1 编译

```bash
mpicc -Wall -O2 -o wave_propagation_ghost \
    wave_propagation_ghost.c \
    mythread/mythread.c \
    -lpthread -lm
```

### 10.2 运行

```bash
mpirun -np 4 ./wave_propagation_ghost
```

## 11. 说明

- 这份示例优先追求“结构正确、便于讲解”，不是极致性能版本。
- 目标是给 `MPI + mythread` 提供一份更稳定的二维波动方程参考实现。
- 波动方程本身的数学背景仍然可以参考 [`WAVE_EXAMPLE.md`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/WAVE_EXAMPLE.md)。
