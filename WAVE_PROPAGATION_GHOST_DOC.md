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
} ThreadTask;
```

每个线程只关心：

- 它属于哪个组 `gid`
- 它在组内的线程号 `tid`
- 它负责的列区间 `[x_begin, x_end)`
- 它负责的行区间 `[y_begin, y_end)`

## 4. 每一步是怎么推进的

每个时间步都按下面顺序执行：

进入时间步循环前，主线程会先对初始 `u_curr` 做一次 halo 交换，用来建立第一步计算和初始能量诊断所需的边界数据。

进入循环后，每个时间步按下面顺序执行：

1. 主线程把 `u_next` 清零。
2. 主线程通知所有组开始当前步计算。
3. 每个组主线程再通知自己的工作线程开始。
4. 所有线程在共享的 `u_curr` 上读数据，在各自的 tile 上写 `u_next`。
5. 组主线程等待本组工作线程完成。
6. 主线程等待所有组完成。
7. 主线程交换 `u_prev / u_curr / u_next` 指针。
8. 主线程对新的 `u_curr` 做 MPI halo 交换。
9. 主线程计算离散总能量，并做 MPI 全局归约。

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

更新核使用二维二阶显式差分：

```text
u_next(y, x) =
  2*u_curr(y, x) - u_prev(y, x)
  + c^2 * dt^2 * (
      (u_curr(y, x-1) - 2*u_curr(y, x) + u_curr(y, x+1)) / dx^2
    + (u_curr(y-1, x) - 2*u_curr(y, x) + u_curr(y+1, x)) / dy^2
    )
```

这里：

- `Y` 向相邻点通过 MPI halo 行获得
- `X` 向相邻点直接从共享本地数组读取
- 左右物理边界固定为 0
- 最上/最下物理边界行也固定为 0，不参与时间推进
- 如果某个 MPI rank 没有上/下邻居，对应 halo 行也固定为 0

## 7. 离散总能量诊断

新版本输出的 `Energy` 已经不是旧版那种单纯的 `sum(u^2)`，而是一个和 leapfrog 时间推进更匹配的半步离散总能量诊断：

```text
E_h^(n+1/2) =
  1/2 Σ ((u^(n+1) - u^n) / dt)^2 dx dy
  + 1/2 c^2 Σ [D_x u^(n+1) · D_x u^n + D_y u^(n+1) · D_y u^n] dx dy
```

其中：

- 动能项使用相邻两个时间层的差分
- 势能项使用相邻两个时间层的空间差分内积
- `Y` 向边差只统计一次，因此不会在 MPI 进程边界重复计数

这比单纯统计 `u^2` 更接近当前离散格式真正应当保持稳定的能量量。
如果后续运行里仍然看到明显单调漂移，就更应该优先检查数值格式、边界条件或同步，而不是先怀疑“能量定义本身不对”。

## 8. 如何通过改分辨率增加计算量

代码顶部现在提供了两种模式：

```c
#define USE_FIXED_DOMAIN 1
```

### 8.1 方法一：固定物理区域，提高分辨率

推荐方式：

- 保持 `USE_FIXED_DOMAIN=1`
- 直接增大 `NX`、`NY`

这时代码会自动使用：

```c
DX = LX / (NX - 1)
DY = LY / (NY - 1)
```

含义是：

- 物理区域大小不变
- 网格更密
- 分辨率更高
- 计算量随网格点数增加

这适合做真正的“网格加密”实验。

### 8.2 方法二：固定网格步长，扩大物理区域

如果改成：

```c
#define USE_FIXED_DOMAIN 0
```

那么 `DX`、`DY` 保持常数，增大 `NX`、`NY` 的效果就会变成：

- 物理区域更大
- 网格点更多
- 计算量也会增加

这适合做“同样的空间步长下，求解更大区域”的实验。

### 8.3 稳定性条件

无论使用哪种模式，都需要满足显式格式的 CFL 条件：

```text
(c*dt/dx)^2 + (c*dt/dy)^2 <= 1
```

代码启动时会检查这一点，如果不满足会直接报错退出。

## 9. 这次重构去掉了什么

和旧版相比，这次设计主动移除了下面这些复杂度：

- 每线程独立波场
- 组内共享边界缓冲区
- 线程级 `X` 向 ghost 交换
- `thread_done / boundary_ready` 之类的旁路状态
- 同一 MPI 进程内组之间的专门 ghost 复制逻辑

新版本只保留一条主线：

> 进程间边界由 MPI 处理，进程内线程只做共享网格上的分块计算。

## 10. 建议阅读顺序

如果你想快速理解新示例，推荐按这个顺序看：

1. 先读主线程 `main_thread()`，看全局时间步顺序。
2. 再读 `group_main_thread()` 和 `worker_thread()`，看分层同步。
3. 再读 `compute_block()`，看真正的差分更新。
4. 最后看 `exchange_y_halos()`，理解 MPI halo 交换。

## 11. 关键代码位置

[`wave_propagation_ghost.c`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c) 里最值得先看的位置有：

- 进程级初始化与内存分配：[`wave_propagation_ghost.c:154`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c#L154)
- MPI halo 交换与物理边界处理：[`wave_propagation_ghost.c:215`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c#L215)
- 差分更新核：[`wave_propagation_ghost.c:253`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c#L253)
- 半步离散总能量诊断：[`wave_propagation_ghost.c:271`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c#L271)
- 线程任务切分：[`wave_propagation_ghost.c:316`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c#L316)
- 工作线程入口：[`wave_propagation_ghost.c:351`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c#L351)
- 组主线程入口：[`wave_propagation_ghost.c:363`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c#L363)
- 主线程入口：[`wave_propagation_ghost.c:377`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c#L377)
- 程序入口：[`wave_propagation_ghost.c:441`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c#L441)

## 12. 编译与运行

### 12.1 编译

```bash
mpicc -Wall -O2 -o wave_propagation_ghost \
    wave_propagation_ghost.c \
    mythread/mythread.c \
    -lpthread -lm
```

### 12.2 运行

```bash
mpirun -np 4 ./wave_propagation_ghost
```

## 13. 说明

- 这份示例优先追求“结构正确、便于讲解”，不是极致性能版本。
- 目标是给 `MPI + mythread` 提供一份更稳定的二维波动方程参考实现。
- 波动方程本身的数学背景仍然可以参考 [`WAVE_EXAMPLE.md`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/WAVE_EXAMPLE.md)。
