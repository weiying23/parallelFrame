# 重新设计的二维波动方程示例说明

## 1. 总览

`wave_propagation_ghost.c` 是一个 `MPI + mythread` 的二维波动方程示例，用来计算：

```text
u_tt = c^2 (u_xx + u_yy)
```

这版实现的目标已经从“先把逻辑做对”进一步推进到“尽量削弱主线程瓶颈”。现在的核心原则是：

1. 每个 MPI 进程只维护一份共享局部网格 `u_prev / u_curr / u_next`。
2. 主线程 `main_thread()` 只做 MPI 通信、阶段调度、指针交换和少量归约。
3. 组主线程 `group_main_thread()` 不再做数值更新，只负责组内同步和小规模归约。
4. 真正的大头计算尽量下放到从线程 `worker_thread()`。

## 2. 并行分解

并行结构仍然是三层：

1. `MPI` 沿 `Y` 方向切分全局网格。
2. 每个 MPI 进程内部，线程组继续沿本地 `Y` 方向切分。
3. 每个组内部，只让从线程沿 `X` 方向切分实际计算区间。

可以把它理解成：

```text
全局网格
  -> MPI rank 拥有一个 Y 子区间
     -> 每个线程组拥有该 rank 中的一段本地 Y
        -> 每个从线程拥有该组中的一段 X
```

和之前版本相比，组主线程不再拥有自己的计算 tile，它现在是纯协调线程。

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

含义：

- `u_prev / u_curr / u_next`：当前 MPI 进程共享的三个时间层。
- `local_y_begin / local_y_end`：本 rank 负责的全局 `Y` 区间。
- `local_ny`：本地真实物理行数。
- `neighbor_up / neighbor_down`：上下相邻 MPI rank。

局部数组带上下两行 halo，所以实际分配行数是 `local_ny + 2`。

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

- `gid / tid`：所属线程组以及组内线程号。
- `x_begin / x_end`：该从线程负责的 `X` 区间。
- `y_begin / y_end`：该线程所在组负责的本地 `Y` 区间。
- `partial_energy`：该从线程负责 tile 的局部离散总能量。

## 4. 线程职责

### 4.1 主线程 `main_thread()`

主线程现在只负责：

1. 初始 halo 交换。
2. 发起各个同步阶段。
3. 在需要时启动 MPI 非阻塞 halo 交换并等待完成。
4. 在每步结束后执行 `swap_fields()`。
5. 在需要输出诊断时做 MPI 全局能量归约。

它不再做整块能量扫描，也不再直接参与差分更新。

### 4.2 组主线程 `group_main_thread()`

组主线程现在只负责：

1. 等待主线程发来的阶段状态。
2. 唤醒本组从线程。
3. 等待本组从线程结束。
4. 在能量阶段把本组从线程的 `partial_energy` 做一次组内小归约。
5. 向主线程回报本组完成。

它不再调用差分更新核。

### 4.3 从线程 `worker_thread()`

从线程负责全部大头计算：

1. 计算内部行的差分更新。
2. 计算边界行的差分更新。
3. 在需要输出诊断时，计算本线程负责 tile 的离散总能量。

## 5. 现在的优化重点

这版相对上一版新增了两类重要优化。

### 5.1 能量诊断不再每步计算

之前每一步都会：

1. 启动能量阶段。
2. 让所有从线程再扫一遍自己的数据。
3. 主线程做 `MPI_Allreduce`。

这会让主线程每一步都多出一次全局同步。

现在改成只在“需要输出日志的步数”计算能量：

- 第 1 步
- 每 `ENERGY_REPORT_INTERVAL` 步
- 最后一步

代码默认是：

```c
#define ENERGY_REPORT_INTERVAL 100
```

这样 `NT=400` 时，会输出：

- `Step 1`
- `Step 100`
- `Step 200`
- `Step 300`
- `Step 400`

其他步直接跳过能量阶段和 `MPI_Allreduce`。

### 5.2 halo 交换与内部计算重叠

显式二维波动方程里，只有紧贴 MPI 分界面的本地首行和末行依赖远端 halo。

因此更新被拆成了两部分：

1. 内部行更新：不依赖新的 MPI halo。
2. 边界行更新：依赖新的 MPI halo。

主线程现在的做法是：

1. 如果当前 `u_curr` 的 halo 还没准备好，就先发起一次非阻塞 halo 交换。
2. 在通信进行时，让从线程先算内部行。
3. 内部行结束后，主线程等待通信完成。
4. 再发起边界行阶段，让从线程补算边界行。

这意味着主线程上的通信等待，会尽量被从线程的内部行计算掩盖掉。

## 6. 每个时间步如何推进

### 6.1 初始阶段

在进入时间循环前：

1. 主线程对 `u_curr` 做初始 halo 交换。
2. 主线程对 `u_prev` 做初始 halo 交换。
3. 主线程发起“初始能量阶段”。
4. 所有从线程各自计算自己的 `partial_energy`。
5. 每个组主线程做本组能量小归约。
6. 主线程再做 MPI 全局归约，得到初始能量。

### 6.2 每一步的推进顺序

对第 `step` 步，执行顺序是：

1. 如果当前 `u_curr` 的 halo 尚未就绪，主线程发起非阻塞 halo 交换。
2. 主线程发起“内部行阶段”。
3. 所有从线程只计算不依赖远端 halo 的内部行。
4. 主线程等待 halo 通信完成。
5. 主线程发起“边界行阶段”。
6. 所有从线程计算依赖 halo 的边界行。
7. 主线程等待所有组完成后执行 `swap_fields()`。
8. 如果当前步需要输出能量诊断：主线程先把新的 `u_curr` halo 补齐，再发起“能量阶段”。
9. 从线程计算各自 tile 的 `partial_energy`。
10. 组主线程做组内小归约，主线程再做 MPI 全局归约并输出日志。

## 7. 同步状态设计

每个时间步现在最多有三个阶段：

```c
initial_energy_state() = 1
compute_phase_state(step) = 3 * step + 2
boundary_phase_state(step) = 3 * step + 3
energy_phase_state(step) = 3 * step + 4
```

含义：

- `1`：初始能量阶段
- `3*step+2`：内部行更新阶段
- `3*step+3`：边界行更新阶段
- `3*step+4`：能量阶段

如果某一步不需要输出能量，就不会触发这一轮能量阶段。

## 8. 离散总能量如何归约

现在的能量归约分两层：

1. 从线程先计算本线程 tile 的 `partial_energy`。
2. 组主线程把本组所有从线程的 `partial_energy` 加起来。
3. 主线程只需要把各个组的能量再加起来，然后做一次 `MPI_Allreduce`。

因此主线程不再需要遍历所有从线程的任务结构去做本地归约。

## 9. 数值更新公式

更新核使用二维显式二阶差分：

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
- MPI 之间只交换 `Y` 向 halo。

## 10. 分辨率与计算量

代码提供两种模式：

```c
#define USE_FIXED_DOMAIN 1
```

### 10.1 固定物理区域，提升分辨率

当 `USE_FIXED_DOMAIN=1` 时：

- `LX / LY` 固定。
- 增大 `NX / NY` 会减小 `DX / DY`。
- 这表示在同一物理区域上做网格加密。

### 10.2 固定网格步长，扩大物理区域

当 `USE_FIXED_DOMAIN=0` 时：

- `DX / DY` 固定。
- 增大 `NX / NY` 会扩大物理区域。
- 这表示在同样空间分辨率下求解更大的区域。

## 11. CFL 稳定性检查

显式格式需要满足：

```text
(c*dt/dx)^2 + (c*dt/dy)^2 <= 1
```

程序启动时会检查：

- `CFL_X = C0 * DT / DX`
- `CFL_Y = C0 * DT / DY`
- `CFL_SUM2 = CFL_X^2 + CFL_Y^2`

若条件不满足，程序会直接报错退出。

## 12. 推荐阅读顺序

建议按下面顺序阅读代码：

1. `main_thread()`：看主线程如何只保留通信和调度。
2. `worker_thread()`：看从线程如何承担内部行、边界行和能量三类计算。
3. `group_main_thread()`：看组主线程如何做同步和组内能量归约。
4. `begin_y_halo_exchange_for()` / `end_y_halo_exchange_for()`：看 halo 通信如何改成非阻塞。
5. `compute_interior_block()` / `compute_boundary_block()`：看内部行与边界行的划分。
6. `compute_energy_block()`：看离散总能量如何按 tile 分配。

## 13. 关键代码位置

关键位置都在 [`wave_propagation_ghost.c`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c)：

- 进程域划分与初始化：[`wave_propagation_ghost.c:183`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c#L183)
- 非阻塞 halo 交换：[`wave_propagation_ghost.c:283`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c#L283)
- 内部行和边界行更新：[`wave_propagation_ghost.c:347`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c#L347)
- 从线程局部能量：[`wave_propagation_ghost.c:377`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c#L377)
- 组内能量归约：[`wave_propagation_ghost.c:430`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c#L430)
- 从线程入口：[`wave_propagation_ghost.c:492`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c#L492)
- 组主线程入口：[`wave_propagation_ghost.c:520`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c#L520)
- 主线程入口：[`wave_propagation_ghost.c:552`](/C:/Users/write/Documents/projects/parallelFrame/parallelFrame/wave_propagation_ghost.c#L552)

## 14. 编译与运行

### 14.1 编译

```bash
mpicc -Wall -O2 -o wave_propagation_ghost \
    wave_propagation_ghost.c \
    mythread/mythread.c \
    -lpthread -lm
```

### 14.2 运行

```bash
mpirun -np 4 ./wave_propagation_ghost
```

## 15. 说明

当前版本优先追求的是：

- 主线程尽量轻量
- 通信尽量和从线程计算重叠
- 诊断开销按需触发
- 后续还容易继续做性能调优

如果后面还要继续压主线程瓶颈，可以沿着这几个方向再往下做：

- 把能量诊断进一步改成异步采样或延后采样
- 把边界打包和更多边界处理也交给工作线程
- 继续细化 tile 粒度，改善组内负载均衡

## 16. 本次修改总结

相对最初那版 `wave_propagation_ghost` 示例，当前已经完成了下面这些关键修改：

### 16.1 并行结构重做

- 不再让每个线程维护一份独立波场和线程级 ghost 区。
- 改成“每个 MPI 进程一份共享局部网格”的结构。
- 保留 `MPI(Y) + Group(Y) + Worker(X)` 三层分解，但只让从线程拥有真实计算 tile。
- 组主线程不再占有自己的计算区块。

### 16.2 物理与数值一致性修正

- 修正了上下物理边界行的处理，确保全局最上和最下边界固定为 0。
- 初值构造时不再把高斯波包写到物理边界行。
- 能量输出从早期的 `sum(u^2)` 改成离散总能量诊断。
- 增加了 `USE_FIXED_DOMAIN` 选项，用来区分“固定物理区域的网格加密”和“固定步长的区域放大”。
- 增加了 CFL 稳定性检查，启动时就能拦截不稳定参数。

### 16.3 主线程和组主线程减负

- 组主线程不再直接调用差分更新核。
- 主线程不再做整块局部能量扫描。
- 本地能量先由从线程分块计算，再由组主线程做组内归约，最后才由主线程做跨组汇总。
- 主线程里整平面 `memset(u_next)` 也已经去掉，避免额外的大内存带宽开销。

### 16.4 计算进一步下放到从线程

- 从线程现在承担三类实际计算：
  1. 内部行更新
  2. 边界行更新
  3. 在需要时计算本线程 tile 的离散总能量
- 从线程已经成为主要的数值计算执行者，主线程和组主线程更像调度层。

### 16.5 通信与计算重叠

- halo 交换从阻塞模式扩展成了非阻塞 `begin/end` 两段式。
- 每步更新被拆成“内部行阶段”和“边界行阶段”。
- 主线程在 halo 通信进行时，先让从线程计算不依赖远端 halo 的内部行。
- 等通信完成后，再补算依赖 halo 的边界行。
- 这样可以把一部分 MPI 等待时间隐藏到从线程计算里。

### 16.6 诊断频率优化

- 能量诊断不再每步都做。
- 现在默认只在初始、每 `ENERGY_REPORT_INTERVAL` 步以及最后一步输出。
- 这显著减少了每步都触发的能量阶段和 `MPI_Allreduce` 次数。

### 16.7 文档与实现对齐

- 说明文档已经改成和当前代码结构一致。
- 文档中现在明确区分了：
  - 主线程职责
  - 组主线程职责
  - 从线程职责
  - 非阻塞 halo 交换
  - 内部行/边界行双阶段推进
  - 只按需做能量诊断
