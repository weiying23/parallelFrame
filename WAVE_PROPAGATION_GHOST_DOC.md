# 二维波动方程并行示例说明

## 1. 总览

本目录包含两个基于 `MPI + mythread` 的二维波动方程求解器，计算：

```
u_tt = C0² (u_xx + u_yy)
```

| 文件 | 调度方式 | 说明 |
|------|----------|------|
| `wave_propagation_ghost.c` | **任务池动态调度** | Worker 通过 `mt_taskpool` 竞争 2D chunk；适合 UMA / 负载不均场景 |
| `wave_propagation_ghost_fix.c` | **固定静态划分** | Worker 持有固定矩形 tile；适合多 NUMA / 核心数均匀场景 |

两个实现共享同一底层架构（三级 mythread 模块、配置文件、诊断输出）。

---

## 2. 架构分层

```
┌─────────────────────────────────────────┐
│  应用层 (wave_propagation_ghost*.c)      │
├─────────────────────────────────────────┤
│  Layer 3: Halo 交换 (mythread_halo)      │
│  · mythread_halo_exchange_intra (组间)   │
│  · mythread_halo_exchange_mpi  (MPI)     │
├─────────────────────────────────────────┤
│  Layer 2: NUMA 感知分配 (mythread_field) │
│  · GroupField (每组独立波场平面)          │
│  · group_alloc_field / GFIDX / swap      │
├─────────────────────────────────────────┤
│  Layer 1: 域分解 (mythread_decomp)       │
│  · Y_ONLY / XY_2D 两种策略               │
│  · group_tiles / worker_tiles / 邻居拓扑  │
├─────────────────────────────────────────┤
│  mythread 核心 (sync / pool / thread)    │
└─────────────────────────────────────────┘
```

### 2.1 NUMA 感知内存

每个线程组拥有一份独立的 `GroupField`，包含 3 个波场平面和 halo 缓冲。物理内存在 GMT 线程中分配（`bindcpu` 已生效），通过 `numa_alloc_onnode` 绑定到该组所在 NUMA 节点。无 libnuma 时退化为 `malloc` + first-touch。

### 2.2 组间 halo 交换

MMT 在每步计算前集中执行：
1. `apply_dirichlet_all_groups()` — 物理边界归零
2. `mythread_halo_exchange_intra()` — 组间 memcpy
3. `mythread_halo_exchange_mpi()` — 边界组 MPI 收发
4. `apply_dirichlet_all_groups()` — 再次归零

仅域边界组参与 MPI 通信，内部组完全通过本地 memcpy 完成。

---

## 3. 并行分解

### 3.1 三层分解

```
全局网格 (NX × NY)
  → MPI rank (Px × Py 二维进程网格)
    → 线程组 (Gx × Gy，Y_ONLY 或 XY_2D)
      → Worker 线程 (静态 tile 或 动态 2D chunk)
```

### 3.2 组间分解策略

| 策略 | 宏值 | 说明 |
|------|------|------|
| `Y_ONLY` | 0 | 仅 Y 方向切分组，每组覆盖进程全部 X 范围；组间邻居最多 2 个 |
| `XY_2D` | 1 | X+Y 都切分，组排为 2D 网格；组间邻居最多 4 个；MPI 通信仅由边界组承担 |

配置方式：`config/hardware.cfg` 中 `GROUP_DECOMP = 0` 或 `1`。

### 3.3 两种调度方式

#### 任务池 (`wave_propagation_ghost.c`)

- GMT 将本组子域按 `CHUNK_ROWS × CHUNK_COLS`（默认 16×256）切分为 2D chunk
- 所有 chunk 作为 `RowTaskCtx` 提交到 `mt_taskpool`
- 组内 worker 动态竞争取任务，快核多抢，慢核少做
- 天然负载均衡，适合 P/E 混合架构

#### 固定划分 (`wave_propagation_ghost_fix.c`)

- 每个 worker 在组内持有固定矩形 tile `[x_begin,x_end) × [y_begin,y_end)`
- 通过 `mythread_decomp_worker_tile` 均匀划分
- 零任务池开销（无 mutex 竞争）
- 配合 NUMA 分配在多 NUMA 机器上效果最优

---

## 4. 主要数据结构

### 4.1 `GroupField` (mythread_field.h)

```c
typedef struct GroupField {
  double *u_prev, *u_curr, *u_next;  // 三个波场平面
  int ny_padded, nx_padded, stride;  // 含 halo 的尺寸
  size_t plane_bytes;                // 单平面字节数
  double *send_to_up/down/left/right;  // 组间 halo 缓冲
  double *send_up/down/left/right;     // MPI halo 缓冲
  int numa_node;                     // NUMA 节点 ID
  double group_energy;               // 本组能量累加器
} GroupField;
```

### 4.2 `RowTaskCtx` (仅任务池版)

```c
typedef struct {
  int y_begin, y_end, x_begin, x_end;  // chunk 在 GroupField 内的坐标
  int gid;                              // 所属组 ID
  GroupField *gf;                       // 所属 GroupField
  double *energy_acc, *l2_acc;         // 原子累加器指针
  double *max_amp;                      // 最大振幅 CAS 指针
  int *max_amp_gx, *max_amp_gy;        // 振幅位置
} RowTaskCtx;
```

### 4.3 `SimulationData`

进程级 MPI 拓扑信息（两个版本共用）：

```c
typedef struct {
  int local_x_begin, local_x_end, local_nx;
  int local_y_begin, local_y_end, local_ny;
  int mpi_rank, mpi_size;
  int proc_x, proc_y, proc_px, proc_py;
  int neighbor_left, neighbor_right, neighbor_up, neighbor_down;
  double initial_energy;
} SimulationData;
```

波场平面不再存放在这里，已迁移到 `GroupField` 中。

---

## 5. 线程职责

### 5.1 主管理线程 MMT (`main_thread`)

- 发起阶段同步 (`mSetGrps` / `mSetSubs`)
- 组间 + MPI halo 交换
- 每组独立 `group_field_swap`
- MPI 全局能量/L² 归约
- 诊断输出

### 5.2 组主线程 GMT (`group_main_thread`)

- **分配本组 GroupField**（此时 bindcpu 已生效，NUMA 节点正确）
- 任务池版：管理 `mt_taskpool` 生命周期（attach/begin/submit/close/wait/shutdown）
- 固定划分版：`gSetSubs` → `gWaitSubs` 同步 worker
- 组内能量/L² 归约 (`reduce_group_worker_energy`)

### 5.3 Worker 线程 (`worker_thread`)

- 任务池版：`mt_taskpool_worker_loop` 循环取任务执行
- 固定划分版：等待 GMT 同步信号，在自己的固定 tile 上计算
- 四种回调：`task_init_field` / `task_compute_interior` / `task_compute_boundary` / `task_compute_energy`

---

## 6. 时间推进流程

### 6.1 初始化阶段

1. MMT 发起 `init_fields_state` 同步
2. GMT 分配 GroupField → Worker 并行初始化波场
3. MMT: `group_field_link_buffers` + Dirichlet + halo 交换 + `copy_curr_to_prev`
4. MMT 发起 `initial_energy_state` → Worker 计算能量/L²/振幅 → 归约输出

### 6.2 每步循环

```
for step in 0..NT-1:
  1. MMT: apply_dirichlet + halo_exchange(intra + mpi)
  2. MMT: start_phase(compute) → workers compute_interior → wait_phase
  3. MMT: start_phase(boundary) → workers compute_boundary → wait_phase
  4. MMT: group_field_swap (每组独立)
  5. if need_energy:
       MMT: apply_dirichlet + halo_exchange
       MMT: start_phase(energy) → workers compute_energy → reduce → MPI_Allreduce
```

### 6.3 同步状态编码

```
init_fields_state     = 1
initial_energy_state  = 2
compute_phase(s)      = 3s + 3
boundary_phase(s)     = 3s + 4
energy_phase(s)       = 3s + 5
```

---

## 7. 诊断指标

每步能量诊断输出包含：

| 指标 | 含义 | 守恒性 |
|------|------|--------|
| `E` | 总能量 ½∫(u_t² + C0²\|∇u\|²) | Dirichlet 边界下衰减 |
| `L2` | L² 范数 √(∫u² dA) | 随波扩散递增 |
| `max\|u\|` | 最大振幅 | 不守恒，检测 blow-up |
| `@(gx,gy)` | 振幅峰值全局坐标 | 追踪波前位置 |

示例输出：
```
[Main] Initial: E=1.570796 L2=148.875 max|u|=9.999996@(6999,6999)
[Main] Step   60/480, time 6.493,  E=1.570796 L2=257.860 max|u|=9.999996@(6999,6999)
```

---

## 8. 配置文件

参数通过两个配置文件读入（编译默认值 → 配置文件 → 环境变量，优先级递增）：

### `config/case.cfg` — 算例参数

```ini
NX = 14000
NY = 14000
NT = 480
A = 10.0
DT = 0.01
C0 = 0.1
USE_FIXED_DOMAIN = 0
HALO = 1
ENERGY_REPORT_INTERVAL = 60
```

### `config/hardware.cfg` — 硬件/线程参数

```ini
N_GROUPS = 2
N_WORKERS = 3
GROUP_DECOMP = 0
NCorePClu = 5
NCluPNode = 2
NCorePGrp = 4
ManageCoreId = 4
```

### 环境变量覆盖

```bash
# 自定义配置路径
WAVE_CASE_CFG=my_case.cfg WAVE_HARDWARE_CFG=my_hw.cfg mpirun ...

# 负载不均衡测试
WAVE_IMBALANCE_PCT=50 mpirun -np 1 ./wave_propagation_ghost
WAVE_IMBALANCE_WORKER=0 WAVE_IMBALANCE_PCT=30 mpirun ...
```

---

## 9. 编译与运行

### 9.1 编译

使用 CMake 构建（默认启用 MPI）：

```bash
cmake -S . -B build
cmake --build build -j
```

需要指定 MPI 编译器时（如华为 HPCKit）：

```bash
cmake -S . -B build -DCMAKE_C_COMPILER=$MPI_HOME/bin/mpicc
cmake --build build -j
```

若不需要 MPI，可关闭 MPI 只构建非 MPI 目标：

```bash
cmake -S . -B build -DUSE_MPI=OFF
cmake --build build -j
```

生成目标：`wave_propagation_ghost`（任务池版）、`wave_propagation_ghost_fix`（固定划分版），及全部单测程序。

编译依赖 `mpicc`、`libpthread`、`libm`。可选 `libnuma`（Linux 上 NUMA 感知分配）。

### 9.2 运行

```bash
# 默认参数（从 config/ 读取）
mpirun -np 4 ./build/wave_propagation_ghost
mpirun -np 4 ./build/wave_propagation_ghost_fix

# 不均衡负载测试
WAVE_IMBALANCE_PCT=50 mpirun -np 1 ./build/wave_propagation_ghost
```

### 9.3 性能测试

```bash
# 简易对比（3 次取中位数，固定 1000x1000）
./scripts/benchmark.sh

# 完整矩阵（quick=3 网格 / full=5 网格，各 5 次）
./scripts/benchmark_full.sh quick
./scripts/benchmark_full.sh full
```

日志保存至 `logs/YYYYMMDD_HHMMSS/`，汇总 CSV 自动生成。

---

## 10. 数值方法

### 10.1 离散格式

```
u_next(y,x) = 2*u_curr(y,x) - u_prev(y,x)
  + C0²·dt²·[(u_curr(y,x-1)-2·u_curr(y,x)+u_curr(y,x+1))/dx²
            + (u_curr(y-1,x)-2·u_curr(y,x)+u_curr(y+1,x))/dy²]
```

### 10.2 边界条件

- 全局 X=0、X=NX-1、Y=0、Y=NY-1：Dirichlet 固定为 0
- MPI 进程间：halo 行/列交换
- 组间：本地 memcpy 交换边界行/列

### 10.3 稳定性

显式格式需满足 CFL 条件：`(C0·dt/dx)² + (C0·dt/dy)² ≤ 1`。程序启动时自动检查。

---

## 11. mythread 模块清单

| 模块 | 文件 | 功能 |
|------|------|------|
| 域分解 | `mythread_decomp.h/.c` | Y_ONLY/XY_2D 策略，tile 分配，邻居拓扑 |
| 字段分配 | `mythread_field.h/.c` | GroupField，NUMA 感知分配，GFIDX/swap |
| Halo 交换 | `mythread_halo.h/.c` | 组内 memcpy 交换 + MPI 边界交换 |
| 配置解析 | `mythread_config.h/.c` | INI 格式配置文件解析 |
| 同步原语 | `mythread_sync.h/.c` | 三级屏障 (SM/MS, SG/GS, GM/MG) |
| 任务池 | `mythread_pool.h/.c` | SPMC 任务队列 (epoch-driven) |
| 线程管理 | `mythread_thread.h/.c` | 线程创建/生命周期/CPU 绑定 |
| 计时器 | `mythread_timer.h/.c` | TSC 性能计时 |
| LocV 存储 | `mythread_locv.h/.c` | 线程/组/进程三级指针存储 |
| Fortran 接口 | `mythread_fortran.h/.c` | trailing-underscore 包装 |
| CPU 绑定测试 | `test_bindcpu.c` | 9 个测试用例 |
| 域分解测试 | `test_field.c` | 10 个 GroupField + Halo 测试 |

---

## 12. 性能调优参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `CHUNK_ROWS` | 16 | 任务池版 Y 方向 chunk 行数 |
| `CHUNK_COLS` | 256 | 任务池版 X 方向 chunk 列数；设极大值退化为 Y-only |
| `TP_CAP` | 512 | 任务池容量；BLOCK 模式下 PMT 阻塞等待 |
| `N_GROUPS` | 2 | 线程组数；≤1 自动切换非分组模式 |
| `N_WORKERS` | 3 | 每组 worker 数（不含 GMT） |
| `GROUP_DECOMP` | 0 | 0=Y_ONLY, 1=XY_2D |

---

## 13. 已知限制

- `choose_2d_grid` 中 `span_x/y` 为 `int`，超大网格（>2³¹ 点）可能溢出
- macOS 不支持 `sched_setaffinity`，CPU 绑定退化为空操作
- 无 libnuma 时 NUMA 感知退化为 malloc + first-touch（仍有效但无显式节点绑定）
- 任务池版 `g_n_flat_tasks` 仍为 `int`（非组模式），极端大网格下需改为 `size_t`
- `apply_imbalance_to_group` 仅修改 Y 方向，X 方向始终保持均匀
