# NUMA 感知的线程组独立内存分配方案

## 一、背景与动机

### 1.1 当前问题

```
NUMA 0                          NUMA 1
┌──────────────────────┐        ┌──────────────────────┐
│ group 0 线程 (绑核)    │        │ group 1 线程 (绑核)    │
│ 访问 ──────────────→  │ remote │  ←────────────── 访问 │
│                      │ ←────→ │                      │
└──────────────────────┘        └──────────────────────┘
         ↓                               ↓
┌──────────────────────────────────────────────────────┐
│  进程唯一网格 (u_prev/u_curr/u_next)                     │
│  malloc 在首次 touch 的 NUMA 节点上分配物理页             │
│  实际物理内存只在 NUMA 0 或 NUMA 1，取决于谁先写           │
└──────────────────────────────────────────────────────┘
```

- 进程内网格数组是一整块 `malloc`/`calloc` 分配的内存
- 物理页面由 first-touch 策略决定落在哪个 NUMA 节点
- 若数据全在 NUMA 0，NUMA 1 上的线程组全部远端访存（~1.5-2× 延迟惩罚）
- mythread 已保证每个线程组的线程绑定到同一 NUMA 节点的核心上，但**数据不随组迁移**

### 1.2 前提条件

| 条件 | 状态 |
|------|------|
| 线程组内的线程绑定到同一 NUMA 节点 | mythread 的 `initmd()` 已保证（通过 `indg` 连续分配） |
| 线程组间无计算依赖 | 组间仅 halo 行有数据依赖，其他区域完全并行 |
| 每组在 Y 方向覆盖进程子域的一个连续子区间 | `split_range` 已保证 |
| 每组在 X 方向的覆盖范围 | 可配置：全部覆盖（Y_ONLY）或子区间覆盖（XY_2D），见 §1.4 |

### 1.3 设计目标

1. **每组独立分配**：每个线程组拥有自己的 3 平面波场数组 + halo 缓冲
2. **NUMA 本地化**：每组的物理内存在该组绑定的 NUMA 节点上分配
3. **组间 halo 交换**：相邻组之间通过 MMT 协调的显式数据拷贝完成
4. **MPI 通信收拢**：仅边界组参与 MPI 收发，内部组的组间 halo 走本地内存拷贝
5. **代码向后兼容**：非 NUMA 平台（macOS / 未装 libnuma）降级为原始单一分配

### 1.4 组间分解策略：Y_ONLY vs XY_2D

"每组是否覆盖 X 方向的全部进程子域" 不是硬约束，而是可配置的设计选择。

**方案 Y_ONLY（仅 Y 方向切分，默认）**：

```
进程本地网格
┌──────────────────────────────────────┐  ← 完整 X 范围 (local_nx)
│  group 0   y: [HALO, y1)             │
├──────────────────────────────────────┤
│  group 1   y: [y1,   y2)             │
├──────────────────────────────────────┤
│  group 2   y: [y2,   ny+HALO)        │
└──────────────────────────────────────┘
```

- 每组覆盖进程 **X 方向全部**，Y 方向一个连续子区间
- 组间邻居：**最多 2 个**（上 Y / 下 Y）
- 组间 halo：**2 次 memcpy / 组对**（仅 Y 方向）
- 所有组都接触 X 边界 → 每个组都要参与 X 方向 MPI halo

**方案 XY_2D（X 和 Y 方向都切分）**：

```
进程本地网格
┌──────────────┬──────────────┐
│  group 0     │  group 1     │  ← Y 区间 A
│  x: [0, x1)  │  x: [x1, nx) │
├──────────────┼──────────────┤
│  group 2     │  group 3     │  ← Y 区间 B
│  x: [0, x1)  │  x: [x1, nx) │
└──────────────┴──────────────┘
```

- 每组覆盖进程 **X + Y 方向的一个矩形子域**
- 组间邻居：**最多 4 个**（上/下/左/右）
- 组间 halo：**4 次 memcpy / 组对**（X + Y 方向，以及角点一致性）
- 仅边界组参与 MPI → MPI 通信总量更少

**选择依据**：

| 条件 | 推荐策略 |
|------|----------|
| 进程本地 `local_nx` 不大（< ~2000），或组数少（≤4） | Y_ONLY（halo 逻辑简单） |
| `local_nx` 很大，Y 区间短 → 组间负载不均衡 | XY_2D（按 X 也切分） |
| NUMA 节点间带宽充足 | Y_ONLY 够用 |
| NUMA 节点间带宽紧张（> 2 组） | XY_2D 优先（减少跨 NUMA 通信面） |

**设计原则**：两种策略使用**相同的 `GroupField` 结构体和代码路径**，仅在 `setup_thread_tasks` 中根据配置计算不同的 `x_begin_global/x_end_global` 和组邻居拓扑。halo 交换函数通过遍历邻居表统一处理，避免 if/else 分支蔓延。



---

## 二、数据结构设计

### 2.1 组级网格字段 `GroupField`

```c
typedef struct {
  /* ── 波场平面 ── */
  double *u_prev;              // (ny+2*HALO) * (nx+2*HALO)
  double *u_curr;
  double *u_next;

  /* ── 尺寸 ── */
  int ny_padded;               // = ny_group + 2*HALO（本组含 halo 的总行数）
  int nx_padded;               // = nx_local  + 2*HALO（本组含 halo 的总列数）
  int stride;                  // = nx_padded（idx 宏用）
  int ny_interior;             // = ny_group（内部计算行数）

  /* ── 全局坐标映射 ── */
  int y_begin_global;          // 本组不含 halo 的第一行在全局网格中的 Y 坐标
  int y_end_global;            // 本组不含 halo 的最后一行+1 的全局 Y 坐标
  int x_begin_global;          // 本组不含 halo 的第一列在全局网格中的 X 坐标
  int x_end_global;

  /* ── 组间 halo 缓冲（Y 方向）── */
  double *send_to_up;          // 向上组发送：本组最上一行内部数据（nx_local 个 double）
  double *recv_from_up;        // 从上级接收：上级组的最下一行内部数据
  double *send_to_down;        // 向下组发送：本组最下一行内部数据
  double *recv_from_down;

  /* ── 组间 halo 缓冲（X 方向，仅 XY_2D 模式使用）── */
  double *send_to_left;        // 向左组发送：本组最左一列（ny_group 个 double）
  double *recv_from_left;
  double *send_to_right;       // 向右组发送：本组最右一列
  double *recv_from_right;

  /* ── MPI halo 缓冲（X / Y 边界）── */
  double *send_up;             // MPI 邻居 halo（仅边界组使用）
  double *recv_up;
  double *send_down;
  double *recv_down;
  double *send_left;
  double *recv_left;
  double *send_right;
  double *recv_right;

  /* ── NUMA 信息 ── */
  int numa_node;               // 分配所在的 NUMA 节点
  int is_boundary_top;         // 本组是否为进程的 Y 上边界组
  int is_boundary_bottom;      // 本组是否为进程的 Y 下边界组
  int neighbor_group_up;       // 相邻上组 ID（-1 表示无）
  int neighbor_group_down;     // 相邻下组 ID（-1 表示无）
  int neighbor_group_left;     // 相邻左组 ID（-1 表示无，仅 XY_2D）
  int neighbor_group_right;    // 相邻右组 ID（-1 表示无，仅 XY_2D）

  /* ── 分解模式 ── */
  int decomp_y_only;           // 1 = Y_ONLY, 0 = XY_2D

  /* ── 能量累加器 ── */
  double group_energy;         // 本组内部能量（本组 worker 累加）
} GroupField;

/* 全局注册表，MMT 持有，所有组的 GMT 可读 */
static GroupField *g_gfields = NULL;  /* [md.ngrp] */
```

### 2.2 `ThreadTask` 扩展

```c
typedef struct {
  int gid, tid, cpu_id;
  int x_begin, x_end;          // 相对于 GroupField 内部坐标（以 HALO 为原点）
  int y_begin, y_end;
  GroupField *gf;              // 【新增】直接指针，避免每次通过 gid 索引
  double partial_energy;
  /* ... 计时字段不变 ... */
} ThreadTask;
```

### 2.3 进程全局结构精简

分配变更为组级后，进程全局只保留：

```c
static SimulationData g_sim = {0};     // MPI 拓扑信息（保留，不再持有波场）
static GroupField *g_gfields = NULL;   // 组级网格数组
static double *g_mpi_x_send = NULL;    // 【可选】X 方向 MPI 收拢缓冲
static double *g_mpi_x_recv = NULL;
```

移除的全局变量：
- `g_sim.u_prev / u_curr / u_next` → 迁移到 `GroupField`
- `g_send_up / down / left / right` → 迁移到边界 `GroupField`
- `g_group_energy` → 替换为 `GroupField.group_energy`

---

## 三、内存分配策略

### 3.1 分配时机

```
main()
  ├─ init_simulation(mpi_rank, mpi_size)   ← 仅初始化拓扑、MPI 邻居
  ├─ InitThreads(...)
  ├─ setup_thread_tasks()                   ← 确定每组 Y 范围
  ├─ StartThreads(thread_run)
  │
  └─ MMT: thread_run() → main_thread()
       ├─ g_gfields = calloc(ngrp, sizeof(GroupField))   ← 分配注册表
       ├─ mSetGrps(init_fields_state)
       │
       └─ 各 GMT 并行执行：
            group_main_thread()
              ├─ 【关键】此时 bindcpu 已完成，线程在目标 NUMA 节点上
              ├─ group_alloc_field(gf, ny, nx)
              │    ├─ numa_node = numa_node_of_cpu(sched_getcpu())
              │    ├─ gf->u_prev = numa_alloc_onnode(plane*8, node)
              │    ├─ gf->u_curr = numa_alloc_onnode(plane*8, node)
              │    ├─ gf->u_next = numa_alloc_onnode(plane*8, node)
              │    └─ halo bufs = numa_alloc_onnode(...)
              └─ 初始化本组 u_curr / u_prev（first-touch 在同一 NUMA）
```

### 3.2 分配函数

```c
static int group_alloc_field(GroupField *gf, int ny_group, int nx_local,
                              int y_begin_global, int y_end_global,
                              int x_begin_global, int x_end_global,
                              int neighbor_up, int neighbor_down) {
  size_t plane = (size_t)(ny_group + 2 * HALO) * (size_t)(nx_local + 2 * HALO);
  int node = numa_node_of_cpu(sched_getcpu());

  gf->ny_padded  = ny_group + 2 * HALO;
  gf->nx_padded  = nx_local + 2 * HALO;
  gf->stride     = gf->nx_padded;
  gf->ny_interior = ny_group;
  gf->y_begin_global = y_begin_global;
  gf->y_end_global   = y_end_global;
  gf->x_begin_global = x_begin_global;
  gf->x_end_global   = x_end_global;
  gf->neighbor_group_up   = neighbor_up;
  gf->neighbor_group_down = neighbor_down;
  gf->numa_node  = node;

  gf->u_prev = numa_alloc_onnode(plane * sizeof(double), node);
  gf->u_curr = numa_alloc_onnode(plane * sizeof(double), node);
  gf->u_next = numa_alloc_onnode(plane * sizeof(double), node);

  /* 组间 halo */
  gf->send_to_up   = numa_alloc_onnode(nx_local * sizeof(double), node);
  gf->recv_from_up = numa_alloc_onnode(nx_local * sizeof(double), node);
  gf->send_to_down = numa_alloc_onnode(nx_local * sizeof(double), node);
  gf->recv_from_down = numa_alloc_onnode(nx_local * sizeof(double), node);

  /* MPI halo — 仅边界组需要分配，内部组指针为 NULL */
  if (gf->neighbor_group_up < 0) {
    gf->send_up = numa_alloc_onnode(nx_local * sizeof(double), node);
    gf->recv_up = numa_alloc_onnode(nx_local * sizeof(double), node);
  }
  if (gf->neighbor_group_down < 0) {
    gf->send_down = numa_alloc_onnode(nx_local * sizeof(double), node);
    gf->recv_down = numa_alloc_onnode(nx_local * sizeof(double), node);
  }
  /* X 方向所有组都可能有 MPI 邻居（NX 边界） */
  /* ... send_left/right/recv_left/right 同理 ... */

  return (gf->u_prev && gf->u_curr && gf->u_next) ? 0 : -1;
}
```

### 3.3 降级路径（无 libnuma）

```c
#ifdef HAS_LIBNUMA
#include <numa.h>
#define numa_alloc(sz, node) numa_alloc_onnode((sz), (node))
#else
#define numa_alloc(sz, node) hmalloc(sz)
#define numa_node_of_cpu(cpu) (-1)
#endif
```

macOS 或不装 libnuma 时退化为普通 `malloc`。first-touch 策略仍然有效（GMT 线程在初始化时首次写数据），但无法显式绑定物理节点。

### 3.4 内存在线量估算

| 组件 | 原方案 | 新方案 | 增量 |
|------|--------|--------|------|
| 波场 3 平面 | `(ny+2)*(nx+2)*8*3` | 每组 `(ny_g+2)*(nx+2)*8*3`，求和 | 重叠 halo 行：`ngrp * 2*HALO * nx * 8 * 3` |
| 组间 halo 缓冲 | 无 | 每组 `nx * 8 * 4` | 新增 |
| MPI halo 缓冲 | `nx*8*4 + (ny+2)*8*4` | 同左（仅边界组） | 减少（内部组不分配） |
| GroupField 元数据 | 无 | `sizeof(GroupField) * ngrp` (~200B×ngrp) | 可忽略 |

**示例**：`ngrp=4, nx=1000, ny=250/组, HALO=1`
- 原方案：`(1000+2)*(250*4+2)*8*3 ≈ 72 MB`
- 新方案：每组 `(252)*(1002)*8*3 ≈ 6.06 MB`，4 组合计 `≈ 24.2 MB`
  - 重叠 halo 增量：`4 * 2 * 1000 * 8 * 3 = 192 KB`
- 新方案实际更省内存（因为每组的 ny 更小，少了大量行）

**注意**：总内存不会翻倍。原方案 NY 行是全局尺寸，新方案每组 ny 是子区间尺寸，加上少量 halo 冗余。

---

## 四、线程-组-NUMA 映射关系

### 4.1 映射保证

mythread 的 `initmd()` 已按以下逻辑分配 `indg`：

```
进程级:
  ipr = mpi_rank % NProcPNode
  cbase = ipr * NCorePClu * NGrpPProc

组 i:
  组基址 = cbase + i * NCorePClu

  组内线程 j:
    indg = 组基址 + j + ib

  若组 i 包含 MMT: ib=1, Nthreads--
```

同一组内所有线程的 `indg` 落在 `[cbase+i*NCorePClu, cbase+(i+1)*NCorePClu)` 区间。只要保证 `NCorePClu` 不超过单个 NUMA 节点的核心数，即可保证组内所有线程在同一 NUMA 节点上。

### 4.2 NUMA 节点推断

```c
// GMT 在分配时调用
static int get_current_numa_node(void) {
#ifdef HAS_LIBNUMA
  int cpu = sched_getcpu();
  return numa_node_of_cpu(cpu);
#else
  return -1;
#endif
}
```

`numa_node_of_cpu()` 返回当前核心所在 NUMA 节点 ID。GMT 在 `group_main_thread()` 入口时已通过 `bindcpu(ti->indg)` 绑核，所以 `sched_getcpu()` 返回的正是目标核心。

---

## 五、Halo 交换重新设计

### 5.1 数据流概览

**Y_ONLY 模式**（组间仅 Y 方向邻居）：

```
  group 0 ───send_to_down──→ group 1
           ←──recv_from_down─
                          group 1 ───send_to_down──→ group 2
                                   ←──recv_from_down─

  MPI: group 0 ⇄ neighbor_down,  group N-1 ⇄ neighbor_up
       all groups ⇄ neighbor_left/right (X 方向)
```

**XY_2D 模式**（组间有 X + Y 方向邻居，以 2×2 为例）：

```
            send_to_right                    send_to_right
  group 0 ──────────────→ group 1    group 2 ──────────────→ group 3
           ←──────────────                   ←──────────────
            recv_from_left                    recv_from_left

       ↑  │                              ↑  │
  send_to_down                          send_to_down
       │  ↓                              │  ↓
       │  recv_from_up                   │  recv_from_up
  
  group 0 ⇵ group 2                   group 1 ⇵ group 3

  MPI: 仅边界组参与（上边界、下边界、左边界、右边界）
       内部组完全通过组间 memcpy 完成 halo 交换
```



### 5.2 同步协议（新增一个同步阶段）

原协议（每步 3 阶段）：compute → boundary → energy

新协议（每步 4 阶段）：

```
Step N:
  ┌──────────────────────────────────────────────────────┐
  │  1. MMT: mSetGrps(COMPUTE)                           │
  │     GMT: gWaitMain → 组间 halo 交换（send_to_down/up）│
  │          → gSetMain(COMPUTE)                          │
  │     MMT: mWaitGrps → 确认组间 halo 完毕               │
  ├──────────────────────────────────────────────────────┤
  │  2. MMT: 边界组 MPI halo 交换                        │
  │          → exchange_halos_for_boundary()              │
  ├──────────────────────────────────────────────────────┤
  │  3. MMT: mSetGrps(COMPUTE_INTERIOR)                  │
  │     GMT: gWaitMain → 提交 compute 任务 → gWaitSubs    │
  │          → gSetMain(COMPUTE_INTERIOR)                 │
  │     MMT: mWaitGrps                                    │
  ├──────────────────────────────────────────────────────┤
  │  4. MMT: mSetGrps(BOUNDARY)                          │
  │     GMT: 提交 boundary 任务 → gSetMain(BOUNDARY)      │
  │     MMT: mWaitGrps                                    │
  └──────────────────────────────────────────────────────┘
```

或者更简化的方案：将组间 halo 和 MPI halo 合并到 **MMT 统一执行**，GMT 不参与 halo 通信：

```
MMT 主循环每步：

  1. 对所有组：执行组间 halo 拷贝（memcpy 相邻组的 send/recv buf）
  2. 对边界组：执行 MPI halo 交换
  3. mSetGrps(COMPUTE) → 所有 GMT + workers 并行计算 → mWaitGrps
  4. mSetGrps(BOUNDARY) → 边界行计算 → mWaitGrps
  5. swap_fields（仅指针交换，数据不动）
```

这个简化方案更清晰，MMT 集中管理所有通信。

### 5.3 MMT 集中式 halo 交换

```c
static void intra_process_halo_exchange(void) {
  int g;

  /* ── 第一阶段：从源组拷贝数据到目标组的 recv 缓冲 ── */
  for (g = 0; g < md.ngrp; g++) {
    GroupField *gf = &g_gfields[g];
    int nx_int = gf->nx_padded - 2 * HALO;
    int ny_int = gf->ny_interior;

    /* Y 方向 — 所有模式都需要 */
    if (gf->neighbor_group_up >= 0) {
      GroupField *up = &g_gfields[gf->neighbor_group_up];
      memcpy(up->recv_from_down,
             &gf->u_curr[HALO * gf->stride + HALO],
             (size_t)nx_int * sizeof(double));
    }
    if (gf->neighbor_group_down >= 0) {
      GroupField *down = &g_gfields[gf->neighbor_group_down];
      memcpy(down->recv_from_up,
             &gf->u_curr[ny_int * gf->stride + HALO],
             (size_t)nx_int * sizeof(double));
    }

    /* X 方向 — 仅 XY_2D 模式 */
    if (!gf->decomp_y_only) {
      if (gf->neighbor_group_left >= 0) {
        GroupField *left = &g_gfields[gf->neighbor_group_left];
        for (int row = 0; row < ny_int; row++) {
          left->recv_from_right[row] = gf->u_curr[(HALO + row) * gf->stride + HALO];
        }
      }
      if (gf->neighbor_group_right >= 0) {
        GroupField *right = &g_gfields[gf->neighbor_group_right];
        int src_col = HALO + nx_int - 1;
        for (int row = 0; row < ny_int; row++) {
          right->recv_from_left[row] = gf->u_curr[(HALO + row) * gf->stride + src_col];
        }
      }
    }
  }

  /* ── 第二阶段：将接收缓冲拷贝到目标组的 halo 行/列 ── */
  for (g = 0; g < md.ngrp; g++) {
    GroupField *gf = &g_gfields[g];
    int nx_int = gf->nx_padded - 2 * HALO;
    int ny_int = gf->ny_interior;

    if (gf->neighbor_group_up >= 0) {
      memcpy(&gf->u_curr[0 * gf->stride + HALO],         // 上 halo 行
             gf->recv_from_up,
             (size_t)nx_int * sizeof(double));
    }
    if (gf->neighbor_group_down >= 0) {
      memcpy(&gf->u_curr[(ny_int + HALO) * gf->stride + HALO],  // 下 halo 行
             gf->recv_from_down,
             (size_t)nx_int * sizeof(double));
    }
    if (!gf->decomp_y_only) {
      if (gf->neighbor_group_left >= 0) {
        for (int row = 0; row < ny_int; row++) {
          gf->u_curr[(HALO + row) * gf->stride + 0] = gf->recv_from_left[row];
        }
      }
      if (gf->neighbor_group_right >= 0) {
        int dst_col = HALO + nx_int;
        for (int row = 0; row < ny_int; row++) {
          gf->u_curr[(HALO + row) * gf->stride + dst_col] = gf->recv_from_right[row];
        }
      }
    }
  }
}

static void boundary_mpi_halo_exchange(ThreadTask *task) {
  MPI_Request requests[8];
  int req_count = 0;

  /* Y 方向：只有最上和最下组参与 */
  GroupField *top = &g_gfields[0];
  GroupField *bot = &g_gfields[md.ngrp - 1];

  if (g_sim.neighbor_up >= 0) {
    int src_row = bot->ny_interior;
    int nx_int = bot->nx_padded - 2 * HALO;
    memcpy(bot->send_up, &bot->u_curr[src_row * bot->stride + HALO],
           (size_t)nx_int * sizeof(double));
    MPI_Isend(bot->send_up, nx_int, MPI_DOUBLE,
              g_sim.neighbor_up, 100, MPI_COMM_WORLD, &requests[req_count++]);
    MPI_Recv(bot->recv_up, nx_int, MPI_DOUBLE,
             g_sim.neighbor_up, 101, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    memcpy(&bot->u_curr[(bot->ny_interior + HALO) * bot->stride + HALO],
           bot->recv_up, (size_t)nx_int * sizeof(double));
  }
  if (g_sim.neighbor_down >= 0) {
    int nx_int = top->nx_padded - 2 * HALO;
    memcpy(top->send_down, &top->u_curr[HALO * top->stride + HALO],
           (size_t)nx_int * sizeof(double));
    MPI_Isend(top->send_down, nx_int, MPI_DOUBLE,
              g_sim.neighbor_down, 101, MPI_COMM_WORLD, &requests[req_count++]);
    MPI_Recv(top->recv_down, nx_int, MPI_DOUBLE,
             g_sim.neighbor_down, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    memcpy(&top->u_curr[0 * top->stride + HALO],
           top->recv_down, (size_t)nx_int * sizeof(double));
  }

  /* X 方向：Y_ONLY 模式所有组参与，XY_2D 模式仅边界组参与 */
  for (int g = 0; g < md.ngrp; g++) {
    GroupField *gf = &g_gfields[g];
    if (!gf->decomp_y_only &&
        gf->neighbor_group_left >= 0) continue;  // XY_2D 且内部组，跳过
    /* ... send_left/right MPI 通信，每列拷贝 ... */
    /* （X 方向 MPI halo 逻辑与 Y_ONLY 相同，仅遍历的组不同） */
  }

  MPI_Waitall(req_count, requests, MPI_STATUSES_IGNORE);
}
```

### 5.4 分解配置与邻居拓扑构建

```c
typedef enum { DECOMP_Y_ONLY = 0, DECOMP_XY_2D = 1 } GroupDecomp;

static GroupDecomp g_group_decomp = DECOMP_Y_ONLY;  // 默认 Y_ONLY

static void setup_group_decomposition(int ngroups) {
  if (g_group_decomp == DECOMP_Y_ONLY) {
    /* Y_ONLY 拓扑：线性链表 */
    for (int g = 0; g < ngroups; g++) {
      g_gfields[g].decomp_y_only = 1;
      g_gfields[g].neighbor_group_up   = (g > 0)          ? (g - 1) : -1;
      g_gfields[g].neighbor_group_down = (g < ngroups - 1) ? (g + 1) : -1;
      g_gfields[g].neighbor_group_left  = -1;  /* 不适用 */
      g_gfields[g].neighbor_group_right = -1;
    }
  } else {
    /* XY_2D 拓扑：2D 网格 */
    int gx = 1, gy = 1;
    choose_group_grid(ngroups, &gx, &gy);  // 复用现有 choose_2d_grid
    for (int g = 0; g < ngroups; g++) {
      int gix = g % gx, giy = g / gx;
      g_gfields[g].decomp_y_only = 0;
      g_gfields[g].neighbor_group_up   = (giy > 0)        ? (g - gx) : -1;
      g_gfields[g].neighbor_group_down = (giy < gy - 1)   ? (g + gx) : -1;
      g_gfields[g].neighbor_group_left  = (gix > 0)        ? (g - 1)  : -1;
      g_gfields[g].neighbor_group_right = (gix < gx - 1)   ? (g + 1)  : -1;
    }
  }
}
```

### 5.5 简化为：不拆阶段，并入现有 COMPUTE 阶段

为避免增加新的同步阶段，可以将组间 halo 放在 `start_phase_from_main(compute_state)` 之后、`wait_phase_from_main(compute_state)` 之前这段 MMT 独占时间里完成：

```c
// MMT main_thread() 主循环:
for (int step = 0; step < NT; step++) {
  start_phase_from_main(compute_state);

  /* MMT 独占：组间 halo + MPI halo */
  intra_process_halo_exchange();
  boundary_mpi_halo_exchange(task);

  wait_phase_from_main(compute_state);   // workers 只做计算
  // ...
}
```

这样 GMT 和 worker 不感知 halo 交换，接口不变。

---

## 六、坐标系统与索引变更

### 6.1 `idx()` 宏变更

```c
// 原版
#define idx(y, x) ((size_t)(y) * g_sim.stride + (size_t)(x))

// 新版 — 每组的 stride 不同
#define idx(gf, y, x) ((size_t)(y) * (size_t)(gf)->stride + (size_t)(x))
```

所有波场访问从 `g_sim.u_curr[idx(y, x)]` 变为 `gf->u_curr[idx(gf, y, x)]`。

### 6.2 全局↔局部坐标转换

```c
// 原版
static inline int global_y_from_local(int local_y) {
  return g_sim.local_y_begin + (local_y - HALO);
}

// 新版
static inline int global_y_from_local(const GroupField *gf, int local_y) {
  return gf->y_begin_global + (local_y - HALO);
}
```

### 6.3 任务上下文 `RowTaskCtx` 扩展

```c
typedef struct {
  int y_begin, y_end;
  int x_begin, x_end;
  GroupField *gf;       // 【新增】直接指向本组网格
  double *energy_acc;   // 指向 gf->group_energy
} RowTaskCtx;
```

计算回调中直接使用 `ctx->gf` 访问本组数据，无需全局查找。

---

## 七、各线程角色的职责变更

### 7.1 MMT（`main_thread`）

| 职责 | 原方案 | 新方案 |
|------|--------|--------|
| 分配波场 | `init_simulation` 时分配全局网格 | 不分配波场，分配 `g_gfields` 注册表 |
| halo 交换 | 每个时间步做一次 MPI halo | 组间 memcpy + 边界组 MPI halo |
| 能量汇总 | `accumulate_worker_energy` 遍历 worker | 遍历 `g_gfields[].group_energy` 累加 |

### 7.2 GMT（`group_main_thread`）

| 职责 | 原方案 | 新方案 |
|------|--------|--------|
| 分配本组内存 | 无 | `group_alloc_field()` |
| 初始化本组波场 | 通过 taskpool 提交 init 任务 | 同上，但数据写入本组 `GroupField` |
| 同步 | `gWaitMain → gSetSubs → gWaitSubs → gSetMain` | 不变 |
| 释放本组内存 | 无 | `group_free_field()` |

### 7.3 Worker（`worker_thread`）

| 职责 | 原方案 | 新方案 |
|------|--------|--------|
| 计算 | 通过 `g_sim.u_curr` 访问全局网格 | 通过 `ctx->gf->u_curr` 访问本组网格 |
| 坐标 | 全局 `y_begin` 换算 | 组内局部坐标 + `gf->y_begin_global` 换算 |

Worker 的代码变更量最大（所有 `g_sim.u_*` 访问点），但变更模式机械统一。

---

## 八、非组模式兼容

当 `ThreadG == 0`（非分组模式）时：

```c
// 只有一个逻辑"组"：整个进程域
g_gfields = calloc(1, sizeof(GroupField));
group_alloc_field(&g_gfields[0],
  g_sim.local_ny, g_sim.local_nx,
  g_sim.local_y_begin, g_sim.local_y_end,
  g_sim.local_x_begin, g_sim.local_x_end,
  -1, -1);  // 无相邻组

// 所有 worker 共享这一个 GroupField
// 组间 halo 不触发（neighbor_up/down == -1）
// MPI halo 正常执行
```

非组模式的内存布局与新版组模式对齐，代码路径统一。

---

## 九、实施计划

### 阶段一：引入 `GroupField`，不改计算路径（≈ 兼容过渡）

1. 定义 `GroupField` 结构体
2. `init_simulation` 改为仅分配 MPI 拓扑信息
3. 在 MMT 中分配 `g_gfields[1]`（单组），`GroupField` 指向原来的全局 `g_sim.u_*` 内存
4. 所有计算路径通过 `GroupField*` 间接访问，但物理内存仍是进程共享的一块
5. 验证数值正确性无回归

### 阶段二：真正的每组独立分配（组模式启用）

1. 移除全局 `g_sim.u_prev/curr/next`
2. GMT 在启动时调用 `group_alloc_field`
3. 实现 MMT 集中式 `intra_process_halo_exchange`
4. Worker 计算路径切换到 `ctx->gf->u_*`
5. 验证组间 halo 正确（对比阶段一的单组结果）

### 阶段三：NUMA 感知分配

1. 引入 `<numa.h>`（`#ifdef HAS_LIBNUMA`）
2. `group_alloc_field` 改为 `numa_alloc_onnode`
3. 增加 NUMA 节点信息输出
4. 在 NUMA 机器上验证本地/远端访存比例

### 阶段四：清理与优化

1. 移除旧版全局 MPI halo buffer 的冗余分配
2. `free_simulation` → `free_all_group_fields`
3. 非组模式路径适配
4. 更新文档

---

## 十、风险与缓解

| 风险 | 影响 | 缓解 |
|------|------|------|
| 组间 halo 数据竞争 | MMT 和 GMT 同时访问相邻组的 halo 行 | MMT 在 GMT 进入 gWaitMain 后独占所有组数据，此时 GMT 阻塞，无竞争 |
| NUMA API 不可用 | macOS / 无 libnuma 环境编译失败 | `#ifdef HAS_LIBNUMA` 条件编译 + configure 检测 |
| first-touch 失效 | calloc 不触发物理页分配 | 分配后用 `memset` 或显式写每页的首字节强制触发 page fault |
| 组间 halo 冗余内存 | 每组多 `2*HALO*nx*8*3` 字节 | HALO=1, nx≈1000 时仅 ~48KB/组，忽略不计 |
| 非组模式的代码路径一致性 | 两组代码路径不同导致行为分歧 | 非组模式也走 `GroupField[1]`，统一接口 |
| 调试复杂度 | 每组独立数组，gdb 查看不便 | 提供 `dump_group_field(gid)` 辅助函数 |

---

## 十一、API 变更摘要

| 函数 | 变更类型 | 说明 |
|------|----------|------|
| `idx()` | 签名变更 | `idx(y,x)` → `idx(gf,y,x)` |
| `global_*_from_local()` | 签名变更 | 增加 `const GroupField *gf` |
| `init_simulation()` | 删减 | 只保留拓扑初始化，移除波场分配 |
| `free_simulation()` | 替换 | 新增 `free_all_group_fields()` |
| `task_init_field()` | 上下文变更 | RowTaskCtx 增加 `gf` 指针 |
| `task_compute_interior()` | 同上 | 同上 |
| `task_compute_boundary()` | 同上 | 同上 |
| `task_compute_energy()` | 同上 | 同上 |
| `exchange_halos_for()` | 重写 | 拆为 `intra_process_halo_exchange` + `boundary_mpi_halo_exchange` |
| `group_main_thread()` | 扩展 | 增加 `group_alloc_field` / `group_free_field` |
| `main_thread()` | 扩展 | 增加组间 halo 交换、`g_gfields` 管理 |
| `thread_run()` | 微调 | 不变 |
| `_getgdsize_()` | 可能变更 | 若需要在 group data (`gi->gd`) 中存储组级数据 |
