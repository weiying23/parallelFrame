# Wave Propagation Ghost 代码说明文档

## 1. 项目概述

`wave_propagation_ghost.c` 是一个基于 **MPI + pthread** 的二维波动方程并行求解器，使用幽灵层（Ghost Cells）进行边界交换。

### 主要特性
- **MPI 进程间通信**：Y 方向边界交换（跨节点）
- **组内线程间同步**：X 方向边界交换（共享内存）
- **分层线程架构**：主线程 → 组主线程 → 工作线程
- **幽灵层设计**：宽度为 1 的幽灵层用于有限差分计算

### 波动方程
```
∂²u/∂t² = c²(∂²u/∂x² + ∂²u/∂y²)
```

---

## 2. 架构设计

### 2.1 线程层次结构

```
MPI 进程 0                    MPI 进程 1
├─ 主线程 (tid=0)             ├─ 主线程 (tid=0)
│  └─ 管理所有组               │  └─ 管理所有组
├─ 组 0 主线程 (tid=1)         ├─ 组 0 主线程 (tid=1)
│  ├─ 工作线程 1 (tid=2)       │  ├─ 工作线程 1 (tid=2)
│  └─ 工作线程 2 (tid=3)       │  └─ 工作线程 2 (tid=3)
└─ 组 1 主线程 (tid=1)         └─ 组 1 主线程 (tid=1)
   ├─ 工作线程 1 (tid=2)          ├─ 工作线程 1 (tid=2)
   └─ 工作线程 2 (tid=3)          └─ 工作线程 2 (tid=3)
```

### 2.2 数据域分解

```
全局网格: NX × NY (默认 512 × 512)

Y 方向: 按 MPI 进程划分
        ┌─────────┬─────────┐
        │  Proc 0 │  Proc 1 │  ...
        │ Y:0-255 │ Y:256-511│
        └─────────┴─────────┘

Y 方向内部: 按组划分（每组负责一段 Y 范围）
        ┌─────────┬─────────┐
        │  Group 0│  Group 1│
        │ Y:0-127 │ Y:128-255│
        └─────────┴─────────┘

X 方向: 按组内线程划分
        ┌────┬────┬────┬────┐
        │T0  │T1  │T2  │T3  │
        │X:0-│X:128│X:256│X:384
        │127 │-255 │-383 │-511 │
        └────┴────┴────┴────┘
```

---

## 3. 核心数据结构

### 3.1 WaveField（波场数据）

```c
typedef struct {
    double *u_curr;     // 当前时刻波场（含幽灵层）
    double *u_prev;     // 前一时刻波场（含幽灵层）
    double *u_next;     // 下一时刻波场（含幽灵层）
    
    // 全局坐标范围
    int x_start, x_end; // X 方向负责范围
    int y_start, y_end; // Y 方向负责范围
    
    // 局部网格尺寸
    int nx, ny;             // 含幽灵层
    int nx_inner, ny_inner; // 不含幽灵层
    
    double energy;      // 局部能量
    
    // 邻居信息
    int neighbor_up;    // 上方 MPI rank
    int neighbor_down;  // 下方 MPI rank
    int igrp;           // 所属组 ID
    int itid;           // 组内线程 ID
} WaveField;
```

### 3.2 GroupBoundaryData（组内边界交换数据）

```c
typedef struct {
    double *left_send;   // 发送给左邻居的数据
    double *left_recv;   // 从左邻居接收的数据
    double *right_send;  // 发送给右邻居的数据
    double *right_recv;  // 从右邻居接收的数据
    int buf_size;        // 缓冲区大小（ny_inner）
    
    volatile int *thread_done;    // 各线程完成标记
    volatile int *boundary_ready; // 边界准备好标记
} GroupBoundaryData;
```

### 3.3 内存布局（含幽灵层）

```
局部网格索引 (i: X方向, j: Y方向):

          j=0          j=GHOST        j=ny_inner    j=ny-1
            │              │               │            │
    i=0     ├──────────────┼───────────────┼────────────┤
            │   幽灵层     │               │   幽灵层   │
    i=GHOST ├──────────────┼───────────────┼────────────┤
            │              │   内部计算区域 │            │
            │              │  nx_inner ×   │            │
            │              │  ny_inner     │            │
    i=nx_inner+GHOST ├─────┼───────────────┼────────────┤
            │   幽灵层     │               │   幽灵层   │
    i=nx-1  ├──────────────┼───────────────┼────────────┤

线性索引: idx(i, j, ny) = i * ny + j
```

---

## 4. 同步机制详解

### 4.1 同步宏定义

```c
// 主线程 ↔ 组主线程（进程级）
#define MSG  mSetGrps(RFB)    // 主线程设置所有组的 gstate
#define MWG  mWaitGrps(RFB)   // 主线程等待所有组完成

// 组主线程 ↔ 主线程（进程级）
#define GWM  gWaitMain(RFB)   // 组主等待主线程信号
#define GSM  gSetMain(RFB)    // 组主通知主线程完成

// 组主线程 ↔ 工作线程（组内）
#define GSS  gSetSubs(RFB)    // 组主通知工作线程开始
#define GWS  gWaitSubs(RFB)   // 组主等待工作线程完成

// 工作线程 ↔ 组主线程（组内）
#define SWG  sWaitGrp(RFB)    // 工作线程等待组主信号
#define SSG  sSetGrp(RFB)     // 工作线程通知组主完成
```

### 4.2 同步变量映射

| 层级 | 发送方 | 操作 | 等待方 | 操作 | 变量位置 |
|------|--------|------|--------|------|----------|
| 进程级 | 主线程 MSG | `gi->gstate = 1` | 组主 GWM | 等待 `gi->gstate != 1` | `threadGroup->gstate` |
| 进程级 | 组主 GSM | `md.gstate[grp] = 1` | 主线程 MWG | 等待 `md.gstate[grp] != 1` | `threadProc->gstate[]` |
| 组内 | 组主 GSS | `gi->state = 1` | 工作线程 SWG | 等待 `gi->state != 1` | `threadGroup->state` |
| 组内 | 工作线程 SSG | `gi->sstate[tid] = 1` | 组主 GWS | 等待 `gi->sstate[tid] != 1` | `threadGroup->sstate[]` |

### 4.3 SSG 级联等待机制

```c
void sSetGrp(int state) {
    // 如果 ti->sib > 0，表示这是"代表线程"，需要等待其他线程
    if(ti->sib) {
        for(int i = ti->sib; i < ti->sie; i++) {
            Wait_LGE(PSSTATE(i), state);  // 等待后续线程
        }
    }
    SSSTATE(ti->ind) = state;  // 设置自己的状态
}
```

**NGG = 8** 时的线程行为：
- 线程 ID 为 8 的倍数：`sib > 0`，需要等待后续 7 个线程
- 其他线程：`sib = 0`，直接设置状态

---

## 5. 执行流程

### 5.1 时间步循环流程

```
Step N:
┌─────────────────────────────────────────────────────────────────┐
│ 主线程                                                           │
│  1. MPI Y边界交换 (mpi_exchange_y_boundaries)                    │
│  2. 组间边界交换 (exchange_boundaries_between_groups)            │
│  3. MSG (设置 gi->gstate = RFB, 通知所有组开始)                   │
│  4. MWG (等待所有组的 md.gstate[] != RFB)                         │
│  5. 能量归约 (MPI_Allreduce)                                      │
└─────────────────────────────────────────────────────────────────┘
                              │
        ┌─────────────────────┼─────────────────────┐
        │                     │                     │
        ▼                     ▼                     ▼
┌───────────────┐    ┌───────────────┐    ┌───────────────┐
│   组 0 主线程  │    │   组 1 主线程  │    │   组 2 主线程  │
│  1. GWM (等待  │    │  1. GWM (等待  │    │  1. GWM (等待  │
│     gi->gstate │    │     gi->gstate │    │     gi->gstate │
│     != RFB)    │    │     != RFB)    │    │     != RFB)    │
│  2. GSS (设置  │    │  2. GSS (设置  │    │  2. GSS (设置  │
│     gi->state) │    │     gi->state) │    │     gi->state) │
│  3. 计算波场   │    │  3. 计算波场   │    │  3. 计算波场   │
│  4. GWS (等待  │    │  4. GWS (等待  │    │  4. GWS (等待  │
│     sstate[])  │    │     sstate[])  │    │     sstate[])  │
│  5. X边界交换  │    │  5. X边界交换  │    │  5. X边界交换  │
│  6. GSM (设置  │    │  6. GSM (设置  │    │  6. GSM (设置  │
│     md.gstate) │    │     md.gstate) │    │     md.gstate) │
└───────────────┘    └───────────────┘    └───────────────┘
        │                     │                     │
        ├──────────┬──────────┤                     │
        │          │          │                     │
        ▼          ▼          ▼                     │
┌──────────┐ ┌──────────┐ ┌──────────┐             │
│ 工作线程0 │ │ 工作线程1 │ │ 工作线程2 │             │
│ 1. SWG   │ │ 1. SWG   │ │ 1. SWG   │             │
│ 2. 计算  │ │ 2. 计算  │ │ 2. 计算  │             │
│ 3. SSG   │ │ 3. SSG   │ │ 3. SSG   │             │
└──────────┘ └──────────┘ └──────────┘             │
```

---

## 6. 已知问题

### 6.1 同步状态不重置问题（严重）

**问题描述**：
- 初始状态：`gi->gstate = 0`, `md.gstate[] = 0`, `RFB = 1`
- 第一轮：GWM 和 MWG 看到 `0 != 1`，立即通过，**无任何同步**
- 第一轮后：`gi->gstate = 1`, `md.gstate[] = 1`
- 第二轮：GWM 和 MWG 看到 `1 == 1`，**死锁等待**

**根本原因**：
```c
Wait_LGE(pv, cv)  // 等待 *pv != cv
```
但第一轮前状态已经是 `!= cv`（0 != 1），所以直接通过。

**修复建议**：
每轮开始前重置状态，或使用递增状态值：
```c
// 方案1：递增状态
int state = step + 1;
MSG(state);  // 设置 gi->gstate = step+1
MWG(state);  // 等待 md.gstate[] == step+1

// 方案2：重置状态
for (int g = 0; g < md.ngrp; g++) {
    md.gstate[g * MSB] = 0;  // 重置
}
MSG(RFB);
```

### 6.2 工作线程缺少边界交换等待

**问题描述**：
工作线程在 SSG 后没有等待边界交换完成：
```c
worker_thread():
    SWG;
    update_wavefield(wf);  // 计算
    SSG;                    // 通知完成
    // 缺少：等待边界交换完成
    sleep(2);              // 用 sleep 代替同步（不可靠）
```

**后果**：
- 工作线程可能在组主完成 X 边界交换前就开始下一轮计算
- 读取未更新的幽灵层数据，导致计算错误

**修复建议**：
```c
worker_thread():
    SWG;
    update_wavefield(wf);
    SSG;
    
    // 等待边界交换完成
    while (!gbd->boundary_ready[tid]) {}
    gbd->boundary_ready[tid] = 0;
    
    swap_buffers(wf);  // 工作线程也需要交换缓冲区
```

### 6.3 工作线程缺少缓冲区交换

**问题描述**：
```c
group_main_thread():
    // ...
    exchange_x_boundaries_group(gid);
    swap_buffers(wf);  // 只有组主交换

worker_thread():
    // 没有 swap_buffers！
```

**后果**：工作线程始终使用同一组缓冲区，`u_curr/u_prev/u_next` 指针不轮转。

### 6.4 SSG 级联等待的死锁风险

当 `N_THREADS >= 8` 时，线程 8 会等待线程 9-15，但线程 9-15 的 `sib = 0` 直接设置状态，不会等待线程 8。这可能导致级联等待链不完整。

---

## 7. 编译和运行

### 7.1 编译

```bash
mpicc -Wall -O2 -o wave_propagation_ghost \
    wave_propagation_ghost.c \
    mythread/mythread.c \
    -lpthread -lm
```

### 7.2 运行

```bash
mpirun -np 4 ./wave_propagation_ghost
```

### 7.3 参数配置

在代码顶部修改：
```c
#define NX          512     // X 方向网格数
#define NY          512     // Y 方向网格数
#define NT          1000    // 时间步数
#define N_GROUPS    2       // 每进程线程组数
#define N_THREADS   3       // 每组工作线程数（不含组主）
#define GHOST_WIDTH 1       // 幽灵层宽度
```

---

## 8. 性能指标

| 指标 | 说明 |
|------|------|
| Mgrid/s | 百万网格点/秒 = NX × NY × NT / time / 1e6 |
| 能量守恒 | 全局能量相对变化应接近机器精度 |
| 负载均衡 | 各线程计算量应大致相等 |

---

## 9. 文件结构

```
.
├── wave_propagation_ghost.c    # 主程序（含幽灵层版本）
├── wave_propagation.c          # 基础版本（无幽灵层）
├── wave_visual.c               # 可视化程序
├── visualize_wave.py           # Python 可视化脚本
├── mythread/
│   ├── mythread.h              # 线程库头文件
│   ├── mythread.c              # 线程库实现
│   └── mythread接口说明.md     # 接口说明（中文）
├── WAVE_EXAMPLE.md             # 波动方程示例说明
└── WAVE_PROPAGATION_GHOST_DOC.md  # 本文档
```

---

## 10. 参考资料

1. `mythread.h` - 同步原语定义
2. `mythread.c` - Wait_LGE, sSetGrp, gWaitSubs 等实现
3. `WAVE_EXAMPLE.md` - 波动方程数学原理
4. 有限差分法（FDM）教程
