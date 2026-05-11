
# mythread 接口函数文档

## 目录

1. [线程初始化和控制](#1-线程初始化和控制函数)
2. [线程信息获取](#2-线程信息获取函数)
3. [同步函数](#3-同步函数)
4. [性能计时 (TSC)](#4-性能计时函数-tsc)
5. [任务池](#5-任务池)
6. [域分解](#6-域分解-mythread_decomp)
7. [NUMA 感知字段分配](#7-numa-感知字段分配-mythread_field)
8. [Halo 交换](#8-halo-交换-mythread_halo)
9. [配置文件解析](#9-配置文件解析-mythread_config)
10. [Fortran 接口](#10-fortran-接口函数汇总)
11. [辅助函数](#11-其他辅助函数)
12. [全局变量](#12-全局变量)
13. [宏常量](#13-宏常量)
14. [使用示例](#14-使用示例)

---

## 1. 线程初始化和控制函数

**头文件**: `mythread_thread.h`

| 函数 | 参数 | 返回值 | 说明 |
|------|------|--------|------|
| `InitThreads` | `int mpi_id_, NCorePClu_, NCluPNode_, NCorePGrp_, NThPGrp_, NGrpPProc_, NProcPNode_, int *ManageCoreId_` | `int` (0=成功, <0=参数越界) | 初始化线程系统，内部调用 `initmd()` 建立全局 `threadProc md` 和所有线程组的 `THREADINFO`，计算每个线程的 `indg`（CPU 核心编号） |
| `StartThreads` | `TFunc tfun` | void | 为所有 worker 线程调用 `pthread_create`，入口函数为 `tfun` |
| `EndThreads` | void | void | `pthread_cancel` + `pthread_join` 终止所有 worker |
| `setthread` | `int NCorePClu_, NThPGrp_, NGrpPProc_, NProcPNode_, ManageCoreId_` | void | 在 `InitThreads` 之前单独设置配置参数 |
| `bindcpu` | `int id` | `int` (0=成功, -1=失败) | 将当前线程绑定到 CPU 核心 `id`（Linux `sched_setaffinity`；macOS 空操作） |
| `bindthread` | void | void | 调用 `bindcpu(ti->indg)` 绑定当前线程到预设核心 |
| `thread_run` | void | void | 用户定义的线程入口函数（`__attribute__((weak))`，可覆盖） |

**Fortran 接口**: `inithreads_`, `startthreads_`, `endthreads_`, `bindthread_`

---

## 2. 线程信息获取函数

| 函数 | 参数 | 返回值 | 说明 |
|------|------|--------|------|
| `gettid` | void | `int` | 返回当前线程在组内的索引 `ti->ind` |
| `getnt` | void | `int` | 返回当前线程所在组的线程总数 `ti->Nthreads` |
| `getcpuid` | void | `int` | 返回当前核心编号（Linux `sched_getaffinity`；macOS `sysctl hw.logicalcpu`） |

---

## 3. 同步函数

**头文件**: `mythread_sync.h`

三级屏障同步协议。命名规则：`s`=sub thread, `m`=main thread, `g`=group main thread。

### 3.1 Sub ↔ Main (SM / MS)

| 函数 | 参数 | 方向 |
|------|------|------|
| `sSetState(int state)` | 状态值 | Worker → MMT：设置 `ti->state` |
| `sWaitState(int state)` | 状态值 | Worker 等待 `md.state` 达到指定值 |
| `sWaitStater(int state)` | 状态值 | 同上（反向条件：≤） |
| `mSetSubs(int state)` | 状态值 | MMT → Worker：设置 `md.state` |
| `mWaitSubs(int state)` | 状态值 | MMT 等待所有 `md.threads[i].state` 达到指定值 |
| `mWaitSubsr(int state)` | 状态值 | 同上（反向条件） |

**宏**: `SSM = sSetState(RFB)` / `SWM = sWaitState(RFB)` / `MSS = mSetSubs(RFB)` / `MWS = mWaitSubs(RFB)`

### 3.2 Sub ↔ Group (SG / GS)

| 函数 | 方向 |
|------|------|
| `sSetGrp(int state)` | Worker → GMT：设置 `gi->sstate[ti->ind*MSBG]` |
| `sWaitGrp(int state)` | Worker 等待 `gi->state` |
| `gSetSubs(int state)` | GMT → Worker：设置 `gi->state` |
| `gWaitSubs(int state)` | GMT 等待所有 worker 的 `sstate` |

**宏**: `SSG = sSetGrp(RFB)` / `SWG = sWaitGrp(RFB)` / `GSS = gSetSubs(RFB)` / `GWS = gWaitSubs(RFB)`

### 3.3 Group ↔ Main (GM / MG)

| 函数 | 方向 |
|------|------|
| `gSetMain(int state)` | GMT → MMT：设置 `md.gstate[(ti->igrp)*MSB]` |
| `gWaitMain(int state)` | GMT 等待 `gi->gstate` |
| `mSetGrps(int state)` | MMT → GMT：设置 `md.state` 和所有 `md.grps[i]->gstate` |
| `mWaitGrps(int state)` | MMT 等待所有 `md.gstate[i*MSB]` |

**宏**: `GSM = gSetMain(RFB)` / `GWM = gWaitMain(RFB)` / `MSG = mSetGrps(RFB)` / `MWG = mWaitGrps(RFB)`

> 所有 Wait 函数有 `r` 后缀变体（如 `mWaitSubsr`），使用反向比较条件（`≤` 代替 `≥`）。默认 `RFB=1`。

---

## 4. 性能计时函数 (TSC)

**头文件**: `mythread_timer.h`

| 函数 | 说明 |
|------|------|
| `tscinit()` | 清零当前线程的计时器 |
| `tscb(int id)` | 记录时间戳到 `ti->tscbs[id]`（ARM `cntvct_el0` / x86 `rdtsc`） |
| `tsce(int id)` | 累加 `atsc() - tscbs[id]` 到 `ti->tscds[id]` |
| `tsceb(int id)` | `tsce(id-1)` + `tscb(id)` 合二为一 |
| `prtsc(const char *tag)` | 按秒输出 16 个计时槽的累计值 |
| Fortran: `tscb_`, `tsce_`, `tsceb_`, `prtsc_` | 指针参数解引用后调用 C 版本 |

---

## 5. 任务池

**头文件**: `mythread_pool.h`

SPMC（单生产者多消费者）固定容量环形队列，epoch 驱动。

| 函数 | 参数 | 返回值 | 说明 |
|------|------|--------|------|
| `mt_taskpool_attach` | `int slot, int capacity, int flags` | `int` (0=成功) | 创建任务池，挂载到 `gi->locv[slot]`（组模式）或 `md.locv[slot]`（非组模式） |
| `mt_taskpool_detach` | `int slot` | `int` | 解挂释放（须先 shutdown 且无在途任务） |
| `mt_taskpool_begin` | `int slot` | `int` | 开启新 epoch，重置队列，广播唤醒所有 worker |
| `mt_taskpool_submit` | `int slot, mt_task_fn fn, void *ctx` | `int` | 提交任务；队列满时的行为由 flags 决定 |
| `mt_taskpool_close` | `int slot` | `int` | 关闭本轮提交（worker 仅排空剩余任务） |
| `mt_taskpool_wait` | `int slot` | `int` | 阻塞至所有 worker 上报 `done_count >= nworkers` |
| `mt_taskpool_shutdown` | `int slot` | `int` | 通知所有 worker loop 退出 |
| `mt_taskpool_worker_loop` | `int slot` | `int` | worker 主循环：等 epoch → 取任务执行 → 报完成；收到 shutdown 后退出 |

**flags**:

| 宏 | 值 | 含义 |
|----|----|------|
| `MT_TASKPOOL_BLOCK` | 1 | 队列满时 `pthread_cond_wait` 阻塞 |
| `MT_TASKPOOL_SPIN` | 2 | 队列满时 `ntdelay(1)` 自旋 |
| `MT_TASKPOOL_TRY` | 4 | 队列满时立即返回 -4 |

**典型调用**（组模式）:
```
GMT:  attach → begin → submit*N → close → wait → (repeat) → shutdown → detach
Worker: worker_loop
```

**数据结构**:

```c
typedef void (*mt_task_fn)(void *ctx);
typedef struct { mt_task_fn fn; void *ctx; } mt_task;
typedef struct _mt_taskpool mt_taskpool;  /* 不透明 */
```

---

## 6. 域分解 (mythread_decomp)

**头文件**: `mythread_decomp.h`

提供二维域分解 + 邻居拓扑的通用接口，供 stencil/PDE 类应用使用。

### 6.1 枚举

```c
typedef enum {
  MYTHREAD_DECOMP_Y_ONLY = 0,   // 仅 Y 方向切分组，每组覆盖全部 X
  MYTHREAD_DECOMP_XY_2D  = 1    // X+Y 都切分，组排为 2D 网格
} mythread_decomp_policy;

typedef enum {
  MYTHREAD_NEIGHBOR_UP    = 0,
  MYTHREAD_NEIGHBOR_DOWN  = 1,
  MYTHREAD_NEIGHBOR_LEFT  = 2,
  MYTHREAD_NEIGHBOR_RIGHT = 3
} mythread_neighbor_dir;
```

### 6.2 数据结构

```c
typedef struct {
  int x_begin, x_end;           // 内部区域（不含 halo），左闭右开
  int y_begin, y_end;
  int nx, ny;                   // 内部区域尺寸
} mythread_tile;

typedef struct {
  int domain_nx, domain_ny;     // 全域尺寸
  int halo;                     // halo 宽度
  int n_groups;                 // 组数
  int n_workers_per_group;      // 每组 worker 数（0=不分配 worker tiles）
  mythread_decomp_policy policy;

  int gx, gy;                   // 组的 2D 网格排布
  mythread_tile *group_tiles;   // [n_groups] 每组子域
  int *group_neighbors;         // [n_groups × 4] 邻居拓扑（flat: gid*4+dir）
  mythread_tile *worker_tiles;  // [n_groups × n_workers] worker 子域
} mythread_decomp;
```

### 6.3 API

| 函数 | 说明 |
|------|------|
| `mythread_decomp *mythread_decomp_create(int nx, int ny, int halo, int ngroups, int nworkers, mythread_decomp_policy policy)` | 创建域分解，返回 NULL 表示参数无效 |
| `void mythread_decomp_free(mythread_decomp *dc)` | 释放 |
| `int mythread_decomp_neighbor(const mythread_decomp *dc, int gid, mythread_neighbor_dir dir)` | 查询邻居组 ID（-1=无） |
| `int mythread_decomp_is_domain_boundary(const mythread_decomp *dc, int gid, mythread_neighbor_dir dir)` | 该方向是否为域边界（需 MPI） |
| `const mythread_tile *mythread_decomp_worker_tile(const mythread_decomp *dc, int gid, int wid)` | 获取 worker tile |
| `void mythread_decomp_dump(const mythread_decomp *dc, FILE *out)` | 打印分解信息 |

---

## 7. NUMA 感知字段分配 (mythread_field)

**头文件**: `mythread_field.h`

为每个线程组提供独立的波场平面和 halo 缓冲，物理内存绑定在组所在 NUMA 节点。

### 7.1 数据结构

```c
typedef struct GroupField {
  double *u_prev, *u_curr, *u_next;  // 三个波场平面
  int ny_padded, nx_padded;          // 含 halo 的总行/列数
  int stride;                         // = nx_padded（索引宏用）
  size_t plane_bytes;                 // 单平面字节数

  double *send_to_up/down/left/right;     // 组间 halo 发送（指向邻居 recv）
  double *recv_from_up/down/left/right;   // 组间 halo 接收
  double *send_up/down/left/right;        // MPI halo（仅域边界组）
  double *recv_up/down/left/right;

  int numa_node;            // NUMA 节点 ID（-1 = 无 NUMA）
  double group_energy;      // 本组能量
} GroupField;
```

### 7.2 API

| 函数 | 说明 |
|------|------|
| `int group_alloc_field(GroupField *gf, int gid, const mythread_decomp *dc)` | 为组 gid 分配全部平面 + halo 缓冲。须在 GMT 线程中调用（bindcpu 已生效）。返回 0=成功 |
| `void group_free_field(GroupField *gf)` | 释放 |
| `void group_field_link_buffers(GroupField *gfields, const mythread_decomp *dc)` | 将相邻组的 `send_to_*` 指针互连到对方的 `recv_from_*` |
| `void group_field_fill(GroupField *gf, double val)` | 全部填充为 val |
| `void group_field_swap(GroupField *gf)` | 指针轮转：prev←curr, curr←next, next←prev |
| `void group_field_dump(const GroupField *gf, int gid, FILE *out)` | 打印调试信息 |

### 7.3 索引宏

```c
#define GFIDX(gf, y, x) ((size_t)(y) * (size_t)(gf)->stride + (size_t)(x))
```

等价于二维数组 `gf->u_curr[y][x]` 的线性索引。y 和 x 都以 **GroupField 内部坐标**（含 halo）为基准：`y=0` 为下方 halo，`y=HALO` 为第一个内部行，`y=HALO+ny` 为上方 halo。

### 7.4 NUMA 支持

- **Linux + libnuma**: `numa_alloc_onnode(bytes, numa_node_of_cpu(sched_getcpu()))`
- **macOS / 无 libnuma**: `malloc` + `memset` 强制 first-touch
- 自动检测：`__has_include(<numa.h>)`

---

## 8. Halo 交换 (mythread_halo)

**头文件**: `mythread_halo.h`

基于 `GroupField` 和 `mythread_decomp` 的邻居拓扑自动完成组间和 MPI halo 交换。

### 8.1 MPI 上下文

```c
typedef struct {
  int mpirank_up, mpirank_down;    // MPI 邻居 rank（-1=无）
  int mpirank_left, mpirank_right;
  int mpi_tag_base;                // 默认 100
} mythread_mpi_ctx;
```

### 8.2 API

| 函数 | 说明 |
|------|------|
| `void mythread_halo_exchange_intra(GroupField *gfields, const mythread_decomp *dc)` | 纯本地 memcpy：拷贝相邻组边界行/列到对方 recv 缓冲，然后应用到 halo 区域 |
| `void mythread_halo_exchange_mpi(GroupField *gfields, const mythread_decomp *dc, const mythread_mpi_ctx *ctx, uintptr_t comm)` | 仅域边界组执行 MPI Isend+Recv；`comm` 为 `(uintptr_t)MPI_COMM_WORLD` |
| `void mythread_halo_exchange_all(...)` | 等价于 `intra` + `mpi`（inline 函数） |

必须在 `group_field_link_buffers` 之后调用。

### 8.3 坐标约定

GroupField 中 `y=0` 映射到最小全局 Y（物理下方），`y=ny+HALO` 映射到最大全局 Y（物理上方）。up/down 方向以此为基准：
- `recv_from_up`（来自上方邻居）→ 写入 `row ny+HALO`（最大 Y）
- `recv_from_down`（来自下方邻居）→ 写入 `row 0`（最小 Y）

---

## 9. 配置文件解析 (mythread_config)

**头文件**: `mythread_config.h`

轻量级 INI 风格解析器，无外部依赖。

### 9.1 格式

```ini
# 注释
key = value
[section]
```

### 9.2 API

| 函数 | 说明 |
|------|------|
| `mythread_cfg *mythread_cfg_load(const char *path)` | 加载文件，失败返回 NULL |
| `void mythread_cfg_free(mythread_cfg *cfg)` | 释放 |
| `const char *mythread_cfg_get(cfg, section, key, defval)` | 查询字符串 |
| `int mythread_cfg_get_int(cfg, section, key, defval)` | 查询整数 |
| `double mythread_cfg_get_double(cfg, section, key, defval)` | 查询浮点 |
| `const char *mythread_env_get(name, defval)` | 读环境变量 |
| `int mythread_env_get_int(name, defval)` | 读整数环境变量 |

### 9.3 示例

```c
mythread_cfg *cs = mythread_cfg_load("config/case.cfg");
int nx = mythread_cfg_get_int(cs, "", "NX", 14000);
double dt = mythread_cfg_get_double(cs, "", "DT", 0.01);
mythread_cfg_free(cs);
```

---

## 10. Fortran 接口函数汇总

**头文件**: `mythread_fortran.h`

所有 Fortran 包装函数遵循 trailing-underscore 约定（`_` 后缀），参数为指针（Fortran call-by-reference → C `int *`）。

| C 函数 | Fortran 接口 | C 函数 | Fortran 接口 |
|--------|-------------|--------|-------------|
| `sWaitState` | `swaitstate_` | `sWaitStater` | `swaitstater_` |
| `sSetState` | `ssetstate_` | `mWaitSubs` | `mwaitsubs_` |
| `mWaitSubsr` | `mwaitsubsr_` | `mSetSubs` | `msetsubs_` |
| `gWaitSubs` | `gwaitsubs_` | `gWaitSubsr` | `gwaitsubsr_` |
| `gSetSubs` | `gsetsubs_` | `gWaitMain` | `gwaitmain_` |
| `gWaitMainr` | `gwaitmainr_` | `gSetMain` | `gsetmain_` |
| `mWaitGrps` | `mwaitgrps_` | `mWaitGrpsr` | `mwaitgrpsr_` |
| `mSetGrps` | `msetgrps_` | `sWaitGrp` | `swaitgrp_` |
| `sWaitGrpr` | `swaitgrpr_` | `sSetGrp` | `ssetgrp_` |
| `InitThreads` | `inithreads_` | `StartThreads` | `startthreads_` |
| `EndThreads` | `endthreads_` | `bindthread` | `bindthread_` |
| `ntdelay` | `ntdelay_` | `tscb/prtsc/...` | `tscb_/prtsc_/...` |

---

## 11. 其他辅助函数

**头文件**: `mythread_util.h` / `mythread_locv.h`

| 函数 | 参数 | 说明 |
|------|------|------|
| `SetLocV` | `int typ, int ind, void *p` | 设置局部指针槽位；`typ=0` 线程级（`ti->locv`），`typ=1` 组级（`gi->locv`），`typ=2` 进程级（`md.locv`）；`ind<8` 存直接指针，`ind≥8` 存堆扩展 |
| `GetLocV` | `int typ, int ind, void *p` | 获取局部指针；返回值 NULL 表示索引无效或未设置 |
| `ntdelay` | `int n` | `usleep(n/20)` |
| `opentf` | void | 打开线程调试日志文件（`THLOG` 宏控制） |

---

## 12. 全局变量

**头文件**: `mythread_types.h`（`extern` 声明）

| 变量 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `md` | `threadProc` | {0} | 进程级线程管理结构 |
| `ti` | `__thread THREADINFO*` | NULL | 当前线程信息（TLS） |
| `gi` | `__thread threadGroup*` | NULL | 当前线程组信息（TLS） |
| `ThreadG` | `int` | -1 | 分组标志：`InitThreads` 后 `NGrpPProc>1` 时=1 |
| `mpi_id` | `int` | 0 | MPI rank |
| `NCorePClu` | `int` | 38 | 每核心簇的核心数 |
| `NCluPNode` | `int` | 16 | 每节点核心簇数 |
| `NCorePGrp` | `int` | 38 | 每组核心数 |
| `NGrpPProc` | `int` | 4 | 每进程组数 |
| `NThPGrp` | `int` | 7 | 每组线程数 |
| `NProcPNode` | `int` | 16 | 每节点进程数 |
| `ManageCoreId` | `int` | -1 | 管理核心 ID（-1=自动，-2=末尾） |
| `NThreads` | `int` | 3 | 总线程数（含 MMT） |

---

## 13. 宏常量

**头文件**: `mythread_types.h`

| 宏 | 默认值 | 说明 |
|----|--------|------|
| `RFB` | 1 | 默认同步标志值 |
| `MCOREPC` | 64 | 每核心上下文最大线程数 |
| `MCLUST` | 20 | 最大簇/组数 |
| `MSB` | 16 | 进程级状态数组维度 |
| `MSBG` | 16 | 组级状态数组维度 |

**头文件**: `mythread_pool.h`

| 宏 | 值 | 说明 |
|----|----|------|
| `MT_TASKPOOL_BLOCK` | 1 | 队列满时阻塞 |
| `MT_TASKPOOL_SPIN` | 2 | 队列满时自旋 |
| `MT_TASKPOOL_TRY` | 4 | 队列满时立即返回 |

---

## 14. 使用示例

### 14.1 基础线程同步

```c
#include "mythread/mythread.h"

void thread_run(void) {
  if (ti->igrp < 0) {
    // MMT: 管理线程
    for (int step = 0; step < NT; step++) {
      mSetGrps(compute_phase_state(step));
      mWaitGrps(compute_phase_state(step));
    }
  } else if (ti->ind == 0) {
    // GMT: 组主线程
    gWaitMain(compute_phase_state(step));
    gSetSubs(compute_phase_state(step));
    gWaitSubs(compute_phase_state(step));
    gSetMain(compute_phase_state(step));
  } else {
    // Worker
    sWaitGrp(compute_phase_state(step));
    do_work();
    sSetGrp(compute_phase_state(step));
  }
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int mc = -1;
  InitThreads(0, 5, 2, 4, 4, 2, 1, &mc);
  StartThreads(thread_run);
  thread_run();  // MMT 入口
  EndThreads();
  MPI_Finalize();
}
```

### 14.2 域分解 + 每组独立分配 + 任务池

```c
// 创建域分解
g_decomp = mythread_decomp_create(local_nx, local_ny, HALO,
  n_groups, n_workers, MYTHREAD_DECOMP_Y_ONLY);

// 分配 GroupField 注册表
g_gfields = calloc(g_decomp->n_groups, sizeof(GroupField));

// GMT 线程中分配本组 GroupField（bindcpu 已生效）
group_alloc_field(&g_gfields[gid], gid, g_decomp);

// 组间 halo 缓冲互连
group_field_link_buffers(g_gfields, g_decomp);

// MMT 中：每步做 Dirichlet + halo 交换
apply_dirichlet_all_groups(g_gfields, g_decomp);
mythread_halo_exchange_intra(g_gfields, g_decomp);
mythread_halo_exchange_mpi(g_gfields, g_decomp, &mpi_ctx, (uintptr_t)MPI_COMM_WORLD);
apply_dirichlet_all_groups(g_gfields, g_decomp);

// Worker 访问: GFIDX(gf, y, x)
double v = gf->u_curr[GFIDX(gf, y, x)];
```

### 14.3 配置文件

```c
mythread_cfg *cs = mythread_cfg_load("config/case.cfg");
int nx = mythread_cfg_get_int(cs, "", "NX", 14000);
mythread_cfg_free(cs);

// 环境变量覆盖
const char *path = mythread_env_get("WAVE_CASE_CFG", "config/case.cfg");
```
