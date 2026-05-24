# wave_propagation_ghost_omp.c → wave_propagation_ghost_fix.c 转换指南

## 速查：OMP 语法 → mythread 语法

| 你要做的事 | OMP 版写法 | mythread-fix 版写法 |
|-----------|-----------|---------------------|
| 并行执行 | `#pragma omp parallel` | `StartThreads(thread_run)` |
| 等待所有线程 | `#pragma omp barrier` | `GWM` / `GSM` / `GWS` / `SWG` 等宏 |
| 单线程执行 | `#pragma omp single` | 放在 `main_thread()` 中（天然单线程） |
| 循环并行化 | `#pragma omp for schedule(static)` | `setup_thread_tasks()` 预切 tile + Worker 在自己的 tile 上 `for` |
| 获取线程 ID | `int tid = omp_get_thread_num()` | `ti->ind` （TLS，框架提供） |
| 线程私有数据 | `g_times[tid * N]` | `ThreadTask *t = (ThreadTask*)ti->td` |
| 跨线程求和 | `#pragma omp for reduction(+:pe)` | Worker 写 `task->partial_energy`，GMT 汇总 `reduce_group_worker_energy()` |
| 跨线程互斥 | `#pragma omp critical` | **不需要**（Worker 写私有 tile，无冲突） |

---

## 总览

OMP 版 789 行，mythread 版 1301 行。多出来的 ~500 行分布在**线程角色分发**（~60 行）、**同步代码**（~100 行）、**手动 tile 分配**（~80 行）、**NUMA 感知分配**（~100 行）、**分组模式支持**（~100 行）、**计时/诊断**（~60 行）。

两个版本的计算核心（`compute_region`/`compute_block_region`、`compute_energy_block`/`compute_energy`、`init_field`/`init_field_block`）**完全相同**，不需要改任何一行。

---

## 两类差异一览

| | OMP 版（从） | mythread-fix 版（到） |
|---|---|---|
| 头文件 | `<omp.h>` + `mythread_config.h` + `mythread_decomp.h` | `mythread/mythread.h`（聚合全部） |
| 线程入口 | 无（`#pragma omp parallel` 内联） | `thread_run()` — 显式函数 |
| 同步 | `#pragma omp barrier` / `single` | `MSG/GWM/GSS/SWG` 等宏 |
| 线程局部数据 | `g_times[omp_get_thread_num() * N]` | `ti->td`（`_gettdsize_()` 指定大小） |
| 工作划分 | `#pragma omp for schedule(static)` | `setup_thread_tasks()` 手动切 Y 方向 |
| GroupField | 简化版（malloc） | 完整版（`group_alloc_field`，NUMA） |
| Halo 交换 | 内联 `halo_exchange_mpi()` | `mythread_halo_exchange_intra` + `_mpi` |
| 分组支持 | 无（永远单组） | `uses_group_threads()` 双路径 |
| 不均衡测试 | 无 | `apply_imbalance_to_group()` |
| 线程生命周期 | 无（`omp parallel` 自动管理） | `InitThreads` + `StartThreads` + `EndThreads` |

---

## 12 步逐项对照修改

### 第 1 步：头文件替换

```diff
- #include <omp.h>
- #include "mythread/mythread_config.h"   /* 被 mythread.h 间接包含 */
- #include "mythread/mythread_decomp.h"   /* 被 mythread.h 间接包含 */
+ #include "mythread/mythread.h"
```

`mythread.h` 已经聚合了所有子模块。

### 第 2 步：恢复完整的 cfg 变量

OMP 版删掉了一些 mythread 专用配置：

```diff
  static int cfg_HALO=1, ...
+ static int cfg_N_GROUPS=2;           /* 新增：线程组数 */
+ static int cfg_THREADS_PER_GROUP=4;  /* 新增：每组线程数（含 GMT） */
+ static int cfg_CoreOffset=0;         /* 新增：CPU 核心偏移 */
+ static int cfg_ClustOffset=0;        /* 新增：NUMA 簇偏移 */
```

`cfg_compute_derived()` 末尾补充：
```diff
+ cfg_THREADS_PER_GROUP = cfg_N_WORKERS + 1;
```

### 第 3 步：GroupField 从简化版切到完整版

OMP 版自己定义了简化版 `GroupField`（73-83 行）和 `GFIDX`（85 行）。全部删掉，改用 `mythread_field.h` 提供的版本。

```diff
- /* ── 简化版 GroupField（去 NUMA 字段，无 mythread_field 依赖）── */
- typedef struct GroupField {
-   double *u_prev, *u_curr, *u_next;
-   int ny_padded, nx_padded, stride;
-   size_t plane_bytes;
-   double *send_up, *recv_up, *send_down, *recv_down;
-   double *send_left, *recv_left, *send_right, *recv_right;
-   double group_energy;
- } GroupField;
- #define GFIDX(gf, y, x) ((size_t)(y) * (size_t)(gf)->stride + (size_t)(x))
```

完整版 `GroupField` 多了这些字段（新增功能，旧访问方式不变）：
- `numa_node, numa_ok` — NUMA 节点信息
- `send_to_up/down/left/right` — 组间发送缓冲（指向邻居 recv 缓冲的指针）
- `recv_from_up/down/left/right` — 组间接收缓冲

### 第 4 步：计时从 OMP 数组切到 ThreadTask

OMP 版的计时：
```c
/* OMP 版 */
double *g_times;   /* [max_threads * TM_N_SLOTS] */
int tid = omp_get_thread_num();
g_times[tid * TM_N_SLOTS + TM_COMPUTE] += wall_time() - t0;
```

mythread 版的计时：
```c
/* mythread 版 */
typedef struct {
    int gid, tid, cpu_id;
    int x_begin, x_end, y_begin, y_end;
    GroupField *gf;
    double t_wait_compute, t_work_compute;
    double t_wait_boundary, t_work_boundary;
    double t_wait_energy, t_work_energy;
    double t_comm, t_allreduce;
    int energy_steps;
} ThreadTask;

int _gettdsize_() { return (int)sizeof(ThreadTask); }  /* 框架自动为每个线程分配 */

/* 使用时 */
ThreadTask *task = (ThreadTask*)ti->td;            /* ti 是 TLS 变量，框架提供 */
task->t_work_compute += wall_time() - t0;
```

**关键点**：`_gettdsize_()` 告诉 mythread 每个线程需要多少字节私有数据，框架在创建线程时自动分配，通过 `ti->td` 访问。

### 第 5 步：同步机制替换（改动最集中的区域）

这是整个转换中改动量最大的地方。OMP 的 `#pragma omp barrier/single` 全部替换为 mythread 的三级屏障宏。

```
OMP 版                              mythread-fix 版
─────────────────────────────────────────────────────────
#pragma omp parallel               thread_run() 中的角色分发
#pragma omp barrier                GWM / GSM / GSS / GWS / SWG / SSG
#pragma omp single                 MMT 中的顺序执行（天然单线程）
#pragma omp for schedule(static)   setup_thread_tasks() 预分配 tile
#pragma omp critical               不需要！每个线程只写自己的 tile 区域
```

**具体的逐行替换**：

| OMP 版行 | 替换为 | 位置 |
|----------|--------|------|
| `#pragma omp single { apply_dirichlet_all(); ... }` | MMT 中直接调用（MMT 天然单线程） | `main_thread()` |
| `#pragma omp barrier`（compute 前） | `start_phase_from_main(cs)` → `wait_phase_from_main(cs)` | `main_thread()` |
| `#pragma omp for`（compute） | Worker 中 `wait_phase_from_worker(cs)` → `compute_*_block()` → `finish_phase_from_worker(cs)` | `worker_thread()` |
| `#pragma omp single { group_field_swap(gf); }` | MMT 中直接调用 | `main_thread()` |
| `#pragma omp for reduction(+:pe)`（energy） | Worker 中各自计算 `partial_energy`，GMT 中汇总 `reduce_group_worker_energy()` | `worker_thread()` / `group_main_thread()` |
| `#pragma omp single { MPI_Allreduce(...); }` | MMT 中直接调用 | `main_thread()` |
| `#pragma omp critical { gf->group_energy += pe; }` | 不需要！Worker 写入自己的 `partial_energy`，最后由 GMT/MMT 汇总 | — |

只需记住这 6 个宏：
```
GWM = gWaitMain      GSM = gSetMain      ← GMT ↔ MMT
GSS = gSetSubs       GWS = gWaitSubs     ← GMT ↔ Worker
SWG = sWaitGrp       SSG = sSetGrp       ← Worker ↔ GMT
```

### 第 6 步：从 OMP 并行区域到 thread_run 角色分发

删除整个 `run_simulation()` 中的 `#pragma omp parallel` 块，改为三个独立函数 + 一个分发器。

OMP 版的结构（单函数内嵌并行）：
```c
void run_simulation(int nthreads) {
    gf_alloc(...);
    #pragma omp parallel num_threads(nthreads) {
        int tid = omp_get_thread_num();
        // 所有线程在这里执行相同代码
        #pragma omp single { /* 串行部分 */ }
        #pragma omp for { /* 并行部分 */ }
    }
}
```

mythread 版的结构（三个函数 + 分发器）：
```c
void worker_thread(void) {     /* 每个 Worker 执行 */
    ThreadTask *task = (ThreadTask*)ti->td;
    /* 等待 GMT 信号 → 在自己的 tile 上计算 → 通知 GMT */
}

void group_main_thread(void) { /* 每个 GMT 执行 */
    ThreadTask *task = (ThreadTask*)ti->td;
    group_alloc_field(...);    /* NUMA 感知分配 */
    /* GWM → GSS → GWS → GSM 循环 */
}

void main_thread(void) {       /* MMT 执行（每个进程一个） */
    ThreadTask *task = (ThreadTask*)ti->td;
    /* Halo 交换 + 阶段同步 + MPI 通信 */
}

void thread_run(void) {        /* 框架回调，每个线程执行一次 */
    if (ti->igrp < 0)      main_thread();       /* MMT */
    else if (ti->ind == 0)  group_main_thread(); /* GMT */
    else                    worker_thread();      /* Worker */
}
```

### 第 7 步：main() 中增加线程生命周期管理

```diff
  int main(int argc, char **argv) {
      MPI_Init(&argc, &argv);
      /* ... 配置加载、MPI 拓扑、域分解 ... */

+     int ManageCoreId = cfg_ManageCoreId;
+     int err = InitThreads(mpi_rank, NCorePClu, NCluPNode, NCorePGrp,
+                           NThPGrp, NGrpPProc, NProcPNode, &ManageCoreId);
+     if (err != 0) { /* 错误处理 */ }

+     setup_thread_tasks();       /* ← 第 8 步新增的函数 */
+     StartThreads(thread_run);   /* 启动所有 Worker + GMT */
+     thread_run();                /* MMT 也进入 */
+     EndThreads();                /* 等待全部退出 */

-     run_simulation(nthreads);
      /* ... 清理 ... */
  }
```

### 第 8 步：增加手动 tile 分配

OMP 版用 `#pragma omp for schedule(static)` 自动分配迭代，mythread 版需要手动切分。

mythread 版的 `setup_thread_tasks()` 做的事情：
```
对每个线程组 g:
   对组内每个 worker t:
       该 worker 的 Y 范围 = HALO + t * (ny / nworkers) ~ HALO + (t+1) * (ny / nworkers)
       该 worker 的 X 范围 = HALO ~ HALO + nx       （覆盖全部 X）
       写入 ThreadTask.y_begin / y_end / x_begin / x_end
```

**注意**：如果在分组模式下，group 0 可能少一个 worker（MMT 占用该组一个槽位），代码中有 `ib++` 的处理逻辑。

### 第 9 步：GroupField 分配改为 group_alloc_field

```diff
- gf_alloc(gf, gid, g_decomp)       /* OMP 版：简化 malloc */
+ group_alloc_field(gf, gid, g_decomp)  /* mythread 版：NUMA 感知 */
```

**关键**：`group_alloc_field` 必须在 GMT 线程中调用（此时 `bindcpu` 已生效，`sched_getcpu()` 返回目标核心，`numa_node_of_cpu()` 可获得正确的 NUMA 节点）。

OMP 版的 `gf_alloc` 在 `main()` 中调用（单线程），mythread 版的 `group_alloc_field` 在 `group_main_thread()` 中调用（GMT 线程）。

### 第 10 步：Halo 交换从内联改为 mythread_halo

```diff
- halo_exchange_mpi(gf, gid);     /* OMP 版：内联 MPI Isend/Recv */
+ mythread_halo_exchange_intra(g_gfields, g_decomp);  /* 组间 memcpy */
+ mythread_halo_exchange_mpi(g_gfields, g_decomp, &g_mpi_ctx, (uintptr_t)MPI_COMM_WORLD);
```

新增了一步 `_intra`（组间 memcpy）。这是因为 mythread 版支持**多线程组**，组间边界需要本地拷贝。如果应用只用单组（`NGrpPProc=1`），`_intra` 自动为空操作。

### 第 11 步：增加分组/非分组的双路径支持

`uses_group_threads()` 几乎出现在每个汇总函数中：

```c
static double accumulate_worker_energy(void) {
    double e = 0.0;
    if (uses_group_threads()) {
        for (int g = 0; g < g_decomp->n_groups; g++)
            e += g_gfields[g].group_energy;       /* 分组：遍历 GroupField */
    } else {
        for (int t = 0; t < md.Nthreads - 1; t++)
            e += ((ThreadTask*)md.threads[t].td)->partial_energy;  /* 非分组：遍历 ThreadTask */
    }
    return e;
}
```

**简化技巧**：如果不需要分组模式，可以始终设 `NGrpPProc=1`，删除所有 `uses_group_threads()` 分支，只保留 `else` 路径。

### 第 12 步：增加计时报告函数

OMP 版用数组索引打印，mythread 版遍历 `md.threads[].td` 和 `md.grps[].threads[].td`：

```c
void print_timing_report(int mpi_rank, int mpi_size, int node_size) {
    // 遍历 MMT
    print_one((ThreadTask*)md.tm.td);
    // 遍历各组各线程
    if (uses_group_threads()) {
        for (int g = 0; g < g_decomp->n_groups; g++)
            for (int t = 0; t < md.grps[g]->Nthreads; t++)
                print_one((ThreadTask*)md.grps[g]->threads[t].td);
    } else {
        for (int t = 0; t < md.Nthreads - 1; t++)
            print_one((ThreadTask*)md.threads[t].td);
    }
}
```

---

## 转换检查清单

- [ ] `#include <omp.h>` → `#include "mythread/mythread.h"`
- [ ] 恢复 `cfg_N_GROUPS`, `cfg_THREADS_PER_GROUP`, `cfg_CoreOffset`, `cfg_ClustOffset`
- [ ] 删除自定义 `GroupField` 结构体，改用 `mythread_field.h` 版本
- [ ] 增加 `ThreadTask` 结构体和 `_gettdsize_()`
- [ ] 所有 `g_times[omp_get_thread_num() * N]` → `((ThreadTask*)ti->td)->t_xxx`
- [ ] `#pragma omp parallel` → `thread_run()` + `StartThreads/EndThreads`
- [ ] `#pragma omp barrier` → `GWM/GSM/GSS/GWS/SWG/SSG`
- [ ] `#pragma omp single` → 直接放在 MMT 函数中
- [ ] `#pragma omp for` → `setup_thread_tasks()` 预切 tile + Worker 在自己的 tile 上循环
- [ ] `#pragma omp critical` → 删除（Worker 写私有 tile、GMT/MMT 汇总）
- [ ] `gf_alloc` → `group_alloc_field`（并在 GMT 线程中调用）
- [ ] `halo_exchange_mpi` → `mythread_halo_exchange_intra` + `_mpi`
- [ ] main 中增加 `InitThreads` / `StartThreads` / `EndThreads`
- [ ] 增加 `main_thread()` / `group_main_thread()` / `worker_thread()` 三个函数
- [ ] 增加 `uses_group_threads()` 双路径（或设为单组模式简化）

---

## 注意事项

**1. `ti` 和 `gi` 只在 `StartThreads` 之后才有效**

在 `main()` 中，`StartThreads` 调用之前，`ti` 指向 MMT 的伪线程信息。`StartThreads` 之后创建的线程才有完整的 `ti` 和 `gi`。所以 `setup_thread_tasks()` 中访问 `ti->td` 需要在对应的线程入口函数中做，不能在 `main()` 中做。

**2. GroupField 的分配时机必须在 bindcpu 之后**

`group_alloc_field` → `numa_alloc_onnode(node)` → `numa_node_of_cpu(sched_getcpu())`。如果分配时线程还没绑核，NUMA 节点会错误。

在 mythread 中，绑核发生在 `threadMain()` → `_threadmain_()` → `bindcpu(ti->indg)`，而 `group_alloc_field` 在 GMT 的 `group_main_thread()` 中调用（此时已绑核）。所以 GMT 中分配 GroupField 是安全的。

**3. 不要跳过 `_gettdsize_()`**

如果忘了定义这个函数，`ti->td = NULL`，后续 `(ThreadTask*)ti->td` 解引用直接崩溃。框架的默认弱符号返回 `sizeof(float)`，不足以容纳 `ThreadTask`。

**4. 同步宏依赖 ThreadG 的值**

`GWM`/`SWG` 等宏内部根据 `ThreadG` 自动选择分组/非分组路径。`ThreadG` 由 `InitThreads` 根据 `NGrpPProc` 自动设置——`NGrpPProc > 1` 时 `ThreadG=1`，否则 `ThreadG=0`。

所以如果只想用非分组模式，设置 `NGrpPProc=1` 即可，同步宏自动走 `SM/MS` 路径。

**5. 不要并行调用 `StartThreads` / `EndThreads`**

这两个函数操作全局线程状态，必须在单线程中调用。MPI 多进程各自调用没问题。

**6. OMP 的 reduction 在 mythread 中怎么替代**

OMP 的 `reduction(+:pe)` 自动完成跨线程求和。mythread 中需要手动汇总：

```c
/* OMP 版：单行搞定 */
#pragma omp for reduction(+:pe)

/* mythread 版：两步 */
// Step 1: 每个 Worker 计算自己的部分，存入 task->partial_energy
// Step 2: GMT 调用 reduce_group_worker_energy(gid) 汇总
```

这是因为 mythread 的 Worker 之间没有隐式同步，需要显式通过 GMT 做归约（在 `gWaitSubs` 之后）。

**7. Halo 交换需要 `group_field_link_buffers` 先执行**

这个函数在初始化时调用一次，把相邻组的 `send_to_*` 指针互连到对方的 `recv_from_*` 缓冲。必须在 `halo_exchange_intra` 之前调用，否则 `send_to_*` 全为 NULL，不交换任何数据。

在代码中正确的位置是：所有组的 `group_alloc_field` 完成之后、主循环之前。
