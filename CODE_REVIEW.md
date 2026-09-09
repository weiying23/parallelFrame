# parallelFrame 代码审查报告

**审查日期**: 2026-09-09
**分支**: simplify
**审查方式**: 5 路并行审查（API 一致性 / QA 实测 / 代码质量 / 安全审计 / 历史挖掘），其中 4 项经 ASan/TSan/实际运行动态验证
**测试环境**: macOS (ARM) + 鲲鹏 608 核 aarch64 (f1 节点, HPCKit / bisheng + hMPI)

---

## 总体结论

**FAILED** — 5 路审查全部未通过，存在 P0 级已实证的内存损坏和数据竞争。

| # | 审查维度 | 结果 | 置信度 |
|---|---------|------|--------|
| 1 | 目标与 API 一致性 | FAIL | HIGH |
| 2 | QA 实测（构建+运行全部测试） | FAIL | HIGH |
| 3 | 代码质量（mythread 库深审） | FAIL | HIGH |
| 4 | 安全/内存安全（ASan/TSan 验证） | FAIL | CRITICAL |
| 5 | 历史与上下文挖掘 | FAIL | HIGH |

---

## 测试结果摘要

### 编译

| 项目 | 结果 |
|------|------|
| `make all`（11 个目标：2 主程序 + 9 单测） | ✅ 0 错误，0 警告 |
| `make debug`（-DDEBUG） | ❌ wave_propagation_ghost.c:652-655 编译失败 |

### 单元测试（f1 鲲鹏节点 + macOS 双平台一致）

| 测试套件 | 用例数 | 通过 | 失败 |
|----------|--------|------|------|
| bindcpu | 8 | 8 | 0 |
| locv | 9 | 9 | 0 |
| sync | 12 | 12 | 0 |
| field | 10 | 10 | 0 |
| config | 11 | 11 | 0 |
| decomp | 11 | 11 | 0 |
| timer | 6 | 6 | 0 |
| util | 4 | 4 | 0 |
| **taskpool** | **8** | **5** | **3** ❌ |
| **合计** | **79** | **76** | **3** |

### 端到端运行（f1 鲲鹏节点, 4000×4000×200 benchmark）

| 程序 | 进程数 | 结果 | 性能 |
|------|--------|------|------|
| `wave_propagation_ghost` | 1 | ✅ 9.518s | 336.22 Mpoint-updates/s |
| `wave_propagation_ghost_fix` | 1 | ✅ 8.107s | 394.73 Mpoint-updates/s |
| `wave_propagation_ghost` | 2 MPI | ✅ 4.807s | 665.73 Mpoint-updates/s |
| `wave_propagation_ghost_fix` | 2 MPI | ✅ 4.143s | 772.42 Mpoint-updates/s |
| `wave_propagation_ghost_fix` | 4 MPI | ❌ 死锁 | — |

---

## 阻塞性问题（按优先级）

### P0 — 已实证的内存损坏 / 数据竞争

#### 1. `_gettdsize_` 被注释掉 → td 仅 4 字节 → 堆溢出

**文件**: `wave_propagation.c:58-64`、`examples/example_mpi_mythread.c:37-44`

弱符号默认返回 `sizeof(float)=4`，随后向 td 写入 ~88 字节的 WaveField（ASan 实证 WRITE）。example_mpi_mythread 还因 `thread_run` 角色误判触发 `md.gstate[-16]` 结构体越界写。

**修复**: 恢复 `_gettdsize_` 定义；example_mpi_mythread 修正角色路由。

---

#### 2. 同步层全部基于裸 `volatile int` 自旋，无 acquire/release 内存序

**文件**: `mythread/mythread_sync.c:11-37`、`mythread_types.h:41-63`

TSan 实测 3 处数据竞争（含波场数据本体）。x86-TSO 掩盖问题，ARM/POWER 等弱序平台上是指针发布可读到陈旧值 → 野指针。

**修复**: 同步状态改为 C11 `_Atomic int` + `memory_order_acquire/release`；或 set/wait 两侧加 `atomic_thread_fence`。

---

#### 3. HALO 配置无校验 → 堆越界写 / 崩溃

**文件**: `wave_propagation_ghost.c:291`、`mythread_decomp_create` 返回值未判空

`HALO=0` → heap-buffer-overflow WRITE（ASan 实证）；`HALO<0` → 空指针 SIGSEGV。同一模式存在于全部 7 个 wave 变体（ghost.c、ghost_fix.c、ghost_mpi.c、ghost_omp.c、mpi_sync.c、mpi_mythread.c、mpi_omp.c）。

**修复**: cfg_load 后统一校验 `HALO >= 1`；`decomp_create` 返回值判空。

---

#### 4. locv 扩展数组 realloc 与并发 GetLocV 存在 use-after-free

**文件**: `mythread/mythread_locv.c:39-53` + `mythread_pool.c:49,78`

多池/中途 attach 场景下 worker 可解引用已释放指针。`locv[7]` 同时作为普通槽和扩展数组指针，双重语义产生冲突。

**修复**: 在 attach/detach 路径加一次性初始化协议，或改为定长直接寻址（`void*locv[104]`），或加锁。

---

### P1 — 正确性 / 可用性

#### 5. taskpool 测试 3/8 确定性失败，且被 `exit(0)` 掩盖

**文件**: `mythread/mythread_thread.c:189-219`、`mythread_sync.c:30`

三个子用例在两个平台（macOS + aarch64 鲲鹏）上完全一致地失败：
- `grouped-block-full`: "wait state timeout 0 < 1" — BLOCK 模式第二次 submit 永久阻塞
- `grouped-invalid-order`: "wait state timeout 0 < 1" — threads_per_group=1 时 group 0 线程数为 0
- `grouped-try-full`: 断言失败 — detach 返回 -2，任务永不执行

**根因**: `ManageCoreId≤0` 时 group 0 被内嵌 MMT（主管理线程）吞掉一个 worker（`pgi->Nthreads--`）。`mythread_sync.c:30` 超时后 `exit(0)` 把死锁报告为成功。

**修复**: InitThreads 拒绝 `NThPGrp<2`；`exit(0)` 改为 `EXIT_FAILURE` 或 `MPI_Abort`；文档化 ManageCoreId 语义。

---

#### 6. 非分组模式下组同步宏解引用 NULL `gi`

**文件**: `mythread_sync.c:113,148,163-176,196,206,223`

`NGrpPProc==1` 时 `ThreadG=0`，`gi==NULL`。但 `sWaitGrp/sSetGrp/gWaitSubs/gSetSubs/gWaitMain/gSetMain` 无任何 NULL 检查，直接 `gi->state`。文档承诺"宏按 ThreadG 自动选路径"与实现不符。

**修复**: 各 gi 族函数入口加 `if(!gi) return;` 或映射到 md 等价物；修正文档。

---

#### 7. halo 库两处静默算错

**文件**: `mythread/mythread_halo.c:27-37,108,118,133,143`

- **halo=1 硬编码**: 所有索引使用 `ny` 和 `0`，仅当 halo==1 时恰好正确。`mythread_decomp_create` 接受任意 halo 但 halo 交换只拷 1 行/列。
- **XY_2D Y 边界 break**: Y 上/下边界循环在第一个匹配组后 `break`（`mythread_halo.c:120,145`），X 方向循环无 break（`mythread_halo.c:150-199`），确认是疏漏。

**修复**: 索引改为 `halo+ny-1`/`halo-1` 等通式，循环拷 `halo` 行；去掉 Y 方向 `break`。

---

#### 8. `make debug` 编译失败

**文件**: `wave_propagation_ghost.c:652-655`

`#ifdef DEBUG` 块引用 ThreadTask 不存在的 5 个成员（`t_wait_energy0/t_wait_energy/t_work_energy0/t_work_energy/t_allreduce`），debug-only 代码未随结构体更新。

**修复**: 删除或改写该打印块。

---

#### 9. Makefile 损坏目标

**文件**: `Makefile:22-29,76-77`

- `run_wave` → `wave_propagation` 无规则（隐式规则链接失败）
- `example_simple/example_mpi_mythread/example_taskpool` 规则引用根目录不存在的 .c（实际在 `examples/`）
- `run_simple/run_full` 连带损坏

**修复**: 修正路径或删除规则，同步更新 WAVE_EXAMPLE.md。

---

#### 10. RFB 默认值为 1 → 便捷宏第 2 步起零同步

**文件**: `mythread_types.h:16-18`、`mythread_sync.h:34-56`

所有 `MSG/MWG/GWM/...` 宏展开为 `xxx(RFB)`，默认 `RFB=1`，而 Wait 全部是 `>=` 单调比较。第 1 步后所有状态位已是 1，从 step 2 开始所有 wait 立即通过，零同步。

**修复**: 默认宏不可用就不该提供；改为基于 `md.nstep` 的表达式，或删掉默认 RFB 强制用户定义。

---

### P2 — 健壮性 / 仓库卫生

#### 11. NUMA 回退路径中 `free()` 错误

**文件**: `mythread/mythread_field.c:31-62,131,143-159`

`fell_back` 是整组级标志：只要任一缓冲回退 malloc，所有缓冲（含成功 `numa_alloc_onnode` 的平面）都用 `free()` 释放 → 堆损坏。halo 缓冲释放时 `numa_free(p,0)` 违反 libnuma 契约。此路径仅 Linux 编译，从未测试过。

**修复**: 每个缓冲单独记录 `numa_ok`；释放时传真实尺寸。

---

#### 12. timer API 无边界检查

**文件**: `mythread/mythread_timer.c:20-51`

`tscb(id)/tsce(id)` 直接写 `tscbs[id]/tscds[id]`（固定 [16]），`id>=16` 越界写；`tsceb(0)` 写 `tscds[-1]` 下溢。Fortran 包装同样问题。

**修复**: 入口处 `if(id<0||id>=16) return;`。

---

#### 13. Fortran 接口缺失实现

**文件**: `mythread/mythread_fortran.h:23-25`

`swaitgrp_`/`swaitgrpr_`/`ssetgrp_` 声明了但无实现，Fortran 用户调用即链接错误。`ntdelay_(int n)` 按值传参，而同文件其他 Fortran 包装全部按 `int*` 引用传参。

**修复**: 补齐 SG 节实现，或从声明中删除；`ntdelay_` 改为 `int*`。

---

#### 14. pthread_create 返回值不检查 + EndThreads 死锁路径

**文件**: `mythread/mythread_thread.c:357,369,312-321`

`pthread_create` 不检查返回值。若第 k 个线程创建失败，主线程后续 `mWaitSubs/mWaitGrps` 永远等一个不存在的线程。`EndThreads`→`pthread_cancel`+`pthread_join`，若线程卡在自旋中（无取消点），`pthread_join` 永久挂起。

**修复**: 检查返回值，失败时回滚；改协作式退出（先广播退出态、join，超时未退出的再 cancel）。

---

#### 15. wave_visual.c 崩溃路径

**文件**: `wave_visual.c:176,179,186,244-261,281`

非分组模式下 setup 循环不执行 → td 带垃圾指针运行；worker 的 gi 为 NULL，`SWG` 解引用 NULL gi → SEGV。

**修复**: 非分组模式用 `md.threads` 遍历初始化 td；td 分配后 memset 0；检查 InitThreads 返回值。

---

#### 16. 仓库卫生

| 问题 | 详情 |
|------|------|
| `.gitignore` 不忽略二进制 | 清理后 28 分钟已重新生成 10 个未跟踪二进制 |
| 3 个二进制仍被跟踪 | `tests/test_sync`（已修改）、`mythread/test_bindcpu`、`mythread/test_field` |
| 重复源文件 | `mythread/test_bindcpu.c` 和 `mythread/test_field.c` 是与 `tests/` 版本发散的陈旧副本 |
| omp 变体未删净 | `wave_propagation_ghost_{omp,mpi,mpi_omp,mpi_sync,mpi_mythread}.c` 共 5 个仍被跟踪，与"simplify"目标相悖 |
| 死配置 | `ENERGY_REPORT_INTERVAL` 被读入却从未使用 |
| 死字段 | `THREADINFO.sync_flag/MainThread/para`、`threadProc.NSpecial/mpi_npe`、`threadGroup.mainGroup` 等只写不读 |

**修复**: 更新 `.gitignore`（忽略二进制）；删除被跟踪的二进制 + 陈旧副本；清理死代码/死配置。

---

#### 17. 文档过时

**文件**: `WAVE_PROPAGATION_GHOST_DOC.md`、`WAVE_EXAMPLE.md`、`scripts/benchmark_full.sh`

| 问题 | 详情 |
|------|------|
| 能量诊断 | 文档大篇幅描述已删除的 energy/l2/max|u| 输出，代码中不存在 |
| 状态编码 | 文档 "3s+3/3s+4/3s+5"，实际 "2s+2/2s+3" |
| RowTaskCtx | 文档有 8 个字段，实际只有 5 个 |
| 编译示例 | WAVE_EXAMPLE.md 指导编译 `wave_propagation` 目标，但 Makefile 无此规则 |
| benchmark 脚本 | `benchmark_full.sh:91` 从输出 grep `E=`/`L2=`，但两个 ghost 二进制已无此输出 → 列恒为 0 |

**修复**: 按 simplify 后代码重写文档；更新 benchmark 脚本。

---

## 非阻塞但值得注意的发现

### 并发正确性

| 问题 | 文件 | 说明 |
|------|------|------|
| StartThreads 复用不清零 | `mythread_thread.c:346` | 只复位 `md.state`，`gi->state`/`sstate` 保留旧值，第二轮 `>=` 单调比较可能被满足 |
| taskpool begin 无防御 | `mythread_pool.c:107-128` | 不查上一轮 `done_count==nworkers`，迟到的 worker 跨 epoch 重复计数 |
| 主线程与 worker tid 重叠 | `mythread_thread.c:272,344` | 主线程 `ti->ind=0`，worker `ind=j` 也从 0 开始 |
| MTM 的 gi 指向最后一组 | `mythread_thread.c:217-220` | `mt==-1` 时每个组都赋值，主线程 gi 最终指向最后一组，留着被误用 |

### 可移植性

| 问题 | 文件 | 说明 |
|------|------|------|
| `asm` vs `__asm__` | `mythread_timer.c:7,11` | `-std=c11` 下编译失败（`mythread_sync.c:21-23` 用的是 `__asm__`，同库两套写法） |
| getcpuid 语义不同 | `mythread_thread.c:47-70` | Linux 返回核号，macOS 返回核心数（如 12），同名函数两平台语义完全不同 |
| 函数指针转换 | `mythread_thread.c:357,369` | `(TSFunc)threadMain` 类型不匹配，clang 已告警 |
| `#include <mpi.h>` 全局 | `mythread_halo.c` | 无条件包含，halo 模块脱离 MPI 无法编译 |

### 性能

| 问题 | 文件 | 说明 |
|------|------|------|
| `MT_TASKPOOL_SPIN` 退避为 `usleep(0)` | `mythread_util.c:4-6`、`mythread_pool.c:147` | `ntdelay(1)` → `usleep(1/20=0)` 立即返回，纯热自旋加锁解锁 |
| MAXCHECK=1e9 自旋 | `mythread_sync.c:15` | 超时实际可能耗时极长（数分钟） |

### 测试

| 问题 | 文件 | 说明 |
|------|------|------|
| test_bindcpu 泄漏 | `mythread/test_bindcpu.c:186-246` | 两次 `memset(&md,0,...)+InitThreads` 丢弃已分配的 grps 指针；恢复时漏存全局变量 |
| `sprintf` 步进错误 | `mythread_timer.c:53-60` | `sprintf(buf+i*10,"|%9.5f     Z")` 每段写 17 字节但步进 10，输出互相覆盖 |
| `strtol` 八进制解析 | `mythread_config.c:53,109` | `strtol(...,0)` 使 `"010"` 按八进制解析 |

---

## 修复建议（按优先级）

### 第一轮 — 修复已实证的内存损坏（P0）

1. 恢复 `wave_propagation.c` 和 `example_mpi_mythread.c` 中注释掉的 `_gettdsize_` 定义
2. 修正 `example_mpi_mythread.c` 的 `thread_run` 角色路由
3. `cfg_load` 后统一校验 `HALO >= 1`，`decomp_create` 返回值判空
4. `mythread_sync.c:30` 的 `exit(0)` → `exit(EXIT_FAILURE)` 或 `MPI_Abort`
5. 同步状态字段改为 `_Atomic int` + `acquire/release` 语义

### 第二轮 — 修复正确性/可用性（P1）

6. InitThreads 文档化/修复 ManageCoreId 语义；拒绝 `NThPGrp<2`
7. 非分组模式给 gi 族函数加 NULL 守卫
8. `mythread_halo.c` 去掉 Y 方向两处 `break`，索引改 halo 通式
9. 删除 `wave_propagation_ghost.c:652-655` 的 `#ifdef DEBUG` 死块
10. 修 Makefile 损坏目标（`run_wave`、`example_*`）

### 第三轮 — 仓库卫生和文档（P2）

11. 更新 `.gitignore` 忽略二进制；`git rm --cached` 被跟踪的 3 个二进制
12. 删除 `mythread/` 下的陈旧副本（`test_bindcpu.c`、`test_field.c`）
13. 清理死代码/死字段/死配置
14. 按 simplify 后代码重写 WAVE_PROPAGATION_GHOST_DOC.md
15. 更新 benchmark 脚本或恢复能量输出
16. 决定是否删除 5 个 omp/mpi 变体源文件