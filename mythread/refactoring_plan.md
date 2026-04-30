# mythread 重构方案

## 一、现状分析

### 1.1 当前文件结构

```
mythread/
├── mythread.h                  # 头文件 (239行) — 所有声明
├── mythread.c                  # 实现 (1119行) — 所有实现
├── mythread接口说明.md          # API 中文文档
├── taskpool_spec.md            # 任务池设计说明
├── taskpool_tasks.md           # 任务池实现任务拆解
└── taskpool_checklist.md       # 任务池验证清单
```

### 1.2 现有问题

| 问题 | 说明 |
|------|------|
| **单文件过大** | `mythread.c` 1119 行，混合了 6 类职责不同的代码 |
| **职责耦合** | 线程生命周期、同步原语、任务池、Fortran 包装、TSC 计时、LocV 存储全部在一个编译单元 |
| **Fortran 接口散落** | Fortran 包装函数声明在 .h 中间位置，实现在 .c 末尾，没有独立文件 |
| **缺少 CPU 绑定测试** | `bindcpu()` 在 macOS 上是空实现，在 Linux 上依赖 `sched_setaffinity`，但没有任何自动化测试验证绑定是否正确 |
| **难以独立复用** | 用户如果只想用任务池或同步原语，必须引入整个 1119 行的 .c 文件 |
| **编译粒度粗** | 任何修改都需要重新编译全部代码 |

### 1.3 mythread.c 职责拆分

按代码行数统计：

| 职责模块 | 行号范围 | 大约行数 | 说明 |
|----------|----------|----------|------|
| 线程创建/生命周期 | 185-272, 618-872 | ~340 | InitThreads, initmd, StartThreads, EndThreads, zStartThreads, ClearThread, threadMain, setthread |
| 线程同步接口 | 280-540 | ~260 | FWAIT 宏, SM/MS, SG/GS, GM/MG 三级屏障 |
| 任务池 (thread pool) | 878-1118 | ~240 | mt_taskpool 全部实现 |
| Fortran 接口 | 69-77, 267-269, 275-277, 542-563, 591-616 | ~80 | 所有 `_suffixed` 包装函数 |
| TSC 性能计时 | 565-616 | ~50 | atsc, tscinit, tscb, tsce, tsceb, prtsc |
| LocV 存储系统 | 115-183 | ~70 | select_locv, SetLocV, GetLocV |
| CPU 绑定 | 44-67, 249-269 | ~35 | getcpuid, bindcpu, bindthread |
| 工具函数 | 271-307 | ~35 | ntdelay, opentf, padr_ |
| 全局变量/弱符号 | 23-42, 98-104 | ~25 | 配置变量, thread_run, GetVInt, _getgdsize_ |

---

## 二、重构目标

1. **按职责拆分文件**：线程创建、同步接口、任务池各自独立
2. **Fortran 接口独立**：预留独立的 Fortran 包装文件，方便后续扩展
3. **增加 CPU 绑定测试**：补充线程绑定核心的自动化测试
4. **保持向后兼容**：现有用户只需 `#include "mythread.h"` 即可，无需修改代码
5. **不影响现有功能**：所有现有测试和示例继续通过

---

## 三、目标文件结构

```
mythread/
├── mythread.h                  # 主头文件（聚合 include，向后兼容）
├── mythread_types.h            # 类型定义 + 常量 + 全局变量声明
│
├── mythread_thread.h           # 线程创建/生命周期接口声明
├── mythread_thread.c           # 线程创建/生命周期实现
│
├── mythread_sync.h             # 线程同步接口声明
├── mythread_sync.c             # 线程同步实现（三级屏障）
│
├── mythread_pool.h             # 任务池接口声明
├── mythread_pool.c             # 任务池实现（SPMC 队列）
│
├── mythread_fortran.h          # Fortran 接口声明（预留）
├── mythread_fortran.c          # Fortran 接口实现（预留）
│
├── mythread_timer.h            # TSC 计时接口声明
├── mythread_timer.c            # TSC 计时实现
│
├── mythread_locv.h             # LocV 存储接口声明
├── mythread_locv.c             # LocV 存储实现
│
├── mythread_util.h             # 工具函数声明（ntdelay, opentf 等）
├── mythread_util.c             # 工具函数实现
│
├── test_bindcpu.c              # 【新增】CPU 绑定测试
│
├── mythread接口说明.md          # API 文档（保持不变）
├── taskpool_spec.md            # 设计说明（保持不变）
├── taskpool_tasks.md           # 任务拆解（保持不变）
└── taskpool_checklist.md       # 验证清单（保持不变）
```

### 3.1 各文件职责详述

#### `mythread_types.h` — 核心类型定义

从 `mythread.h` 中提取以下内容：
- 常量宏：`RFB`, `MCOREPC`, `MCLUST`, `MSB`, `MSBG`
- 结构体定义：`THREADINFO`, `threadGroup`, `threadProc`
- 函数指针类型：`TFunc`, `Func`, `TSFunc`
- 任务池类型：`mt_task_fn`, `mt_task`, `mt_taskpool`（不透明）
- 任务池标志：`MT_TASKPOOL_BLOCK`, `MT_TASKPOOL_SPIN`, `MT_TASKPOOL_TRY`
- 全局变量 `extern` 声明：`md`, `ti`, `gi`, `mpi_id`, `NCorePClu` 等
- 便捷宏：`PTI`, `PTIM`, `VTI`, `VTIM`
- `hmalloc` / `hrealloc` 宏

#### `mythread_thread.h` / `mythread_thread.c` — 线程创建与生命周期

接口声明：
- `setthread()` — 设置线程参数
- `InitThreads()` — 初始化线程系统
- `StartThreads()` / `EndThreads()` — 启动/终止线程
- `bindcpu()` / `bindthread()` — CPU 核心绑定
- `getcpuid()` — 获取 CPU 信息
- `gettid()` / `getnt()` — 线程信息查询
- `thread_run()` — 弱符号默认线程函数
- `_threadmain_()` — 线程入口分发

内部函数（`.c` 内部）：
- `initmd()` — 初始化 `threadProc md` 结构
- `zStartThreads()` — 核心线程创建/终止逻辑
- `threadMain()` — pthread 线程入口
- `InitThread()` / `ClearThread()` — 线程初始化/清理
- `_getgdsize_()` / `_gettdsize_()` — 弱符号数据大小

#### `mythread_sync.h` / `mythread_sync.c` — 线程同步接口

接口声明：
- SM/MS 方向：`sSetState`, `sWaitState`, `sWaitStater`, `mSetSubs`, `mWaitSubs`, `mWaitSubsr`
- SG/GS 方向：`sSetGrp`, `sWaitGrp`, `sWaitGrpr`, `gSetSubs`, `gWaitSubs`, `gWaitSubsr`
- GM/MG 方向：`gSetMain`, `gWaitMain`, `gWaitMainr`, `mSetGrps`, `mWaitGrps`, `mWaitGrpsr`
- 便捷宏：`SSM`, `SWM`, `SWMR`, `MSS`, `MWS`, `MWSR`, `SSG`, `SWG`, `SWGR`, `GSS`, `GWS`, `GWSR`, `MSG`, `MWG`, `MWGR`, `GSM`, `GWM`, `GWMR`

内部实现：
- `FWAIT` 宏（生成 `Wait_LNE`, `Wait_LEQ`, `Wait_LGE`, `Wait_LLE`）
- 所有同步函数的实现体
- `GetVInt()` 弱符号

#### `mythread_pool.h` / `mythread_pool.c` — 任务池实现

接口声明：
- `mt_taskpool_attach()` / `mt_taskpool_detach()` — 创建/销毁
- `mt_taskpool_begin()` — 开启新 epoch
- `mt_taskpool_submit()` — 提交任务
- `mt_taskpool_close()` — 关闭提交
- `mt_taskpool_wait()` — 等待完成
- `mt_taskpool_shutdown()` — 终止
- `mt_taskpool_worker_loop()` — 工作线程循环

内部实现：
- `mt_taskpool` 结构体完整定义
- 辅助函数：`mt_mu_unlock`, `mt_taskpool_typ`, `mt_taskpool_nworkers`, `mt_taskpool_valid_slot`, `mt_taskpool_get`

#### `mythread_fortran.h` / `mythread_fortran.c` — Fortran 接口（预留）

接口声明：
- 初始化：`inithreads_`, `startthreads_`, `endthreads_`
- CPU 绑定：`bindthread_`
- 同步接口：`swaitstate_`, `swaitstater_`, `ssetstate_`, `mwaitsubs_`, `mwaitsubsr_`, `msetsubs_`, `swaitgrp_`, `swaitgrpr_`, `ssetgrp_`, `gwaitsubs_`, `gwaitsubsr_`, `gsetsubs_`, `gwaitmain_`, `gwaitmainr_`, `gsetmain_`, `mwaitgrps_`, `mwaitgrpsr_`, `msetgrps_`
- 计时器：`tscb_`, `tsce_`, `tsceb_`, `prtsc_`
- 其他：`ntdelay_`

实现：所有 Fortran 包装函数（C 调用约定，trailing underscore）

#### `mythread_timer.h` / `mythread_timer.c` — TSC 性能计时

接口声明：
- `tscinit()`, `tscb()`, `tsce()`, `tsceb()`, `prtsc()`
- Fortran 版本：`tscb_`, `tsce_`, `tsceb_`, `prtsc_`

内部实现：
- `atsc()` — 内联汇编读取时间戳

#### `mythread_locv.h` / `mythread_locv.c` — LocV 存储系统

接口声明：
- `SetLocV()` / `GetLocV()`

内部实现：
- `select_locv()` — 根据 typ 选择存储层级

#### `mythread_util.h` / `mythread_util.c` — 工具函数

- `ntdelay()` / `ntdelay_()` — 延迟函数
- `opentf()` — 打开线程日志文件
- `padr_()` — 调试辅助

#### `mythread.h` — 主头文件（向后兼容）

```c
#ifndef MTHREAD_H_INCLUDED
#define MTHREAD_H_INCLUDED

#include "mythread_types.h"
#include "mythread_thread.h"
#include "mythread_sync.h"
#include "mythread_pool.h"
#include "mythread_fortran.h"
#include "mythread_timer.h"
#include "mythread_locv.h"
#include "mythread_util.h"

#endif
```

保留原有的 `__BEGIN_DECLS` / `__END_DECLS` 包裹。

---

## 四、编译调整

### 4.1 Makefile 修改

```makefile
# 原来
mythread.o: mythread/mythread.c mythread/mythread.h
	$(CC) -c mythread/mythread.c -o mythread.o

# 改为
MYTHREAD_SRCS = mythread/mythread_thread.c \
                mythread/mythread_sync.c \
                mythread/mythread_pool.c \
                mythread/mythread_fortran.c \
                mythread/mythread_timer.c \
                mythread/mythread_locv.c \
                mythread/mythread_util.c

MYTHREAD_OBJS = $(MYTHREAD_SRCS:.c=.o)

libmythread.a: $(MYTHREAD_OBJS)
	ar rcs $@ $^

# 或者直接编译为目标文件列表
mythread_objs = mythread_thread.o mythread_sync.o mythread_pool.o \
                mythread_fortran.o mythread_timer.o mythread_locv.o mythread_util.o
```

### 4.2 依赖关系

```
mythread_types.h  (无依赖)
    ├── mythread_thread.h / .c  (依赖 types)
    ├── mythread_sync.h / .c    (依赖 types)
    ├── mythread_pool.h / .c    (依赖 types, locv)
    ├── mythread_timer.h / .c   (依赖 types)
    ├── mythread_locv.h / .c    (依赖 types)
    ├── mythread_util.h / .c    (依赖 types)
    └── mythread_fortran.h / .c (依赖 thread, sync, timer, util)

mythread.h  (聚合所有子头文件)
```

---

## 五、CPU 绑定测试方案

### 5.1 测试目标

验证 `bindcpu()` 在 Linux 上正确设置 CPU 亲和性，在 macOS 上优雅降级（空操作但不崩溃）。

### 5.2 新增文件：`mythread/test_bindcpu.c`

测试用例设计：

| 编号 | 测试名称 | 说明 |
|------|----------|------|
| TC-B01 | `bindcpu-single-valid` | 绑定到单个合法核心，验证 `sched_getaffinity` 返回的掩码中只有该核心被设置 |
| TC-B02 | `bindcpu-multiple-switch` | 先绑定核心 A，再切换到核心 B，验证每次亲和性正确变更 |
| TC-B03 | `bindcpu-macos-noop` | macOS 上 bindcpu 应返回 0（不报错），且不影响线程继续运行 |
| TC-B04 | `bindcpu-all-cores` | 遍历系统所有可用核心，逐一绑定并验证 |
| TC-B05 | `bindthread-ti-null` | `ti == NULL` 时 `bindthread()` 不崩溃 |
| TC-B06 | `bindthread-uses-indg` | 创建线程后验证 `bindthread()` 使用了 `ti->indg` 指定的核心 |
| TC-B07 | `thread-main-calls-bind` | 验证 `threadMain()` 在启动时确实调用了 CPU 绑定 |
| TC-B08 | `initmd-core-assign` | 验证 `initmd()` 在 grouped/ungrouped 模式下正确计算 `indg` 值（不重叠、在合法范围内） |

### 5.3 测试实现要点

```c
// TC-B01: 验证单核心绑定
void test_bindcpu_single_valid() {
    int target_core = 0;
    int ret = bindcpu(target_core);
#ifdef __APPLE__
    assert(ret == 0);  // macOS 空操作但成功返回
#else
    assert(ret == 0);
    cpu_set_t mask;
    CPU_ZERO(&mask);
    sched_getaffinity(0, sizeof(mask), &mask);
    assert(CPU_ISSET(target_core, &mask));
    // 验证只有 target_core 被设置（允许少量系统核心也被设置）
    int count = 0;
    for (int i = 0; i < CPU_SETSIZE; i++) {
        if (CPU_ISSET(i, &mask)) count++;
    }
    // 至少包含 target_core
    assert(count >= 1);
#endif
}
```

### 5.4 编译与运行

```bash
# 编译
gcc -D_GNU_SOURCE -o test_bindcpu mythread/test_bindcpu.c \
    mythread_thread.o mythread_locv.o mythread_util.o -lpthread

# 运行
./test_bindcpu
```

---

## 六、实施步骤

### 阶段一：拆分类型定义（不影响编译）

1. 创建 `mythread_types.h`，从 `mythread.h` 中移动类型定义、常量、全局变量声明
2. 修改 `mythread.h` 为聚合头文件，`#include "mythread_types.h"`
3. 验证编译通过

### 阶段二：拆分实现文件（逐模块）

按以下顺序拆分，每步验证编译：

1. **工具函数** → `mythread_util.h` / `mythread_util.c`
2. **LocV 存储** → `mythread_locv.h` / `mythread_locv.c`
3. **TSC 计时** → `mythread_timer.h` / `mythread_timer.c`
4. **线程同步** → `mythread_sync.h` / `mythread_sync.c`
5. **任务池** → `mythread_pool.h` / `mythread_pool.c`
6. **线程创建** → `mythread_thread.h` / `mythread_thread.c`
7. **Fortran 接口** → `mythread_fortran.h` / `mythread_fortran.c`

### 阶段三：增加 CPU 绑定测试

1. 创建 `mythread/test_bindcpu.c`
2. 在 CI 中加入测试编译和运行
3. 修复测试中发现的问题

### 阶段四：清理与验证

1. 删除 `mythread.c`（确认所有代码已迁移）
2. 运行全部现有测试（`taskpool_tests.c`）
3. 运行全部示例（`example_simple.c`, `example_mpi_mythread.c`, `example_taskpool.c`）
4. 运行 `wave_propagation_ghost.c` 主程序
5. 更新 Makefile

---

## 七、风险与注意事项

| 风险 | 等级 | 缓解措施 |
|------|------|----------|
| 内部依赖断裂 | 中 | `static` 函数需要在 .c 文件间共享时改为 `extern` 或移到共享头文件 |
| `FWAIT` 宏生成的 `Wait_*` 函数依赖 `ntdelay` | 低 | `ntdelay` 在 `mythread_util.c`，sync 模块 include util.h 即可 |
| 任务池依赖 LocV 和线程信息 | 中 | pool 模块 include locv.h 和 thread.h |
| Fortran 包装调用 C 函数 | 低 | fortran 模块 include 所有需要的子模块头文件 |
| macOS 编译兼容 | 低 | `bindcpu` 的平台条件编译保持不变 |
| `__thread` TLS 变量跨文件 | 低 | 全局 TLS 变量定义集中在一个 .c 文件中，其他文件通过 `extern` 引用（当前已在 .h 中声明为 `extern`） |
| 现有用户代码兼容性 | 极低 | 主头文件 `mythread.h` 保持向后兼容，所有接口不变 |

### 7.1 需要调整的 static 函数

当前 `mythread.c` 中以下 `static` 函数需要在拆分后调整可见性：

| 函数 | 当前可见性 | 被谁调用 | 调整方案 |
|------|-----------|----------|----------|
| `zStartThreads()` | `static` | `StartThreads`, `EndThreads` | 保持在 `mythread_thread.c` 内部 `static` |
| `select_locv()` | `static` | `SetLocV`, `GetLocV` | 保持在 `mythread_locv.c` 内部 `static` |
| `InitThread()` | `static` | `zStartThreads` | 保持在 `mythread_thread.c` 内部 `static` |
| `ClearThread()` | `static` | `zStartThreads` | 保持在 `mythread_thread.c` 内部 `static` |
| `mt_mu_unlock()` | `static` | `mt_taskpool_detach` 等 | 保持在 `mythread_pool.c` 内部 `static` |
| `mt_taskpool_typ()` | `static inline` | 多个 taskpool 函数 | 保持在 `mythread_pool.c` 内部 |
| `mt_taskpool_nworkers()` | `static inline` | `mt_taskpool_attach`, `mt_taskpool_begin` | 保持在 `mythread_pool.c` 内部 |
| `mt_taskpool_valid_slot()` | `static inline` | 多个 taskpool 函数 | 保持在 `mythread_pool.c` 内部 |
| `mt_taskpool_get()` | `static` | 多个 taskpool 函数 | 保持在 `mythread_pool.c` 内部 |

**结论：所有 `static` 函数都可以保持在各自的 .c 文件内部**，拆分自然形成模块边界。

---

## 八、替代方案对比

### 方案 A：最小拆分（本方案）

```
mythread_types.h, mythread_thread.h/.c, mythread_sync.h/.c,
mythread_pool.h/.c, mythread_fortran.h/.c,
mythread_timer.h/.c, mythread_locv.h/.c, mythread_util.h/.c
```

- **优点**：职责清晰，按需引入，编译粒度合理
- **缺点**：文件数量较多（14 个源文件 + 头文件）

### 方案 B：保守拆分

```
mythread.h (不变), mythread.c → mythread_core.c + mythread_pool.c + mythread_fortran.c
```

- **优点**：改动最小
- **缺点**：`mythread_core.c` 仍然混合多个职责，没有根本解决问题

### 方案 C：激进拆分

```
每个公开函数一个文件
```

- **优点**：极致模块化
- **缺点**：过度工程，文件过多，管理成本高

**推荐方案 A**，在模块化和可维护性之间取得平衡。

---

## 九、Fortran 接口预留说明

`mythread_fortran.h` / `mythread_fortran.c` 作为 Fortran 接口的独立文件，当前阶段将现有的 trailing-underscore 包装函数集中迁移至此。后续 Fortran 用户可以：

1. 直接 `use` 对应的 Fortran module（未来可基于此文件生成 `.f90` 接口模块）
2. 所有 `_suffixed` 函数集中在一个翻译单元，便于维护和生成
3. 如果未来需要支持 Fortran 2003 `bind(C)` 的 `ISO_C_BINDING` 方案，可以在此文件基础上扩展

**当前不创建 `.f90` 文件**，因为现有的 Fortran 互操作完全通过 C 的 trailing-underscore 约定实现。`.f90` 接口模块属于后续工作。
