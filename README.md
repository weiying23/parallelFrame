# parallelFrame

一个面向 HPC 的轻量级并行框架,核心是自研的 **mythread** 线程运行时,提供「组(group)→ 工作线程(worker)」两级线程模型、NUMA 感知的 CPU 绑定、进程内/跨进程 halo 通信,以及任务池(task pool)。框架以二维波动方程求解器作为端到端示范应用,同时提供 MPI、OpenMP、mythread 多种实现变体的对照。

## 特性

- **两级线程模型**:主线程 → 组主(Group Main)→ 工作线程(Worker),配套多层同步宏(`MSG/MWG`、`GSM/GWM`、`SSG/SWG` 等)。
- **NUMA 感知**:按 cluster / node 拓扑分配核心,支持 `NCorePClu`、`NCluPNode`、`NCorePGrp`、`ManageCoreId` 等参数精细控制线程绑定。
- **域分解**:`mythread_decomp` 支持 `Y_ONLY` 与 `XY_2D` 两种组内分解策略,自动生成带 halo 的子域分块。
- **Halo 通信**:区分进程内(`mythread_halo_exchange_intra`)与跨进程(`mythread_halo_exchange_mpi`)两套路径。
- **任务池**:`mt_taskpool_*` 系列接口,支持阻塞/非阻塞提交与多轮 barrier 同步。
- **配置驱动**:INI 风格的 `case.cfg` / `hardware.cfg`,可由环境变量覆盖路径。
- **示范应用齐全**:同一波动方程提供 NUMA-aware、纯 MPI、MPI+OpenMP、MPI+mythread、MPI+sync 多个变体,便于性能对照。

## 目录结构

```
parallelFrame/
├── CMakeLists.txt              构建入口
├── mythread/                   mythread 运行时库(源码 + 头文件)
│   ├── mythread.h              统一头文件
│   ├── mythread_*.c/.h         各模块(util/locv/timer/sync/pool/thread/fortran/decomp/field/halo/config)
│   └── *.md                    设计文档与接入指南
├── examples/                   使用示例(example_simple / example_mpi_mythread / example_taskpool)
├── tests/                      各模块单元测试(test_*、taskpool_tests)
├── config/                     示例配置(case.cfg / hardware.cfg 及 _bench 版本)
├── scripts/                    基准测试脚本(benchmark.sh / benchmark_full.sh)
├── wave_propagation_ghost*.c  二维波动方程求解器变体
├── wave_visual.c               小网格可视化变体(配合 visualize_wave.py)
└── visualize_wave.py           波场快照可视化脚本
```

## 编译

依赖:CMake ≥ 3.16、C 编译器、Pthreads。启用 MPI 时还需 MPI 实现(如 OpenMPI/MPICH);OpenMP 为可选项,存在则自动编译对应的 `*_omp` 目标。

```bash
# 默认构建(启用 MPI)
cmake -S . -B build
cmake --build build -j

# 不需要 MPI(仅构建非 MPI 目标:taskpool_tests、test_*、example_taskpool)
cmake -S . -B build -DUSE_MPI=OFF
cmake --build build -j
```

构建配置会打印 MPI / OpenMP / Threads 的探测结果。MPI 程序请用 `mpirun` 启动,例如:

```bash
mpirun -np 4 ./build/wave_propagation_ghost
```

## 配置

波动求解器通过两个 INI 文件驱动,路径可用环境变量覆盖:

| 环境变量 | 默认路径 | 说明 |
|---|---|---|
| `WAVE_CASE_CFG` | `config/case.cfg` | 算例参数:网格 `NX/NY`、步数 `NT`、振幅 `A`、步长 `DT`、波速 `C0`、halo 宽、能量诊断间隔等 |
| `WAVE_HARDWARE_CFG` | `config/hardware.cfg` | 硬件/线程参数:`N_GROUPS`、`N_WORKERS`、`NCorePClu`、`NCluPNode`、`NCorePGrp`、`ManageCoreId`、负载不均衡测试参数等 |

关键字段含义见各 `.cfg` 文件内注释。

## 构建目标

### 库
- **`mythread`**(静态库 `libmythread.a`):运行时核心,对外暴露 `mythread/mythread.h`。

### 非 MPI 目标(`USE_MPI=OFF` 也可构建)
`taskpool_tests`、`test_bindcpu`、`test_config`、`test_decomp`、`test_locv`、`test_sync`、`test_timer`、`test_util`、`example_taskpool`。

### MPI 目标(仅 `USE_MPI=ON`)
- 示例:`example_simple`、`example_mpi_mythread`、`test_field`
- 波动求解器变体:

| 目标 | 说明 |
|---|---|
| `wave_propagation_ghost` | NUMA-aware 版,使用完整 mythread API + 任务池分块计算 |
| `wave_propagation_ghost_fix` | 上述版本的重构版,改进 halo/内存校验 |
| `wave_propagation_ghost_mpi` | 最小纯 MPI 实现 |
| `wave_propagation_ghost_mpi_mythread` | MPI + mythread 工作池 |
| `wave_propagation_ghost_mpi_sync` | 多组 MPI,使用 `mythread_sync` 阶段同步 |
| `wave_propagation_ghost_omp` | OpenMP 版(检测到 OpenMP 时构建) |
| `wave_propagation_ghost_mpi_omp` | MPI + OpenMP 版(检测到 OpenMP 时构建) |

> 注:`wave_visual.c` 为独立可视化工具,未纳入 CMake,可手动编译:`cc wave_visual.c -lm -o wave_visual`。

## 快速上手

最小任务池示例(`examples/example_taskpool.c`):用 `InitThreads` 初始化线程组,通过 `mt_taskpool_submit`/`mt_taskpool_wait` 提交并等待任务,`StartThreads`/`EndThreads` 管理生命周期。完整 API 与同步宏语义见 [`examples/USAGE.md`](examples/USAGE.md) 与 [`mythread/应用接入指南.md`](mythread/应用接入指南.md)。

## 文档

- `examples/USAGE.md` — mythread + MPI 使用指南(架构、API、同步宏、示例解析)
- `mythread/应用接入指南.md` — 应用接入流程
- `mythread/mythread接口说明.md` — 接口说明
- `mythread/numa_per_group_design.md` — NUMA 分组设计
- `mythread/taskpool_spec.md` — 任务池规格
- `mythread/OMP到mythread转换指南.md` — 从 OpenMP 迁移

## 许可

详见仓库根目录许可证文件(如未附带,默认按仓库声明执行)。
