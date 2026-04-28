# mythread 组内线程池（任务队列）管理方案设计

## 背景

当前 mythread 的典型使用方式是“轮次/阶段式”同步：主线程（或组主线程）以 `MSG/GSS` 发出开始信号，从线程以 `SWG/SSG` 或 `GWS` 等原语完成一次并行阶段后再汇合。该模型适合每轮任务结构固定、可静态切分的场景；当任务粒度不均或需要动态产生任务（例如自适应网格、稀疏遍历、队列驱动计算）时，静态切分容易出现负载不均、尾部拖慢。

本设计在保留现有同步模型与线程绑定策略的前提下，引入“组内线程池 + 任务队列”：由组主线程维护任务队列，从线程循环从队列取任务并执行，直至该轮任务关闭且队列清空。

## 目标

- 在每个线程组（`threadGroup`）内部提供任务队列与 worker 执行循环，实现动态负载均衡。
- 保持 mythread 的线程亲和力/绑核逻辑不变，队列机制不改变线程创建销毁语义。
- 支持“按轮次提交任务并等待完成”的用法，与 `GSS/GWS`、`SWG/SSG` 等同步风格可组合。
- 提供可控的背压策略（队列满时阻塞/自旋/返回失败），避免无界内存增长。
- 明确任务的内存所有权、生命周期与取消/退出行为，避免死锁与悬挂指针。

## 非目标

- 不做跨组 work-stealing（组间窃取）与全局调度。
- 不做任务优先级、多队列、多生产者并发入队（默认单生产者：组主线程）。
- 不保证任务可抢占或强制中断；取消仅在轮次边界/退出协议下生效。

## 术语与角色

- 进程主线程（Main Main Thread，简称 MMT）：`ThreadG=1` 时负责跨组协调（示例中 `MSG/MWG`）。
- 组主线程（Group Main Thread，简称 GMT）：每组 `tid==0` 的线程，负责本组队列管理与轮次控制。
- 从线程/工作线程（Sub/Worker Thread）：组内 `tid>0` 的线程，执行任务。
- 任务（Task）：可被 worker 执行的最小工作单元，包含函数指针与上下文指针。
- 轮次（Epoch）：一次“提交任务→执行→等待完成”的闭环，用单调递增的 epoch 标识，防止跨轮信号串扰。

## 总体架构

### 每组一条任务队列

- 任务队列绑定到 `threadGroup`，生命周期与线程组一致。
- 推荐把队列指针存放在组共享位置，优先选择：
  - `gi->locv[k]`（通过 `SetLocV(1, k, ptr)`）存放队列对象指针；或
  - `gi->gd` 指向的一段组共享内存里内嵌队列头部结构（适合 Fortran/跨语言共享）。

### 单生产者 + 多消费者（SPMC）

- 生产者：GMT（每组仅一个）负责入队与关闭本轮任务。
- 消费者：组内 worker 多线程并发出队执行。
- 该约束允许使用更轻量的无锁/低锁实现（环形队列 + 原子 head）。

### 轮次驱动的执行协议

每个 epoch 由 GMT 驱动：

1. GMT 开启 epoch（`epoch++`），将任务入队。
2. GMT 广播“epoch 开始”（可复用 `GSS` 或独立的 `queue_epoch` 状态）。
3. Worker 被唤醒后循环出队执行，直到同时满足：
   - 队列为空；且
   - GMT 已关闭本 epoch（`open=0`）。
4. Worker 上报“本 epoch 完成”（可复用 `SSG` 或独立的 `done_epoch` 状态）。
5. GMT 等待所有 worker 完成（可复用 `GWS`），收敛后进入下一个 epoch 或退出。

该协议与现有 mythread 的“阶段式同步”一致，只是把阶段内部从“静态划分”替换为“动态取任务”。

## 数据结构（建议）

### 任务描述

- `fn`：任务函数指针
- `ctx`：任务上下文指针（由提交方分配/管理）
- `tag/arg`（可选）：用于轻量参数传递或调试统计

任务函数签名建议二选一：

1. 最小签名（通过 `ti/gi` TLS 获取线程信息）  
   `typedef void (*mt_task_fn)(void *ctx);`
2. 显式携带线程上下文（便于单测/复用）  
   `typedef void (*mt_task_fn2)(void *ctx, THREADINFO *ti, threadGroup *gi);`

### 环形队列（固定容量）

- `capacity`：固定容量，建议为 2 的幂，方便取模。
- `tail`：生产者写入位置（仅 GMT 修改）。
- `head`：消费者读取位置（worker 并发修改，需原子/锁保护）。
- `open`：本 epoch 是否仍允许继续入队/执行（GMT 关闭）。
- `epoch`：当前轮次号，单调递增。

并发策略建议：

- 低锁策略：使用 C11 原子（`stdatomic.h`）实现 `head` 的 `fetch_add`，`tail` 由单生产者写；读写屏障保证任务可见性。
- 保守策略：`pthread_mutex + pthread_cond` 保护队列（入队 signal，出队 wait），实现简单且避免自旋，但开销略高。

在 HPC 场景下可默认“低锁 + 轻量自旋退避（调用 `ntdelay`）”，并提供编译期开关或运行时 flags 切换策略。

## API（对外行为约定）

下面 API 以“组内线程池”为中心，所有调用默认发生在组内线程（可通过 TLS `ti/gi` 获取上下文）。

### 生命周期

- `mt_taskpool_attach(int slot, int capacity, int flags)`  
  在当前组创建并挂载一个任务池到 `gi->locv[slot]`。GMT 调用一次即可；worker 只读访问。
- `mt_taskpool_detach(int slot)`  
  释放并从 `gi->locv[slot]` 解挂。要求该组无正在执行的 epoch。

### 提交与关闭

- `mt_taskpool_begin(int slot)`  
  GMT 开启新 epoch，重置统计并标记 `open=1`。
- `mt_taskpool_submit(int slot, mt_task_fn fn, void *ctx)`  
  GMT 提交任务。队列满时的行为由 `flags` 决定：
  - 阻塞直到有空间
  - 自旋直到有空间
  - 返回失败（由上层决定拆分/重试）
- `mt_taskpool_close(int slot)`  
  GMT 关闭本 epoch：不再入队，worker 仅需 drain 现有任务。

### 执行与等待

- `mt_taskpool_worker_loop(int slot)`  
  worker 执行循环：等待 epoch 开始信号，然后不断出队执行，直到 `open=0 && empty`，再上报完成并回到等待。
- `mt_taskpool_wait(int slot)`  
  GMT 等待所有 worker 完成本 epoch（内部可用 `GWS` 或自维护计数/状态）。

### 退出

- `mt_taskpool_shutdown(int slot)`  
  GMT 发出退出信号，使 worker loop 结束并返回到 `thread_run` 上层逻辑（便于 `EndThreads` 取消前干净退出）。

## 与现有 mythread 同步原语的结合建议

为了最大化复用已有同步机制，建议把“epoch 开始/完成”的信号复用为现有状态流：

- epoch 开始：GMT 在 `begin` 后调用一次 `GSS`（组主给从线程发信号），worker 的 `worker_loop` 用 `SWG` 等待。
- epoch 完成：worker 在完成 drain 后调用 `SSG`，GMT 用 `GWS` 等待全体完成。

这样，任务队列仅负责“阶段内部的动态分配”，阶段边界依然使用 mythread 的既有协议，避免引入第二套同步系统导致交织死锁。

## 任务上下文与内存所有权

- `ctx` 的内存由提交方负责：
  - 若 `ctx` 指向共享数组元素或结构体，需保证任务执行期间有效。
  - 若 `ctx` 指向临时栈变量，禁止提交（除非等待完成前栈帧始终有效）。
- 任务函数内部若访问组共享数据，需遵守上层的同步约束（例如 epoch 内可以无锁写不同片段，或自行加锁）。
- 推荐提供可选的统计字段（提交数、执行数、失败数、空转次数）用于调优，但不作为正确性依赖。

## 错误处理与边界条件

- 队列容量不足：
  - 若选择返回失败，上层可以把“大任务”拆小或改用静态切分。
  - 若选择阻塞/自旋，必须保证 worker 会持续出队，否则可能死锁（典型场景：GMT 在单线程里提交任务同时等待空间，但没有唤醒 worker 的 epoch 开始信号）。
- 关闭顺序：
  - `close` 必须在提交完所有任务后调用。
  - `wait` 必须在 `close` 之后调用，否则 wait 条件不成立（worker 可能一直等待 open 关闭）。
- 退出与取消：
  - 推荐先走 `shutdown` 让 worker loop 自然退出，再由上层 `EndThreads` 进入清理；避免 `pthread_cancel` 在持锁区/半更新队列时中断。

## 性能与资源约束

- 默认固定容量环形队列，避免每个任务入队都动态分配。
- 任务粒度建议保持在“每个任务至少数百到数千条指令”级别，避免调度开销占比过高。
- 若采用自旋退避，退避函数复用 `ntdelay`，避免把 CPU 100% 占满。

