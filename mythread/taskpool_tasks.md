# mythread 组内线程池（任务队列）实现任务拆分

## 目标

按 [taskpool_spec.md](file:///c:/Users/write/Documents/projects/parallelFrame/parallelFrame/mythread/taskpool_spec.md) 交付一套可用的“组主维护任务队列、从线程取任务执行”的线程池能力，并可与现有 `GSS/GWS`、`SWG/SSG` 同步原语组合。

## 任务拆分

### 1) 代码现状梳理与落点选择

- 阅读并确认 `threadGroup/THREADINFO` 的共享数据存储方式（`gd/locv`）与 TLS 访问路径（`ti/gi`）。
- 选定任务池挂载位置：
  - 首选：`gi->locv[slot]` 存放 `mt_taskpool*`
  - 备选：`gi->gd` 内嵌头部（需要保证 `_getgdsize_()` 足够大）

### 2) 新增任务池数据结构与 API（头文件）

修改文件：
- [mythread.h](file:///c:/Users/write/Documents/projects/parallelFrame/parallelFrame/mythread/mythread.h)

新增内容（不改变既有接口）：
- `typedef`：任务函数类型、任务结构体、flags
- API 声明：
  - `mt_taskpool_attach/detach`
  - `mt_taskpool_begin/submit/close/wait`
  - `mt_taskpool_worker_loop`
  - `mt_taskpool_shutdown`

约束：
- 不引入 C++ 依赖，保持 C ABI；必要时用 `__BEGIN_DECLS/__END_DECLS` 包裹。

### 3) 任务池实现（.c）

修改文件：
- [mythread.c](file:///c:/Users/write/Documents/projects/parallelFrame/parallelFrame/mythread/mythread.c)

实现要点：
- 环形队列（固定容量）
- 并发模型 SPMC（单生产者 GMT，多消费者 worker）
- 同步策略选择其一并固化：
  - 保守：`pthread_mutex + pthread_cond`
  - 低锁：C11 原子 + 自旋退避（若编译器/平台支持）
- epoch 状态机与 shutdown 退出逻辑

注意：
- 避免在 `pthread_cancel` 可打断的路径里持锁；必要时对 worker loop 设定取消点策略（或完全不依赖 cancel）。
- 不记录日志、不新增注释（遵守仓库风格约束）。

### 4) 示例与用法文档

修改/新增文件：
- 新增 `examples/example_taskpool.c`（或扩展 `example_mpi_mythread.c`）
- 更新 [mythread接口说明.md](file:///c:/Users/write/Documents/projects/parallelFrame/parallelFrame/mythread/mythread接口说明.md) 增加“任务队列线程池”章节

示例场景：
- GMT 在每轮提交 N 个任务（比如处理不均匀工作量的 block 列表）
- worker 在 `SWG` 后进入 `mt_taskpool_worker_loop` 的“drain once”模式或循环模式
- GMT `close+wait` 收敛，然后进入下一轮

### 5) 回归与验证

- 编译验证：保证 examples 能编译链接（`-lpthread`）
- 行为验证：
  - 任务总数=执行总数（无丢失/重复）
  - 关闭后能正常退出/进入下一轮
  - 队列满时 flags 行为符合约定

## 风险点与缓解

- 死锁风险：GMT 在未唤醒 worker 前提交任务且队列满 → 通过约束“begin 后先广播开始，再允许阻塞 submit”或 submit 返回失败避免。
- 生命周期风险：ctx 指向栈内存 → 文档明确禁止；示例使用堆或全局数组。
- 性能风险：cond 过多唤醒 → 通过批量提交、一次 broadcast、worker 自旋退避混合减少系统调用。

