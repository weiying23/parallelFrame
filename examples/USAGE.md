# mythread + MPI 使用指南

## 目录
1. [架构概述](#架构概述)
2. [编译运行](#编译运行)
3. [API 说明](#api-说明)
4. [示例代码解析](#示例代码解析)
5. [注意事项](#注意事项)

---

## 架构概述

```
┌─────────────────────────────────────────────────────────────┐
│                        节点 (Node)                           │
│  ┌─────────────────────────────────────────────────────┐   │
│  │                    MPI 进程 0                        │   │
│  │  ┌─────────────────────────────────────────────┐   │   │
│  │  │              线程组 0 (Group 0)              │   │   │
│  │  │  ┌─────┐ ┌─────┐ ┌─────┐ ┌─────┐ ┌─────┐  │   │   │
│  │  │  │Main │ │ T1  │ │ T2  │ │ T3  │ │ T4  │  │   │   │
│  │  │  │(GM) │ │(SWG)│ │(SWG)│ │(SWG)│ │(SWG)│  │   │   │
│  │  │  └─────┘ └─────┘ └─────┘ └─────┘ └─────┘  │   │   │
│  │  └─────────────────────────────────────────────┘   │   │
│  │  ┌─────────────────────────────────────────────┐   │   │
│  │  │              线程组 1 (Group 1)              │   │   │
│  │  │  ┌─────┐ ┌─────┐ ┌─────┐ ┌─────┐ ┌─────┐  │   │   │
│  │  │  │Main │ │ T1  │ │ T2  │ │ T3  │ │ T4  │  │   │   │
│  │  │  │(GM) │ │(SWG)│ │(SWG)│ │(SWG)│ │(SWG)│  │   │   │
│  │  │  └─────┘ └─────┘ └─────┘ └─────┘ └─────┘  │   │   │
│  │  └─────────────────────────────────────────────┘   │   │
│  │              [进程内通信: MSG/MWG]                  │   │
│  └─────────────────────────────────────────────────────┘   │
│                            │                                │
│                      [MPI 通信]                             │
│                            │                                │
│  ┌─────────────────────────────────────────────────────┐   │
│  │                    MPI 进程 1                        │   │
│  │                    [同上结构]                        │   │
│  └─────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### 线程层级

1. **Main Thread** (tid=0, igrp=0): 主线程，管理所有组
2. **Group Main** (tid=0, igrp>0): 组主线程，管理组内工作线程
3. **Worker Thread** (tid>0): 工作线程，执行实际计算

### 同步机制

| 层级 | 方向 | 宏 | 说明 |
|------|------|-----|------|
| SM/MS | 子线程↔主线程 | SSM/SWM | 单组模式下的同步 |
| SG/GS | 子线程↔组主 | SSG/SWG | 组内工作线程同步 |
| GM/MG | 组主↔主线程 | GSM/GWM | 跨组同步 |

---

## 编译运行

### 编译

```bash
make all
```

或手动编译：

```bash
# 简化示例
mpicc -Wall -O2 -o example_simple example_simple.c mythread/mythread.c -lpthread -lm

# 完整示例
mpicc -Wall -O2 -o example_mpi_mythread example_mpi_mythread.c mythread/mythread.c -lpthread -lm
```

### 运行

```bash
# 运行简化示例（2个MPI进程）
make run_simple
# 或
mpirun -np 2 ./example_simple

# 运行完整示例（4个MPI进程）
make run_full
# 或
mpirun -np 4 ./example_mpi_mythread
```

---

## API 说明

### 必需函数

在使用 mythread 前，必须定义以下两个函数：

```c
// 返回线程私有数据大小
int _gettdsize_() {
    return sizeof(MyThreadData);
}

// 返回组共享数据大小
int _getgdsize_() {
    return sizeof(MyGroupData);
}
```

### 初始化函数

```c
int InitThreads(int mpi_id_,       // MPI 进程ID
                int NCorePClu_,    // 每簇CPU核心数
                int NCluPNode_,    // 每节点簇数
                int NCorePGrp_,    // 每组CPU核心数
                int NThPGrp_,      // 每组线程数（含主线程）
                int NGrpPProc_,    // 每进程组数
                int NProcPNode_,   // 每节点进程数
                int *ManageCoreId_ // 管理核心ID（输出）
               );
```

### 线程控制

```c
// 启动所有工作线程
void StartThreads(TFunc tfun);

// 结束所有线程
void EndThreads();
```

### 信息获取

```c
// 获取当前线程ID (0=主线程)
int gettid();

// 获取组内线程总数
int getnt();

// 获取MPI进程ID（全局变量）
extern int mpi_id;
```

### 同步宏

```c
// 主线程与所有组同步
MSG;    // 主线程设置组状态（通知所有组开始）
MWG;    // 主线程等待所有组完成

// 组主线程与主线程同步
GWM;    // 组主线程等待主线程信号
GSM;    // 组主线程通知主线程完成

// 组内工作线程与组主线程同步
SWG;    // 工作线程等待组主信号
SSG;    // 工作线程通知组主完成

// 组主线程与组内工作线程同步
GSS;    // 组主通知所有工作线程开始
GWS;    // 组主等待所有工作线程完成
```

---

## 示例代码解析

### 基本使用流程

```c
#include <mpi.h>
#include "mythread/mythread.h"

// 1. 定义数据大小函数
int _gettdsize_() { return 64; }
int _getgdsize_() { return 64; }

// 2. 定义线程函数
void worker_func() {
    int tid = gettid();
    
    // 等待组主信号
    SWG;
    
    // 执行计算...
    printf("Worker %d working\n", tid);
    
    // 通知完成
    SSG;
}

void thread_entry() {
    int tid = gettid();
    int gid = ti->igrp;
    
    if (tid == 0) {
        if (ThreadG && gid > 0) {
            // 组主线程
            GWM;    // 等待主线程
            GSS;    // 通知工作线程
            GWS;    // 等待工作线程
            GSM;    // 通知主线程
        } else {
            // 主线程
            MSG;    // 通知所有组
            MWG;    // 等待所有组
        }
    } else {
        // 工作线程
        worker_func();
    }
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    // 3. 初始化线程
    int manage_id = -1;
    InitThreads(rank, 16, 2, 16, 5, 2, size, &manage_id);
    
    // 4. 启动线程
    StartThreads(thread_entry);
    
    // 5. 主线程也执行入口函数
    thread_entry();
    
    // 6. 结束线程
    EndThreads();
    
    // 7. MPI 通信
    MPI_Barrier(MPI_COMM_WORLD);
    
    MPI_Finalize();
    return 0;
}
```

### 数据分配示例

```c
void parallel_compute(double *data, int n, double *result) {
    int tid = gettid();
    int nt = getnt();
    
    // 计算本线程的数据范围
    int chunk = n / nt;
    int start = tid * chunk;
    int end = (tid == nt - 1) ? n : start + chunk;
    
    // 执行计算
    for (int i = start; i < end; i++) {
        result[i] = compute(data[i]);
    }
}
```

### 进程间通信

```c
// 在所有线程完成后进行MPI通信
void mpi_communication() {
    // 准备本地数据
    double local_sum = compute_local_sum();
    double global_sum = 0;
    
    // 归约到所有进程
    MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, 
                  MPI_SUM, MPI_COMM_WORLD);
    
    // 或使用 Gather
    double all_sums[size];
    MPI_Gather(&local_sum, 1, MPI_DOUBLE, all_sums, 1, MPI_DOUBLE,
               0, MPI_COMM_WORLD);
}
```

---

## 注意事项

### 1. 线程安全

- **MPI 调用**: 应在主线程中进行 MPI 调用，或使用 MPI_THREAD_MULTIPLE 初始化
- **全局变量**: 注意线程间的数据竞争，使用线程本地存储

### 2. 线程绑定

```c
// 线程会自动绑定到指定CPU核心
// 可通过 ManageCoreId 参数控制管理线程的位置
```

### 3. 调试

```c
// 编译时添加 -DDEBUG 启用调试日志
make debug

// 查看生成的日志文件
// MTW_XX_XX_XX.log  (MPI进程_组ID_线程ID)
```

### 4. 性能建议

1. **负载均衡**: 确保各线程工作量相近
2. **减少同步**: 尽量减少同步点数量
3. **数据局部性**: 利用线程本地数据减少竞争
4. **CPU亲和性**: 合理的核心分配可减少缓存失效

### 5. 常见问题

| 问题 | 可能原因 | 解决方案 |
|------|----------|----------|
| InitThreads 返回错误 | 参数超出范围 | 检查 NThPGrp <= MCOREPC |
| 死锁 | 同步宏不匹配 | 确保每个 wait 都有对应的 set |
| 段错误 | 数据大小函数未定义 | 实现 _gettdsize_ 和 _getgdsize_ |

---

## 高级用法

### 自定义数据

```c
typedef struct {
    double *input;
    double *output;
    int start_idx;
    int end_idx;
} ThreadData;

int _gettdsize_() {
    return sizeof(ThreadData);
}

void thread_func() {
    ThreadData *td = (ThreadData*)ti->td;  // 获取线程私有数据
    // 使用 td->input 等...
}
```

### 多轮迭代

```c
for (int iter = 0; iter < max_iter; iter++) {
    // 通知所有组开始
    MSG;
    MWG;
    
    // 处理结果
    process_results();
    
    // 准备下一轮
    MSS;  // 重置子线程状态
}
```
