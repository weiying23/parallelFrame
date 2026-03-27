/**
 * MPI + mythread 示例程序
 * 
 * 功能：
 * 1. 通过 MPI 启动多进程
 * 2. 每个进程中启动线程组（支持分组模式）
 * 3. 线程组内并行计算
 * 4. 进程间通过 MPI 进行通信
 * 
 * 编译：mpicc -o example_mpi_mythread example_mpi_mythread.c mythread/mythread.c -lpthread -lm
 * 运行：mpirun -np 4 ./example_mpi_mythread
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "mythread/mythread.h"

// 全局数据
#define DATA_SIZE 1024
#define NGROUPS 2       // 每进程的线程组数
#define NTHREADS_PER_GROUP 4  // 每组线程数

typedef struct {
    double *data;
    double *result;
    int start;
    int end;
    int iter;
} ThreadData;

ThreadData g_thread_data;

// 线程本地数据大小
// int _gettdsize_() {
//     return sizeof(ThreadData);
// }

// // 组共享数据大小
// int _getgdsize_() {
//     return sizeof(double) * 16; // 组内共享数据
// }

/**
 * 工作线程函数 - 执行实际的并行计算
 */
void worker_thread_func() {
    int tid = gettid();           // 获取线程ID (0 是主线程)
    int nt = getnt();             // 获取组内线程总数
    
    // 获取线程本地数据
    ThreadData *td = (ThreadData*)ti->td;
    
    printf("[Worker] MPI=%d, Group=%d, Thread=%d/%d started\n", 
           mpi_id, ti->igrp, tid, nt);
    
    // 迭代计算示例
    for (int iter = 0; iter < 3; iter++) {
        // 等待组主线程发出开始信号 (SG: SubThread wait GroupMain)
        SWG;
        
        // 计算分配给本线程的数据范围
        int chunk_size = (td->end - td->start) / nt;
        int my_start = td->start + tid * chunk_size;
        int my_end = (tid == nt - 1) ? td->end : my_start + chunk_size;
        
        // 执行计算任务
        double local_sum = 0.0;
        for (int i = my_start; i < my_end; i++) {
            td->result[i] = sin(td->data[i]) * cos(td->data[i]);
            local_sum += td->result[i];
        }
        
        printf("[Worker] MPI=%d, Group=%d, Thread=%d: iter=%d, range=[%d,%d), local_sum=%.6f\n",
               mpi_id, ti->igrp, tid, iter, my_start, my_end, local_sum);
        
        // 通知组主线程计算完成 (SG: SubThread Set Group state)
        SSG;
    }
    
    printf("[Worker] MPI=%d, Group=%d, Thread=%d finished\n", 
           mpi_id, ti->igrp, tid);
}

/**
 * 组主线程函数 - 管理组内子线程
 */
void group_main_func() {
    int my_group = ti->igrp;
    
    printf("[GroupMain] MPI=%d, Group=%d started, managing %d threads\n", 
           mpi_id, my_group, ti->Nthreads);
    
    for (int iter = 0; iter < 3; iter++) {
        printf("[GroupMain] MPI=%d, Group=%d: Starting iteration %d\n", 
               mpi_id, my_group, iter);
        
        // 通知所有子线程开始计算 (GS: GroupMain Set Subs)
        GSS;
        
        // 等待所有子线程完成 (GS: GroupMain Wait Subs)
        GWS;
        
        printf("[GroupMain] MPI=%d, Group=%d: Iteration %d completed\n", 
               mpi_id, my_group, iter);
        
        // 通知主线程本组已完成 (GM: GroupMain Set Main)
        GSM;
        
        // 等待主线程的下一次迭代信号 (GM: GroupMain Wait Main)
        GWM;
    }
    
    printf("[GroupMain] MPI=%d, Group=%d finished\n", mpi_id, my_group);
}

/**
 * 主线程函数 - 管理所有组
 */
void main_thread_func() {
    printf("[Main] MPI=%d main thread started, managing %d groups\n", 
           mpi_id, md.ngrp);
    
    for (int iter = 0; iter < 3; iter++) {
        printf("[Main] MPI=%d: Starting iteration %d\n", mpi_id, iter);
        
        // 通知所有组主线程开始迭代 (MG: Main Set Groups)
        MSG;
        
        // 等待所有组完成 (MG: Main Wait Groups)
        MWG;
        
        printf("[Main] MPI=%d: All groups completed iteration %d\n", mpi_id, iter);
        
        // 组内归约：汇总各组结果
        double group_sums[NGROUPS] = {0};
        for (int g = 0; g < NGROUPS; g++) {
            // 这里可以读取各组的计算结果
            group_sums[g] = 1.0 * g; // 简化示例
        }
        
        // 进程内归约
        double local_total = 0;
        for (int g = 0; g < NGROUPS; g++) {
            local_total += group_sums[g];
        }
        
        printf("[Main] MPI=%d: Local total = %.6f\n", mpi_id, local_total);
        
        // 进程间通信：MPI 归约
        double global_total = 0;
        MPI_Allreduce(&local_total, &global_total, 1, MPI_DOUBLE, 
                      MPI_SUM, MPI_COMM_WORLD);
        
        printf("[Main] MPI=%d: Global total from all processes = %.6f\n", 
               mpi_id, global_total);
        
        // 设置子线程进行下一轮迭代 (MS: Main Set Subs)
        MSS;
    }
    
    printf("[Main] MPI=%d: All iterations completed\n", mpi_id);
}

/**
 * 线程入口函数 - 根据线程类型分发
 */
void thread_run() {
    int tid = gettid();
    
    if (tid == 0) {
        // 主线程 (tid=0)
        if (ThreadG && ti->igrp == 0) {
            // 分组模式下，igrp=0 的组主线程作为主管理线程
            main_thread_func();
        } else if (ThreadG) {
            // 其他组的主线程
            group_main_func();
        } else {
            // 非分组模式，主线程直接管理
            main_thread_func();
        }
    } else {
        // 工作线程
        worker_thread_func();
    }
}

int main(int argc, char *argv[]) {
    int mpi_rank, mpi_size;
    
    // MPI 初始化
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    
    printf("[MPI] Process %d/%d started\n", mpi_rank, mpi_size);
    
    // 线程配置参数
    int NCorePClu = 16;     // 每簇CPU核心数
    int NCluPNode = 2;      // 每节点簇数
    int NCorePGrp = 16;     // 每组CPU核心数
    int NThPGrp = NTHREADS_PER_GROUP + 1;  // 每组线程数（+1 包含主线程）
    int NGrpPProc = NGROUPS; // 每进程组数
    int NProcPNode = mpi_size; // 每节点进程数
    int ManageCoreId = -1;  // 管理核心ID (-1 表示自动选择)
    
    // 初始化数据
    g_thread_data.data = (double*)malloc(sizeof(double) * DATA_SIZE);
    g_thread_data.result = (double*)malloc(sizeof(double) * DATA_SIZE);
    g_thread_data.start = 0;
    g_thread_data.end = DATA_SIZE;
    
    // 初始化数据
    for (int i = 0; i < DATA_SIZE; i++) {
        g_thread_data.data[i] = (double)i / DATA_SIZE * 3.14159;
    }
    
    // 初始化线程系统
    int err = InitThreads(mpi_rank, NCorePClu, NCluPNode, NCorePGrp, 
                          NThPGrp, NGrpPProc, NProcPNode, &ManageCoreId);
    
    if (err != 0) {
        printf("[Error] MPI=%d: InitThreads failed with error %d\n", mpi_rank, err);
        MPI_Finalize();
        return 1;
    }
    
    printf("[MPI] Process %d: Thread system initialized\n", mpi_rank);
    printf("[MPI] Process %d: ManageCoreId = %d, ThreadG = %d, ngrp = %d\n", 
           mpi_rank, ManageCoreId, ThreadG, md.ngrp);
    
    // 启动线程
    StartThreads(thread_run);
    
    // 主进程也执行工作（主线程会调用 thread_entry）
    thread_run();
    
    // 结束线程
    EndThreads();
    
    printf("[MPI] Process %d: Threads ended\n", mpi_rank);
    
    // 进程间同步
    MPI_Barrier(MPI_COMM_WORLD);
    
    // 示例：进程间数据交换
    if (mpi_rank == 0) {
        printf("\n[MPI] All processes completed computation\n");
        
        // 收集各进程的部分结果
        double *all_results = NULL;
        if (mpi_rank == 0) {
            all_results = (double*)malloc(sizeof(double) * mpi_size);
        }
        
        double my_result = 3.14159 * (mpi_rank + 1);
        MPI_Gather(&my_result, 1, MPI_DOUBLE, all_results, 1, MPI_DOUBLE, 
                   0, MPI_COMM_WORLD);
        
        if (mpi_rank == 0) {
            printf("[MPI] Gathered results from all processes:\n");
            for (int i = 0; i < mpi_size; i++) {
                printf("  Process %d: %.6f\n", i, all_results[i]);
            }
            free(all_results);
        }
    } else {
        double my_result = 3.14159 * (mpi_rank + 1);
        MPI_Gather(&my_result, 1, MPI_DOUBLE, NULL, 1, MPI_DOUBLE, 
                   0, MPI_COMM_WORLD);
    }
    
    // 释放资源
    free(g_thread_data.data);
    free(g_thread_data.result);
    
    printf("[MPI] Process %d: Exiting\n", mpi_rank);
    
    // MPI 结束
    MPI_Finalize();
    
    return 0;
}
