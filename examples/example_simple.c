/**
 * 简化版 MPI + mythread 示例
 * 
 * 功能：
 * - 基本的 MPI 多进程 + 线程组并行
 * - 展示三层同步机制
 * - 进程间 MPI 通信
 * 
 * 编译：mpicc -o example_simple example_simple.c mythread/mythread.c -lpthread -lm
 * 运行：mpirun -np 2 ./example_simple
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <unistd.h>
#include "mythread/mythread.h"

#define N_GROUPS 4
#define N_THREADS_PER_GROUP 5

// 线程数据大小函数（必需）
int _gettdsize_() { return 64; }
int _getgdsize_() { return 64; }

int GetVInt(int volatile *volatile p){
    int tid = ti->ind;
    int gid = ti->igrp;
    // printf("Group=%d, Thread=%d: get func runed!\n", gid, tid);
    return *p;
}

// 全局计数器
static volatile int work_counter = 0;

/**
 * 工作线程函数
 */
void do_work() {
    int tid = ti->ind;
    int gid = ti->igrp;
    int s = 1;
    
    printf("  [Worker] MPI=%d, Group=%d, Thread=%d: waiting for work...\n", 
           mpi_id, gid, tid);
    
    // 等待组主线程信号 (SubThread wait GroupMain)
    sWaitGrp(s);

    // 执行工作
    work_counter++;
    printf("  [Worker] MPI=%d, Group=%d, Thread=%d: working (counter=%d)\n", 
           mpi_id, gid, tid, work_counter);
    sleep(10);
    // 模拟计算
    for (int i = 0; i < 1000000; i++) {
        work_counter += (i % 2 == 0) ? 1 : -1;
    }
    
    // 通知组主线程完成 (SubThread Set Group)
    sSetGrp(s);
    
    printf("  [Worker] MPI=%d, Group=%d, Thread=%d: work done\n", 
           mpi_id, gid, tid);
}

/**
 * 组主线程
 */
void group_main() {
    int gid = ti->igrp;
    int s = 1;
    printf("[GroupMain] MPI=%d, Group=%d: started\n", mpi_id, gid);
    
    // 等待主线程信号
    gWaitMain(s);
    
    printf("[GroupMain] MPI=%d, Group=%d: got signal from main, assigning work\n", 
           mpi_id, gid);
    
    // 通知子线程开始工作
    gSetSubs(s);

    gWaitSubs(s);
    // sleep(2);
    // 等待所有子线程完成
    
    printf("[GroupMain] MPI=%d, Group=%d: all workers done\n", mpi_id, gid);
    
    // 通知主线程本组完成
    gSetMain(s);

    s++;
    sleep(4);
    gWaitMain(s);
    gSetMain(s);
}

/**
 * 主线程（管理所有组）
 */
void main_thread() {
    int s = 1;
    printf("[Main] MPI=%d: main thread started\n", mpi_id);
    
    // 第一轮：通知所有组开始工作
    printf("[Main] MPI=%d: signaling all groups to start...\n", mpi_id);
    mSetGrps(s);  // Main Set Groups
    
    // 等待所有组完成
    mWaitGrps(s);  // Main Wait Groups
    printf("[Main] MPI=%d: all groups completed\n", mpi_id);
    
    s++;
    // 第二轮：重复
    printf("[Main] MPI=%d: starting second round...\n", mpi_id);
    // MSS;  // Main Set Subs (用于下一轮)
    mSetGrps(s);
    mWaitGrps(s);
    
    printf("[Main] MPI=%d: second round completed\n", mpi_id);
}

/**
 * 线程入口
 */
void thread_run() {
    int tid = gettid();
    
    printf("[INFO] thread info, gid=%d, tid=%d\n", ti->igrp, ti->ind);
    if (ti->igrp == -1) {
        // 主线程
        main_thread();
    } else if (ti->ind == 0){
        group_main();
    } else {
        // 工作线程
        do_work();
    }
}

int main(int argc, char *argv[]) {
    int rank, size;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    printf("[MPI] Process %d/%d starting...\n", rank+1, size);
    
    // 配置线程参数
    int NCorePClu = 16;
    int NCluPNode = 6;
    int NCorePGrp = 16;
    int NThPGrp = N_THREADS_PER_GROUP;  // 包含主线程
    int NGrpPProc = N_GROUPS;
    int NProcPNode = size;
    int ManageCoreId = -1;
    
    // 初始化线程
    int err = InitThreads(rank, NCorePClu, NCluPNode, NCorePGrp, 
                          NThPGrp, NGrpPProc, NProcPNode, &ManageCoreId);
    
    if (err != 0) {
        printf("[Error] Process %d: InitThreads failed: %d\n", rank, err);
        MPI_Finalize();
        return 1;
    }
    
    printf("[MPI] Process %d: Initialized %d groups, ThreadG=%d\n", 
           rank, N_GROUPS, ThreadG);
    
    // 启动线程
    StartThreads(thread_run);
    
    // 主线程也执行
    thread_run();
    
    // 结束线程
    EndThreads();
    
    printf("[MPI] Process %d: Threads ended\n", rank);
    
    // 进程间同步
    MPI_Barrier(MPI_COMM_WORLD);
    
    // 进程间通信示例
    int local_value = work_counter;
    int global_sum = 0;
    
    MPI_Reduce(&local_value, &global_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    
    if (rank == 0) {
        printf("\n[MPI] Global work counter sum: %d\n", global_sum);
    }
    
    MPI_Finalize();
    
    printf("[MPI] Process %d: Exiting\n", rank);
    return 0;
}
