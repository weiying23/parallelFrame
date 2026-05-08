/**
 * 波能量空间传播并行计算示例
 * 
 * 使用 mythread + MPI 实现 2D 波动方程的并行求解
 * 
 * 波动方程：∂²u/∂t² = c²(∂²u/∂x² + ∂²u/∂y²)
 * 
 * 编译：mpicc -Wall -O2 -o wave_propagation wave_propagation.c mythread/mythread.c -lpthread -lm
 * 运行：mpirun -np 4 ./wave_propagation
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h> 
#include <math.h>
#include <mpi.h>
#include "mythread/mythread.h"

/* ==================== 参数配置 ==================== */
#define NX          512     // X方向网格数
#define NY          512     // Y方向网格数
#define NT          1000    // 时间步数
#define DX          0.01    // 空间步长
#define DT          0.001   // 时间步长
#define C           1.0     // 波速
#define N_GROUPS    2       // 每进程线程组数
#define N_THREADS   3       // 每组工作线程数

//  Courant 条件：c*dt/dx <= 1/sqrt(2) 对于2D
#define COURANT     (C * DT / DX)

/* ==================== 数据结构 ==================== */

typedef struct {
    double *u_curr;     // 当前时刻波场
    double *u_prev;     // 前一时刻波场
    double *u_next;     // 下一时刻波场
    int x_start, x_end; // X方向负责范围
    int y_start, y_end; // Y方向负责范围
    int nx, ny;         // 局部网格尺寸
    double energy;      // 局部能量
} WaveField;

// 全局数据
typedef struct {
    WaveField *fields;  // 每个线程的波场数据
    double *u_global;   // 全局波场（仅在主线程）
    int total_nx, total_ny;
    double dx, dt, c;
    int n_steps;
    int current_step;
} SimulationData;

static SimulationData g_sim;

// /* ==================== 必需函数 ==================== */
// int _gettdsize_() {
//     return sizeof(WaveField);
// }

// int _getgdsize_() {
//     return sizeof(double) * 16; // 组共享缓冲区
// }

/* ==================== 工具函数 ==================== */

// 获取波场值（带边界检查）
static inline double get_u(double *u, int i, int j, int nx, int ny) {
    if (i < 0 || i >= nx || j < 0 || j >= ny) return 0.0;
    return u[i * ny + j];
}

// 设置波场值
static inline void set_u(double *u, int i, int j, int ny, double val) {
    u[i * ny + j] = val;
}

// 计算局部区域能量
static double compute_energy(WaveField *wf) {
    double energy = 0.0;
    for (int i = wf->x_start; i < wf->x_end; i++) {
        for (int j = wf->y_start; j < wf->y_end; j++) {
            double u = wf->u_curr[i * wf->ny + j];
            // 能量密度 = u² + (du/dt)² (简化)
            energy += u * u;
        }
    }
    return energy * DX * DX;
}

// 初始化高斯波包
static void init_gaussian_pulse(WaveField *wf, double cx, double cy, double sigma) {
    for (int i = wf->x_start; i < wf->x_end; i++) {
        for (int j = wf->y_start; j < wf->y_end; j++) {
            double x = i * DX;
            double y = j * DX;
            double dist2 = (x - cx) * (x - cx) + (y - cy) * (y - cy);
            double val = exp(-dist2 / (2 * sigma * sigma));
            set_u(wf->u_curr, i, j, wf->ny, val);
            set_u(wf->u_prev, i, j, wf->ny, val); // 初始速度为0
        }
    }
}

/* ==================== 核心计算 ==================== */

// 更新波场（有限差分）
static void update_wavefield(WaveField *wf) {
    double courant2 = COURANT * COURANT;
    
    for (int i = wf->x_start; i < wf->x_end; i++) {
        for (int j = wf->y_start; j < wf->y_end; j++) {
            // 获取周围点（处理边界）
            double u_ij = wf->u_curr[i * wf->ny + j];
            double u_im1 = get_u(wf->u_curr, i-1, j, wf->nx, wf->ny);
            double u_ip1 = get_u(wf->u_curr, i+1, j, wf->nx, wf->ny);
            double u_jm1 = get_u(wf->u_curr, i, j-1, wf->nx, wf->ny);
            double u_jp1 = get_u(wf->u_curr, i, j+1, wf->nx, wf->ny);
            
            // 有限差分公式
            double laplacian = (u_ip1 + u_im1 + u_jp1 + u_jm1 - 4 * u_ij);
            double u_next = 2 * u_ij - wf->u_prev[i * wf->ny + j] + courant2 * laplacian;
            
            set_u(wf->u_next, i, j, wf->ny, u_next);
        }
    }
}

// 交换时间层
static void swap_buffers(WaveField *wf) {
    double *temp = wf->u_prev;
    wf->u_prev = wf->u_curr;
    wf->u_curr = wf->u_next;
    wf->u_next = temp;
}

/* ==================== 线程函数 ==================== */

// 工作线程：执行实际计算
static void worker_thread() {
    int tid = gettid();
    int gid = ti->igrp;
    WaveField *wf = (WaveField*)ti->td;
    
    printf("[Worker] MPI=%d, Group=%d, Thread=%d ready (region: x=[%d,%d), y=[%d,%d))\n",
           mpi_id, gid, tid, wf->x_start, wf->x_end, wf->y_start, wf->y_end);
    
    for (int step = 0; step < g_sim.n_steps; step++) {
        // 等待组主线程信号
        SWG;
        
        // 更新波场
        update_wavefield(wf);
        
        // 计算局部能量
        wf->energy = compute_energy(wf);
        
        // 通知组主线程完成
        SSG;
    }
    
    printf("[Worker] MPI=%d, Group=%d, Thread=%d finished\n", mpi_id, gid, tid);
}

// 组主线程：管理组内工作线程
static void group_main_thread() {
    int gid = ti->igrp;
    WaveField *wf = (WaveField*)ti->td;
    
    printf("[GroupMain] MPI=%d, Group=%d started\n", mpi_id, gid);
    
    // 初始化波场（高斯波包，不同组不同位置）
    double cx = 0.25 + gid * 0.25;  // 不同组的波源位置
    double cy = 0.5;
    init_gaussian_pulse(wf, cx, cy, 0.05);
    
    for (int step = 0; step < g_sim.n_steps; step++) {
        // 等待主线程信号
        GWM;
        
        // 通知工作线程开始计算
        GSS;
        
        // 组主也参与计算（最后一行，避免边界问题）
        update_wavefield(wf);
        
        // 等待所有工作线程完成
        GWS;
        
        // 交换时间层
        swap_buffers(wf);
        
        // 计算组总能量
        double group_energy = wf->energy;
        // 收集工作线程的能量（简化，实际应该汇总）
        
        // 通知主线程本组完成
        GSM;
    }
    
    printf("[GroupMain] MPI=%d, Group=%d finished\n", mpi_id, gid);
}

// 主线程：协调所有组，处理MPI通信
static void main_thread() {
    printf("[Main] MPI=%d: Main thread started, managing %d groups\n", 
           mpi_id, md.ngrp);
    
    double start_time = MPI_Wtime();
    
    for (int step = 0; step < g_sim.n_steps; step++) {
        g_sim.current_step = step;
        
        // 每100步输出一次
        if (step % 100 == 0 && mpi_id == 0) {
            printf("[Main] Step %d/%d...\n", step, g_sim.n_steps);
        }
        
        // 通知所有组开始计算
        MSG;

        // GSS;

        // GWS;
        
        // 等待所有组完成
        MWG;
        
        // 进程内能量归约（简化版）
        double local_energy = 0.0;
        // 实际应该收集各组能量
        
        // 进程间能量归约
        double global_energy = 0.0;
        MPI_Allreduce(&local_energy, &global_energy, 1, MPI_DOUBLE, 
                      MPI_SUM, MPI_COMM_WORLD);
        
        // 准备下一时间步
        MSS;
    }
    
    double end_time = MPI_Wtime();
    
    if (mpi_id == 0) {
        printf("\n[Main] Simulation completed in %.3f seconds\n", end_time - start_time);
        printf("[Main] Performance: %.2f Mgrid/s\n", 
               (double)NX * NY * g_sim.n_steps / (end_time - start_time) / 1e6);
    }
}

// 线程入口
void thread_run() {
    int tid = gettid();
    int gid = ti->igrp;
    printf("[INFO] MPI=%d, Group=%d, Thread=%d finished\n", mpi_id, gid, tid);

    // sleep(30);
    
    if (gid == -1) {
        main_thread();
    } else if (tid==1){
        group_main_thread();
    }
    else {
        worker_thread();
    }
}

/* ==================== 初始化 ==================== */

static void init_simulation(int mpi_rank, int mpi_size) {
    g_sim.total_nx = NX;
    g_sim.total_ny = NY;
    g_sim.dx = DX;
    g_sim.dt = DT;
    g_sim.c = C;
    g_sim.n_steps = NT;
    g_sim.current_step = 0;
    
    // 分配全局波场（仅在主线程使用）
    if (mpi_rank == 0) {
        g_sim.u_global = (double*)calloc(NX * NY, sizeof(double));
    }
}

// 分配波场内存
static void alloc_wavefield(WaveField *wf, int nx, int ny) {
    wf->nx = nx;
    wf->ny = ny;
    wf->u_curr = (double*)calloc(nx * ny, sizeof(double));
    wf->u_prev = (double*)calloc(nx * ny, sizeof(double));
    wf->u_next = (double*)calloc(nx * ny, sizeof(double));
    wf->energy = 0.0;
}

// 释放波场内存
static void free_wavefield(WaveField *wf) {
    free(wf->u_curr);
    free(wf->u_prev);
    free(wf->u_next);
}

// 设置线程数据区域（在 InitThreads 之后，StartThreads 之前）
static void setup_thread_data() {
    int threads_per_group = N_THREADS;
    
    // Y方向划分给各组
    int ny_per_group = NY / md.ngrp;
    
    // X方向划分给组内线程
    int nx_per_thread = NX / threads_per_group;
    
    // 为每个线程设置数据区域
    for (int g = 0; g < md.ngrp; g++) {
        threadGroup *pg = md.grps[g];
        
        for (int t = 0; t < pg->Nthreads; t++) {
            THREADINFO *pti = &pg->threads[t];
            WaveField *wf = (WaveField*)pti->td;
            
            // Y方向：按组分
            wf->y_start = g * ny_per_group;
            wf->y_end = (g == md.ngrp - 1) ? NY : (g + 1) * ny_per_group;
            
            // X方向：按线程分（包含边界缓冲区）
            wf->x_start = t * nx_per_thread;
            wf->x_end = (t == pg->Nthreads - 1) ? NX : (t + 1) * nx_per_thread;
            
            // 分配内存（包含边界）
            alloc_wavefield(wf, NX, wf->y_end - wf->y_start);
            
            printf("[Init] Group=%d, Thread=%d: region [%d,%d) x [%d,%d)\n",
                   g, t, wf->x_start, wf->x_end, wf->y_start, wf->y_end);
        }
    }
}

/* ==================== 主函数 ==================== */

int main(int argc, char *argv[]) {
    int mpi_rank, mpi_size;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    
    if (mpi_rank == 0) {
        printf("============================================\n");
        printf("    2D Wave Propagation Simulation\n");
        printf("============================================\n");
        printf("Grid: %d x %d\n", NX, NY);
        printf("Steps: %d\n", NT);
        printf("MPI Processes: %d\n", mpi_size);
        printf("Thread Groups per Process: %d\n", N_GROUPS);
        printf("Threads per Group: %d\n", N_THREADS + 1);
        printf("Courant Number: %.4f\n", COURANT);
        printf("============================================\n\n");
    }
    
    // 初始化模拟参数
    init_simulation(mpi_rank, mpi_size);
    
    // 线程配置
    int NCorePClu = 16;
    int NCluPNode = 2;
    int NCorePGrp = 16;
    int NThPGrp = N_THREADS ;
    int NGrpPProc = N_GROUPS;
    int NProcPNode = mpi_size;
    int ManageCoreId = -2;
    
    // 初始化线程系统
    int err = InitThreads(mpi_rank, NCorePClu, NCluPNode, NCorePGrp, 
                          NThPGrp, NGrpPProc, NProcPNode, &ManageCoreId);
    
    if (err != 0) {
        printf("[Error] MPI=%d: InitThreads failed: %d\n", mpi_rank, err);
        MPI_Finalize();
        return 1;
    }
    
    // 设置线程数据区域
    setup_thread_data();
    
    // 启动线程
    StartThreads(thread_run);
    
    // 主线程也参与模拟
    thread_run();
    
    // 结束线程
    EndThreads();
    
    // 清理
    if (mpi_rank == 0) {
        free(g_sim.u_global);
    }
    
    MPI_Finalize();
    return 0;
}
