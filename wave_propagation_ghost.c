/**
 * 波能量空间传播并行计算示例 - 带幽灵层版本
 * 
 * 使用 mythread + MPI 实现 2D 波动方程的并行求解
 * 支持幽灵层（Ghost Cells）边界交换，包括：
 * - Y方向：MPI进程间 + 组间交换
 * - X方向：组内线程间交换
 * 
 * 波动方程：∂²u/∂t² = c²(∂²u/∂x² + ∂²u/∂y²)
 * 
 * 编译：mpicc -Wall -O2 -o wave_propagation_ghost wave_propagation_ghost.c mythread/mythread.c -lpthread -lm
 * 运行：mpirun -np 4 ./wave_propagation_ghost
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
#define NT          400    // 时间步数
#define DX          0.01    // 空间步长
#define DT          0.001   // 时间步长
#define C           1.0     // 波速
#define N_GROUPS    2       // 每进程线程组数
#define N_THREADS   3       // 每组工作线程数（不含主线程）
#define GHOST_WIDTH 1       // 幽灵层宽度（有限差分需要1层）

//  Courant 条件：c*dt/dx <= 1/sqrt(2) 对于2D
#define COURANT     (C * DT / DX)

/* ==================== 数据结构 ==================== */

typedef struct {
    double *u_curr;     // 当前时刻波场（包含幽灵层）
    double *u_prev;     // 前一时刻波场（包含幽灵层）
    double *u_next;     // 下一时刻波场（包含幽灵层）
    
    // 实际计算区域（不含幽灵层）
    int x_start, x_end; // X方向负责范围（全局坐标）
    int y_start, y_end; // Y方向负责范围（全局坐标）
    int nx, ny;         // 局部网格尺寸（含幽灵层）
    int nx_inner, ny_inner; // 内部网格尺寸（不含幽灵层）
    double energy;      // 局部能量
    
    // 邻居信息
    int neighbor_up;    // 上方邻居MPI rank (-1表示无)
    int neighbor_down;  // 下方邻居MPI rank (-1表示无)
    int igrp;           // 所属组ID
    int itid;           // 组内线程ID
} WaveField;

// 组共享边界数据（用于组内线程间交换）
typedef struct {
    // 左右边界缓冲区：每个线程需要左右各一个缓冲区
    // 格式：buffer[线程ID][边][数据]
    double *left_send;   // 发送给左邻居的数据
    double *left_recv;   // 从左邻居接收的数据
    double *right_send;  // 发送给右邻居的数据
    double *right_recv;  // 从右邻居接收的数据
    int buf_size;        // 缓冲区大小（ny_inner）
    
    // 同步标记
    volatile int *thread_done;    // 各线程完成标记
    volatile int *boundary_ready; // 边界准备好标记
} GroupBoundaryData;

// 全局数据
typedef struct {
    WaveField *fields;  // 每个线程的波场数据
    double *u_global;   // 全局波场（仅在主线程）
    int total_nx, total_ny;
    double dx, dt, c;
    int n_steps;
    int current_step;
    int mpi_size;       // MPI进程数
} SimulationData;

static SimulationData g_sim;

// MPI边界交换缓冲区
static double *send_buf_up, *send_buf_down;
static double *recv_buf_up, *recv_buf_down;

int GetVInt(int volatile *volatile p){
    int tid = ti->ind;
    int gid = ti->igrp;
    // printf("Group=%d, Thread=%d: get func runed!\n", gid, tid);
    return *p;
}
/* ==================== 必需函数 ==================== */
int _gettdsize_() {
    return sizeof(WaveField);
}

int _getgdsize_() {
    // 组共享数据：边界缓冲区 + 同步标记
    return sizeof(GroupBoundaryData) + sizeof(double) * NY * 4 + sizeof(int) * (N_THREADS + 1) * 2;
}

/* ==================== 索引工具函数 ==================== */

// 局部索引转线性偏移
static inline int idx(int i, int j, int ny) {
    return i * ny + j;
}

/* ==================== 核心计算 ==================== */

// 计算局部区域能量（仅内部区域）
static double compute_energy(WaveField *wf) {
    double energy = 0.0;
    for (int i = 0; i < wf->nx_inner; i++) {
        for (int j = 0; j < wf->ny_inner; j++) {
            int li = i + GHOST_WIDTH;
            int lj = j + GHOST_WIDTH;
            double u = wf->u_curr[idx(li, lj, wf->ny)];
            energy += u * u;
        }
    }
    return energy * DX * DX;
}

// 初始化高斯波包
static void init_gaussian_pulse(WaveField *wf, double cx, double cy, double sigma) {
    for (int i = 0; i < wf->nx_inner; i++) {
        for (int j = 0; j < wf->ny_inner; j++) {
            int gi = wf->x_start + i;  // 全局坐标
            int gj = wf->y_start + j;
            double x = gi * DX;
            double y = gj * DX;
            double dist2 = (x - cx) * (x - cx) + (y - cy) * (y - cy);
            double val = exp(-dist2 / (2 * sigma * sigma));
            
            int li = i + GHOST_WIDTH;
            int lj = j + GHOST_WIDTH;
            wf->u_curr[idx(li, lj, wf->ny)] = val;
            wf->u_prev[idx(li, lj, wf->ny)] = val;
        }
    }
}

// 更新波场（有限差分）- 使用幽灵层
static void update_wavefield(WaveField *wf) {
    double courant2 = COURANT * COURANT;
    
    // 仅更新内部区域
    for (int i = 0; i < wf->nx_inner; i++) {
        for (int j = 0; j < wf->ny_inner; j++) {
            int li = i + GHOST_WIDTH;
            int lj = j + GHOST_WIDTH;
            
            // 从幽灵层读取边界数据
            double u_ij = wf->u_curr[idx(li, lj, wf->ny)];
            double u_im1 = wf->u_curr[idx(li-1, lj, wf->ny)];
            double u_ip1 = wf->u_curr[idx(li+1, lj, wf->ny)];
            double u_jm1 = wf->u_curr[idx(li, lj-1, wf->ny)];
            double u_jp1 = wf->u_curr[idx(li, lj+1, wf->ny)];
            
            double laplacian = (u_ip1 + u_im1 + u_jp1 + u_jm1 - 4 * u_ij);
            double u_next_val = 2 * u_ij - wf->u_prev[idx(li, lj, wf->ny)] + courant2 * laplacian;
            
            wf->u_next[idx(li, lj, wf->ny)] = u_next_val;
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
 
/* ==================== Y方向边界交换（MPI）==================== */

// 提取Y方向边界（上下边界是X方向的整条线）
static void pack_y_boundary(WaveField *wf, double *buf, int is_top) {
    // is_top: 1表示上边界(最大i), 0表示下边界(最小i)
    int i = is_top ? wf->nx_inner : 1;  // 内部区域的边界行
    for (int j = 0; j < wf->ny_inner; j++) {
        int lj = j + GHOST_WIDTH;
        buf[j] = wf->u_curr[idx(i, lj, wf->ny)];
    }
}

// 解包Y方向边界到幽灵层
static void unpack_y_boundary(WaveField *wf, double *buf, int is_top) {
    int i = is_top ? wf->nx_inner + GHOST_WIDTH : 0;  // 幽灵层位置
    for (int j = 0; j < wf->ny_inner; j++) {
        int lj = j + GHOST_WIDTH;
        wf->u_curr[idx(i, lj, wf->ny)] = buf[j];
    }
}

// MPI进程间Y方向边界交换
static void mpi_exchange_y_boundaries(WaveField *wf) {
    MPI_Status status;
    
    // 与上方进程交换（发送上边界，接收上幽灵层）
    if (wf->neighbor_up >= 0) {
        pack_y_boundary(wf, send_buf_up, 1);
        MPI_Sendrecv(send_buf_up, wf->ny_inner, MPI_DOUBLE, wf->neighbor_up, 0,
                     recv_buf_up, wf->ny_inner, MPI_DOUBLE, wf->neighbor_up, 1,
                     MPI_COMM_WORLD, &status);
        unpack_y_boundary(wf, recv_buf_up, 1);
    }
    
    // 与下方进程交换（发送下边界，接收下幽灵层）
    if (wf->neighbor_down >= 0) {
        pack_y_boundary(wf, send_buf_down, 0);
        MPI_Sendrecv(send_buf_down, wf->ny_inner, MPI_DOUBLE, wf->neighbor_down, 1,
                     recv_buf_down, wf->ny_inner, MPI_DOUBLE, wf->neighbor_down, 0,
                     MPI_COMM_WORLD, &status);
        unpack_y_boundary(wf, recv_buf_down, 0);
    }
}

// 同一MPI进程内组间边界交换（共享内存）
static void exchange_boundaries_between_groups() {
    // 遍历相邻的组对
    for (int g = 0; g < md.ngrp - 1; g++) {
        threadGroup *g0 = md.grps[g];      // 当前组（下方）
        threadGroup *g1 = md.grps[g+1];    // 下一组（上方）
        
        WaveField *wf0 = (WaveField*)g0->threads[0].td;  // 组0的组主
        WaveField *wf1 = (WaveField*)g1->threads[0].td;  // 组1的组主
        
        int ny_inner = wf0->ny_inner;
        
        // 组0的下边界（内部最后一行）→ 组1的下幽灵层
        // 组1的上边界（内部第一行）→ 组0的上幽灵层
        for (int j = 0; j < ny_inner; j++) {
            int lj = j + GHOST_WIDTH;
            
            // 组0的底行（i = nx_inner，内部最后一行）
            double bottom_val = wf0->u_curr[idx(wf0->nx_inner, lj, wf0->ny)];
            // 复制到组1的下幽灵层（i = 0）
            wf1->u_curr[idx(0, lj, wf1->ny)] = bottom_val;
            
            // 组1的顶行（i = GHOST_WIDTH，内部第一行）
            double top_val = wf1->u_curr[idx(GHOST_WIDTH, lj, wf1->ny)];
            // 复制到组0的上幽灵层（i = nx_inner + GHOST_WIDTH）
            wf0->u_curr[idx(wf0->nx_inner + GHOST_WIDTH, lj, wf0->ny)] = top_val;
        }
    }
}

/* ==================== X方向边界交换（线程间）==================== */

// 获取组的边界数据结构
static GroupBoundaryData* get_group_boundary_data(int gid) {
    threadGroup *pg = md.grps[gid];
    char *gdata = (char*)pg->gd;  // gd 是 threadGroup 的成员
    return (GroupBoundaryData*)gdata;
}

// 初始化组边界数据结构
static void init_group_boundary_data(int gid) {
    threadGroup *pg = md.grps[gid];
    char *gdata_raw = (char*)pg->gd;  // gd 是 threadGroup 的成员
    GroupBoundaryData *gbd = (GroupBoundaryData*)gdata_raw;
    
    WaveField *wf0 = (WaveField*)pg->threads[0].td;
    int buf_size = wf0->ny_inner;
    gbd->buf_size = buf_size;
    
    // 分配缓冲区（在gd内存区域内）
    char *ptr = gdata_raw + sizeof(GroupBoundaryData);
    gbd->left_send = (double*)ptr;  ptr += sizeof(double) * buf_size;
    gbd->left_recv = (double*)ptr;  ptr += sizeof(double) * buf_size;
    gbd->right_send = (double*)ptr; ptr += sizeof(double) * buf_size;
    gbd->right_recv = (double*)ptr; ptr += sizeof(double) * buf_size;
    
    gbd->thread_done = (volatile int*)ptr; ptr += sizeof(int) * (N_THREADS + 1);
    gbd->boundary_ready = (volatile int*)ptr;
    
    // 初始化标记
    for (int i = 0; i <= N_THREADS; i++) {
        gbd->thread_done[i] = 0;
        gbd->boundary_ready[i] = 0;
    }
}

// 提取X方向边界（左右边界是Y方向的整条线）
static void pack_x_boundary(WaveField *wf, double *buf, int is_right) {
    // is_right: 1表示右边界(最大j), 0表示左边界(最小j)
    int j = is_right ? wf->ny_inner : 1;  // 内部区域的边界列
    for (int i = 0; i < wf->nx_inner; i++) {
        int li = i + GHOST_WIDTH;
        buf[i] = wf->u_curr[idx(li, j, wf->ny)];
    }
}

// 解包X方向边界到幽灵层
static void unpack_x_boundary(WaveField *wf, double *buf, int is_right) {
    int j = is_right ? wf->ny_inner + GHOST_WIDTH : 0;  // 幽灵层位置
    for (int i = 0; i < wf->nx_inner; i++) {
        int li = i + GHOST_WIDTH;
        wf->u_curr[idx(li, j, wf->ny)] = buf[i];
    }
}

// 组内X方向边界交换（由组主线程协调）
static void exchange_x_boundaries_group(int gid) {
    threadGroup *pg = md.grps[gid];
    GroupBoundaryData *gbd = get_group_boundary_data(gid);
    int nthreads = pg->Nthreads;
    
    // 等待所有线程完成计算
    // for (int t = 1; t < nthreads; t++) {
    //     while (!gbd->thread_done[t]) {
    //         // 自旋等待（可用更高效的同步）
    //     }
    // }
    
    // 收集各线程的边界
    for (int t = 0; t < nthreads; t++) {
        THREADINFO *pti = &pg->threads[t];
        WaveField *wf = (WaveField*)pti->td;
        
        // 提取左边界
        if (t > 0) {  // 不是最左线程
            pack_x_boundary(wf, &gbd->left_send[t * gbd->buf_size], 0);
        }
        // 提取右边界
        if (t < nthreads - 1) {  // 不是最右线程
            pack_x_boundary(wf, &gbd->right_send[t * gbd->buf_size], 1);
        }
    }
    
    // 交换边界：左线程的右边界 → 右线程的左幽灵层
    for (int t = 0; t < nthreads - 1; t++) {
        double *src = &gbd->right_send[t * gbd->buf_size];      // 线程t的右边界
        double *dst = &gbd->left_recv[(t+1) * gbd->buf_size];   // 线程t+1的左接收
        memcpy(dst, src, sizeof(double) * gbd->buf_size);
        
        src = &gbd->left_send[(t+1) * gbd->buf_size];           // 线程t+1的左边界
        dst = &gbd->right_recv[t * gbd->buf_size];              // 线程t的右接收
        memcpy(dst, src, sizeof(double) * gbd->buf_size);
    }
    
    // 分发边界到各线程的幽灵层
    for (int t = 0; t < nthreads; t++) {
        THREADINFO *pti = &pg->threads[t];
        WaveField *wf = (WaveField*)pti->td;
        
        if (t > 0) {  // 有左邻居
            unpack_x_boundary(wf, &gbd->left_recv[t * gbd->buf_size], 0);
        }
        if (t < nthreads - 1) {  // 有右邻居
            unpack_x_boundary(wf, &gbd->right_recv[t * gbd->buf_size], 1);
        }
    }
    
    // 重置标记，通知线程可以继续
    for (int t = 1; t < nthreads; t++) {
        gbd->thread_done[t] = 0;
        gbd->boundary_ready[t] = 1;
    }
}

/* ==================== 线程函数 ==================== */

// 工作线程
static void worker_thread() {
    int tid = gettid();
    int gid = ti->igrp;
    int s = 1;
    WaveField *wf = (WaveField*)ti->td;
    GroupBoundaryData *gbd = get_group_boundary_data(gid);
    
    printf("[Worker] MPI=%d, Group=%d, Thread=%d ready (region: x=[%d,%d), y=[%d,%d))\n",
           mpi_id, gid, tid, wf->x_start, wf->x_end, wf->y_start, wf->y_end);
    
    for (int step = 0; step < g_sim.n_steps; step++) {
        // 等待组主信号（边界已交换好）
        s = step + 2;
        sWaitGrp(s);
        // printf("[Worker] MPI=%d, Group=%d, Thread=%d SWG\n", mpi_id, gid, tid);
        
        // 更新波场
        update_wavefield(wf);
        
        // 计算局部能量
        wf->energy = compute_energy(wf);
        
        // 标记计算完成
        // gbd->thread_done[tid] = 1;
        
        // 等待边界准备好
        // while (!gbd->boundary_ready[tid]) {}
        // gbd->boundary_ready[tid] = 0;
        
        // 通知组主线程完成
        sSetGrp(s);
        // printf("[Worker] MPI=%d, Group=%d, Thread=%d SSG\n", mpi_id, gid, tid);
    }
    
    printf("[Worker] MPI=%d, Group=%d, Thread=%d finished\n", mpi_id, gid, tid);
}

// 组主线程
static void group_main_thread() {
    int gid = ti->igrp;
    int tid = gettid();
    WaveField *wf = (WaveField*)ti->td;
    
    // printf("[GroupMain] MPI=%d, Group=%d, %d started\n", mpi_id, gid, tid);
    
    // 初始化组边界数据
    init_group_boundary_data(gid);
    
    // 初始化波场
    double cx = 0.25 + gid * 0.25;
    double cy = 0.5;
    int s;
    init_gaussian_pulse(wf, cx, cy, 0.05);
    
    for (int step = 0; step < g_sim.n_steps; step++) {
        // 等待主线程信号（Y方向边界已交换）
        s = step + 2;
        gWaitMain(s);
        
        // printf("[GroupMain]MPI=%d, Group=%d, %d GWM\n",mpi_id, gid, tid);
        // 通知工作线程开始
        gSetSubs(s);
        // printf("[GroupMain]MPI=%d, Group=%d, %d GSS\n",mpi_id, gid, tid);
        // 组主也参与计算
        update_wavefield(wf);
        
        // 标记完成并等待工作线程
        GroupBoundaryData *gbd = get_group_boundary_data(gid);
        gbd->thread_done[0] = 1;
        
        gWaitSubs(s);
        // sleep(1);
        // 协调X方向边界交换
        exchange_x_boundaries_group(gid);
        
        // 交换时间层
        swap_buffers(wf);
        
        // 重置组主标记
        gbd->thread_done[0] = 0;
        // SSG;

        // printf("[GroupMain]MPI=%d, Group=%d, %d GWS\n", mpi_id, gid, tid);
        // 通知主线程本组完成
        gSetMain(s);
        // printf("[GroupMain]MPI=%d, Group=%d, %d GSM\n", mpi_id, gid, tid);
    }
    
    // printf("[GroupMain] MPI=%d, Group=%d finished\n", mpi_id, gid);
}

// 主线程
static void main_thread() {
    printf("[Main] MPI=%d: Main thread started, managing %d groups\n", 
           mpi_id, md.ngrp);
    
    double start_time = MPI_Wtime();
    double initial_energy = -1.0;
    int s;
    
    for (int step = 0; step < g_sim.n_steps; step++) {
        s = step+2;
        g_sim.current_step = step;
        
        if (step % 100 == 0 && mpi_id == 0) {
            printf("[Main] Step %d/%d...\n", step, g_sim.n_steps);
        }
        
        // Y方向：MPI边界交换（为每个组，跨进程）
        for (int g = 0; g < md.ngrp; g++) {
            threadGroup *pg = md.grps[g];
            WaveField *wf = (WaveField*)pg->threads[0].td;
            mpi_exchange_y_boundaries(wf);
        }
        
        // Y方向：进程内组间边界交换（共享内存）
        exchange_boundaries_between_groups();
        
        // 通知所有组开始
        mSetGrps(s);
        // printf("[MAIN] MSG\n");
        
        // 等待所有组完成
        mWaitGrps(s);
        // printf("[MAIN] MWG\n");
        
        // 能量统计
        double local_energy = 0.0;
        for (int g = 0; g < md.ngrp; g++) {
            threadGroup *pg = md.grps[g];
            for (int t = 0; t < pg->Nthreads; t++) {
                WaveField *wf = (WaveField*)pg->threads[t].td;
                local_energy += wf->energy;
            }
        }
        
        double global_energy = 0.0;
        MPI_Allreduce(&local_energy, &global_energy, 1, MPI_DOUBLE, 
                      MPI_SUM, MPI_COMM_WORLD);
        
        if (step == 0) {
            initial_energy = global_energy;
            if (mpi_id == 0) {
                printf("[Main] Initial energy: %.6f\n", initial_energy);
            }
        } else if (step % 100 == 0 && mpi_id == 0) {
            double rel_diff = fabs(global_energy - initial_energy) / initial_energy;
            printf("[Main] Step %d: Energy = %.6f, Rel. diff = %.2e\n", 
                   step, global_energy, rel_diff);
        }
        
        // 准备下一时间步
        // MSS;
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
    int tid = ti->ind;
    int gid = ti->igrp;
    int nt = getnt();
    printf("[Info] MPI=%d, Group=%d, Thread=%d of %d start run!\n",
           mpi_id, gid, tid, getnt());
    if (gid == -1 ) {
        main_thread();
    } else if (tid == 0) {
        group_main_thread();
    } else {
        worker_thread();
    }
    printf("[Info] MPI=%d, Group=%d, Thread=%d finish!\n",
        mpi_id, gid, tid);
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
    g_sim.mpi_size = mpi_size;
    
    if (mpi_rank == 0) {
        g_sim.u_global = (double*)calloc(NX * NY, sizeof(double));
    }
    
    // 分配MPI边界缓冲区
    send_buf_up = (double*)malloc(sizeof(double) * NY);
    send_buf_down = (double*)malloc(sizeof(double) * NY);
    recv_buf_up = (double*)malloc(sizeof(double) * NY);
    recv_buf_down = (double*)malloc(sizeof(double) * NY);
}

// 分配波场内存（含幽灵层）
static void alloc_wavefield(WaveField *wf, int nx_inner, int ny_inner) {
    wf->nx_inner = nx_inner;
    wf->ny_inner = ny_inner;
    wf->nx = nx_inner + 2 * GHOST_WIDTH;
    wf->ny = ny_inner + 2 * GHOST_WIDTH;
    
    int total_size = wf->nx * wf->ny;
    wf->u_curr = (double*)calloc(total_size, sizeof(double));
    wf->u_prev = (double*)calloc(total_size, sizeof(double));
    wf->u_next = (double*)calloc(total_size, sizeof(double));
    wf->energy = 0.0;
}

// 释放波场内存
static void free_wavefield(WaveField *wf) {
    free(wf->u_curr);
    free(wf->u_prev);
    free(wf->u_next);
}

// 计算MPI邻居
static void setup_mpi_neighbors(int mpi_rank, int mpi_size, WaveField *wf) {
    int ny_per_proc = NY / mpi_size;
    int y_start_proc = mpi_rank * ny_per_proc;
    int y_end_proc = (mpi_rank == mpi_size - 1) ? NY : (mpi_rank + 1) * ny_per_proc;
    
    wf->neighbor_up = (y_end_proc < NY) ? mpi_rank + 1 : -1;
    wf->neighbor_down = (y_start_proc > 0) ? mpi_rank - 1 : -1;
}

// 设置线程数据区域
static void setup_thread_data(int mpi_rank, int mpi_size) {
    int nthreads = N_THREADS + 1;  // 包含主线程
    
    // 本进程负责的Y范围
    int ny_per_proc = NY / mpi_size;
    int y_start_proc = mpi_rank * ny_per_proc;
    int y_end_proc = (mpi_rank == mpi_size - 1) ? NY : (mpi_rank + 1) * ny_per_proc;
    int ny_proc = y_end_proc - y_start_proc;
    
    // Y方向划分给各组
    int ny_per_group = ny_proc / md.ngrp;
    
    // X方向划分给组内线程
    int nx_per_thread = NX / nthreads;
    
    for (int g = 0; g < md.ngrp; g++) {
        threadGroup *pg = md.grps[g];
        
        int group_y_start = y_start_proc + g * ny_per_group;
        int group_y_end = (g == md.ngrp - 1) ? y_end_proc : y_start_proc + (g + 1) * ny_per_group;
        int group_ny = group_y_end - group_y_start;
        
        for (int t = 0; t < pg->Nthreads; t++) {
            THREADINFO *pti = &pg->threads[t];
            WaveField *wf = (WaveField*)pti->td;
            
            int thread_x_start = t * nx_per_thread;
            int thread_x_end = (t == pg->Nthreads - 1) ? NX : (t + 1) * nx_per_thread;
            int thread_nx = thread_x_end - thread_x_start;
            
            wf->x_start = thread_x_start;
            wf->x_end = thread_x_end;
            wf->y_start = group_y_start;
            wf->y_end = group_y_end;
            wf->igrp = g;
            wf->itid = t;
            
            alloc_wavefield(wf, thread_nx, group_ny);
            
            if (t == 0) {
                setup_mpi_neighbors(mpi_rank, mpi_size, wf);
            }
            
            printf("[Init] MPI=%d, G=%d, T=%d: region [%d,%d)x[%d,%d), buf=%dx%d\n",
                   mpi_id, g, t, wf->x_start, wf->x_end, wf->y_start, wf->y_end,
                   wf->nx, wf->ny);
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
        printf("    2D Wave Propagation with Ghost Cells\n");
        printf("    X/Y Direction Boundary Exchange\n");
        printf("============================================\n");
        printf("Grid: %d x %d\n", NX, NY);
        printf("Steps: %d\n", NT);
        printf("Ghost Width: %d\n", GHOST_WIDTH);
        printf("MPI Processes: %d\n", mpi_size);
        printf("Thread Groups per Process: %d\n", N_GROUPS);
        printf("Threads per Group: %d\n", N_THREADS + 1);
        printf("Courant Number: %.4f\n", COURANT);
        printf("============================================\n\n");
    }
    
    init_simulation(mpi_rank, mpi_size);
    
    int NCorePClu = 16;
    int NCluPNode = 2;
    int NCorePGrp = 16;
    int NThPGrp = N_THREADS + 1;
    int NGrpPProc = N_GROUPS;
    int NProcPNode = mpi_size;
    int ManageCoreId = 0;
    
    int err = InitThreads(mpi_rank, NCorePClu, NCluPNode, NCorePGrp, 
                          NThPGrp, NGrpPProc, NProcPNode, &ManageCoreId);
    
    if (err != 0) {
        printf("[Error] MPI=%d: InitThreads failed: %d\n", mpi_rank, err);
        MPI_Finalize();
        return 1;
    }
    
    setup_thread_data(mpi_rank, mpi_size);
    
    StartThreads(thread_run);
    thread_run();
    EndThreads();
    
    // 清理
    if (mpi_rank == 0) {
        free(g_sim.u_global);
    }
    free(send_buf_up);
    free(send_buf_down);
    free(recv_buf_up);
    free(recv_buf_down);
    
    MPI_Finalize();
    return 0;
}
