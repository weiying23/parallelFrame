/**
 * 波能量传播 - 可视化版本
 * 
 * 输出波场快照到文件，可用于后续可视化
 * 
 * 编译：mpicc -Wall -O2 -o wave_visual wave_visual.c mythread/mythread.c -lpthread -lm
 * 运行：mpirun -np 2 ./wave_visual
 * 可视化：python3 visualize_wave.py
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "mythread/mythread.h"

/* ==================== 参数配置 ==================== */
#define NX          256     // X方向网格数
#define NY          256     // Y方向网格数
#define NT          500     // 时间步数
#define SAVE_EVERY  50      // 每多少步保存一次
#define DX          0.01
#define DT          0.001
#define C           2.0

#define N_GROUPS    1
#define N_THREADS   3

#define COURANT     (C * DT / DX)

/* ==================== 数据结构 ==================== */

typedef struct {
    double *u_curr;
    double *u_prev;
    double *u_next;
    int x_start, x_end;
    int y_start, y_end;
    int nx, ny;
    double max_amplitude;
} WaveField;

typedef struct {
    int nx, ny, nt;
    double dx, dt;
    int current_step;
} SimParams;

static SimParams g_params;

/* ==================== 必需函数 ==================== */
int _gettdsize_() { return sizeof(WaveField); }
int _getgdsize_() { return 64; }

/* ==================== 核心计算 ==================== */

static inline double get_u(double *u, int i, int j, int ny) {
    return u[i * ny + j];
}

static inline void set_u(double *u, int i, int j, int ny, double val) {
    u[i * ny + j] = val;
}

static void update_wavefield(WaveField *wf) {
    double courant2 = COURANT * COURANT;
    
    for (int i = wf->x_start; i < wf->x_end; i++) {
        for (int j = wf->y_start; j < wf->y_end; j++) {
            // 内部点计算（简化，不考虑边界）
            if (i <= 0 || i >= wf->nx-1 || j <= 0 || j >= wf->ny-1) continue;
            
            double u_ij = get_u(wf->u_curr, i, j, wf->ny);
            double u_im1 = get_u(wf->u_curr, i-1, j, wf->ny);
            double u_ip1 = get_u(wf->u_curr, i+1, j, wf->ny);
            double u_jm1 = get_u(wf->u_curr, i, j-1, wf->ny);
            double u_jp1 = get_u(wf->u_curr, i, j+1, wf->ny);
            
            double laplacian = u_ip1 + u_im1 + u_jp1 + u_jm1 - 4 * u_ij;
            double u_next = 2 * u_ij - get_u(wf->u_prev, i, j, wf->ny) + courant2 * laplacian;
            
            set_u(wf->u_next, i, j, wf->ny, u_next);
        }
    }
}

static void swap_buffers(WaveField *wf) {
    double *temp = wf->u_prev;
    wf->u_prev = wf->u_curr;
    wf->u_curr = wf->u_next;
    wf->u_next = temp;
}

static void init_wave(WaveField *wf, double cx, double cy) {
    for (int i = wf->x_start; i < wf->x_end; i++) {
        for (int j = wf->y_start; j < wf->y_end; j++) {
            double x = i * DX;
            double y = j * DX;
            double dist = sqrt((x-cx)*(x-cx) + (y-cy)*(y-cy));
            // 高斯波包
            double val = exp(-dist*dist / 0.002) * cos(10 * dist);
            set_u(wf->u_curr, i, j, wf->ny, val);
            set_u(wf->u_prev, i, j, wf->ny, val);
        }
    }
}

/* ==================== 输出功能 ==================== */

static void save_slice(const char *filename, WaveField *wf) {
    FILE *fp = fopen(filename, "w");
    if (!fp) return;
    
    // 输出中间切片的1D波形
    int mid_y = NY / 2;
    if (mid_y >= wf->y_start && mid_y < wf->y_end) {
        for (int i = 0; i < wf->nx; i++) {
            fprintf(fp, "%f %f\n", i * DX, get_u(wf->u_curr, i, mid_y, wf->ny));
        }
    }
    fclose(fp);
}

static void save_binary(const char *filename, WaveField *wf) {
    FILE *fp = fopen(filename, "wb");
    if (!fp) return;
    
    // 保存元数据
    int dims[2] = {wf->nx, wf->ny};
    fwrite(dims, sizeof(int), 2, fp);
    
    // 保存波场数据
    for (int i = 0; i < wf->nx; i++) {
        for (int j = 0; j < wf->ny; j++) {
            double val = get_u(wf->u_curr, i, j, wf->ny);
            fwrite(&val, sizeof(double), 1, fp);
        }
    }
    fclose(fp);
}

static void print_ascii_art(WaveField *wf) {
    // 简单的ASCII可视化
    int step = NY / 20;
    printf("\nWave Field (step %d):\n", g_params.current_step);
    printf("+");
    for (int i = 0; i < 40; i++) printf("-");
    printf("+\n");
    
    for (int j = 0; j < NY; j += step) {
        printf("|");
        for (int i = 0; i < NX; i += step) {
            double val = get_u(wf->u_curr, i, j, wf->ny);
            char c = ' ';
            if (val > 0.5) c = '#';
            else if (val > 0.2) c = '*';
            else if (val > 0.0) c = '.';
            else if (val > -0.2) c = ',';
            else if (val > -0.5) c = '+';
            else c = '-';
            printf("%c", c);
        }
        printf("|\n");
    }
    
    printf("+");
    for (int i = 0; i < 40; i++) printf("-");
    printf("+\n");
}

/* ==================== 线程函数 ==================== */

static void worker_thread() {
    int tid = gettid();
    WaveField *wf = (WaveField*)ti->td;
    
    for (int step = 0; step < g_params.nt; step++) {
        SWG;
        update_wavefield(wf);
        SSG;
    }
}

static void main_thread() {
    WaveField *wf = (WaveField*)ti->td;
    
    // 初始化波源
    init_wave(wf, 0.5, 0.5);
    
    printf("[Main] Starting wave propagation simulation...\n");
    
    for (int step = 0; step < g_params.nt; step++) {
        g_params.current_step = step;
        
        MSG;
        
        // 组主也参与计算
        update_wavefield(wf);
        
        MWG;
        swap_buffers(wf);
        
        // 输出
        if (step % SAVE_EVERY == 0) {
            char fname[256];
            snprintf(fname, sizeof(fname), "wave_step_%04d.txt", step);
            save_slice(fname, wf);
            
            if (mpi_id == 0 && step % (SAVE_EVERY * 2) == 0) {
                print_ascii_art(wf);
            }
        }
        
        MSS;
    }
    
    // 保存最终结果
    save_binary("wave_final.bin", wf);
    printf("[Main] Simulation complete. Results saved.\n");
}

static void wave_simulation() {
    int tid = gettid();
    
    if (tid == 0) {
        main_thread();
    } else {
        worker_thread();
    }
}

/* ==================== 初始化 ==================== */

static void alloc_wavefield(WaveField *wf, int nx, int ny) {
    wf->nx = nx;
    wf->ny = ny;
    wf->u_curr = (double*)calloc(nx * ny, sizeof(double));
    wf->u_prev = (double*)calloc(nx * ny, sizeof(double));
    wf->u_next = (double*)calloc(nx * ny, sizeof(double));
    wf->max_amplitude = 0.0;
}

static void setup_thread_data() {
    int nx_per_thread = NX / (N_THREADS + 1);
    
    for (int g = 0; g < md.ngrp; g++) {
        threadGroup *pg = md.grps[g];
        for (int t = 0; t <= N_THREADS; t++) {  // 包含组主
            THREADINFO *pti = &pg->threads[t];
            WaveField *wf = (WaveField*)pti->td;
            
            wf->x_start = t * nx_per_thread;
            wf->x_end = (t == N_THREADS) ? NX : (t + 1) * nx_per_thread;
            wf->y_start = 0;
            wf->y_end = NY;
            
            alloc_wavefield(wf, NX, NY);
        }
    }
}

/* ==================== 主函数 ==================== */

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);
    
    g_params.nx = NX;
    g_params.ny = NY;
    g_params.nt = NT;
    g_params.dx = DX;
    g_params.dt = DT;
    
    if (mpi_id == 0) {
        printf("2D Wave Propagation Visualization\n");
        printf("Grid: %d x %d, Steps: %d\n", NX, NY, NT);
    }
    
    int manage_id = -1;
    InitThreads(mpi_id, 16, 2, 16, N_THREADS + 1, N_GROUPS, 1, &manage_id);
    setup_thread_data();
    
    StartThreads(wave_simulation);
    wave_simulation();
    EndThreads();
    
    MPI_Finalize();
    return 0;
}
