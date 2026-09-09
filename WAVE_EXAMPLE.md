# 波能量传播并行计算示例

这个示例展示了如何使用 mythread + MPI 实现 2D 波动方程的并行求解。

## 物理模型

### 波动方程

```
∂²u/∂t² = c²(∂²u/∂x² + ∂²u/∂y²)
```

其中：
- `u(x,y,t)`：波场振幅
- `c`：波速
- `t`：时间
- `x, y`：空间坐标

### 有限差分离散

使用中心差分格式：

```
u(t+1,i,j) = 2*u(t,i,j) - u(t-1,i,j) + 
             (c*dt/dx)² * [u(t,i+1,j) + u(t,i-1,j) + 
                          u(t,i,j+1) + u(t,i,j-1) - 4*u(t,i,j)]
```

## 文件说明

| 文件 | 说明 |
|------|------|
| `wave_propagation.c` | 完整的波传播模拟（高性能版本） |
| `wave_visual.c` | 带可视化输出的简化版本 |
| `visualize_wave.py` | Python 可视化脚本 |
| `CMakeLists.txt` | CMake 构建脚本 |

## 编译运行

### 编译

使用 CMake 构建（默认启用 MPI）：

```bash
cmake -S . -B build
cmake --build build -j

# 需要指定 MPI 编译器时
cmake -S . -B build -DCMAKE_C_COMPILER=mpicc
cmake --build build -j
```

若不需要 MPI，可关闭 MPI 只构建非 MPI 目标：

```bash
cmake -S . -B build -DUSE_MPI=OFF
cmake --build build -j
```

### 运行

```bash
# 运行高性能版本（4个MPI进程）
mpirun -np 4 ./build/wave_propagation

# 运行可视化版本（1个进程，输出文件）
mpirun -np 1 ./build/wave_visual

# 波动方程 ghost 版本
mpirun -np 4 ./build/wave_propagation_ghost
mpirun -np 4 ./build/wave_propagation_ghost_fix
```

### 可视化

```bash
# 需要安装依赖
pip install numpy matplotlib

# 生成所有可视化
python3 visualize_wave.py --all

# 只生成动画
python3 visualize_wave.py --animate

# 3D表面图
python3 visualize_wave.py --3d

# 能量分析
python3 visualize_wave.py --energy
```

## 并行架构

```
┌─────────────────────────────────────────────┐
│           MPI 进程 0                          │
│  ┌─────────────────────────────────────┐   │
│  │         线程组 0 (Group 0)           │   │
│  │  ┌─────┐ ┌─────┐ ┌─────┐ ┌─────┐  │   │
│  │  │Main │ │ W1  │ │ W2  │ │ W3  │  │   │
│  │  │(GWM)│ │(SWG)│ │(SWG)│ │(SWG)│  │   │
│  │  │(GSS)│ │(SSG)│ │(SSG)│ │(SSG)│  │   │
│  │  │(GSM)│ │     │ │     │ │     │  │   │
│  │  └─────┘ └─────┘ └─────┘ └─────┘  │   │
│  │      X方向划分：每个线程负责一段     │   │
│  └─────────────────────────────────────┘   │
│           ↑ MSG/MWG ↓                      │
│  ┌─────────────────────────────────────┐   │
│  │         线程组 1 (Group 1)           │   │
│  │  ┌─────┐ ┌─────┐ ┌─────┐ ┌─────┐  │   │
│  │  │Main │ │ W1  │ │ W2  │ │ W3  │  │   │
│  │  │(GWM)│ │(SWG)│ │(SWG)│ │(SWG)│  │   │
│  │  └─────┘ └─────┘ └─────┘ └─────┘  │   │
│  └─────────────────────────────────────┘   │
└─────────────────────────────────────────────┘
                      │ MPI通信
                      ↓
┌─────────────────────────────────────────────┐
│           MPI 进程 1                          │
│              [相同结构]                       │
└─────────────────────────────────────────────┘
```

### 数据划分

- **Y方向**：按线程组划分（不同组处理不同的Y范围）
- **X方向**：组内按线程划分（组内线程分担X方向计算）

```
┌─────────────────────────────────┐
│  Group 0: Y=[0, NY/2)          │
│  ┌─────────┬─────────┬─────────┐│
│  │Thread 0 │Thread 1 │Thread 2 ││
│  │X:[0,n)  │X:[n,2n) │X:[2n,NX)││
│  └─────────┴─────────┴─────────┘│
├─────────────────────────────────┤
│  Group 1: Y=[NY/2, NY)         │
│  ┌─────────┬─────────┬─────────┐│
│  │Thread 0 │Thread 1 │Thread 2 ││
│  │X:[0,n)  │X:[n,2n) │X:[2n,NX)││
│  └─────────┴─────────┴─────────┘│
└─────────────────────────────────┘
```

## 同步流程

每个时间步的同步流程：

```
时间步 n:
    
主线程                    组主线程                  工作线程
   │                          │                        │
   │──── MSG ────────────────►│                        │
   │                          │──── GSS ──────────────►│
   │                          │                        │──┐
   │                          │                        │  │ 计算
   │                          │◄─── SSG ──────────────│──┘
   │                          │                        │
   │                          │──┐                     │
   │                          │  │ 计算                │
   │                          │◄─┘                     │
   │                          │                        │
   │◄─── GSM ─────────────────│                        │
   │                          │                        │
◄──┴──► MPI_Allreduce         │                        │
   │                          │                        │
   │──── MSS ────────────────►│                        │
   │                          │                        │
   ▼                          ▼                        ▼
时间步 n+1
```

## 关键代码解释

### 1. 线程数据设置

```c
static void setup_thread_data() {
    for (int g = 0; g < md.ngrp; g++) {
        for (int t = 0; t <= N_THREADS; t++) {
            THREADINFO *pti = &pg->threads[t];
            WaveField *wf = (WaveField*)pti->td;
            
            // 设置计算区域
            wf->x_start = t * nx_per_thread;
            wf->x_end = (t + 1) * nx_per_thread;
            wf->y_start = g * ny_per_group;
            wf->y_end = (g + 1) * ny_per_group;
            
            alloc_wavefield(wf, NX, NY);
        }
    }
}
```

### 2. 工作线程函数

```c
static void worker_thread() {
    WaveField *wf = (WaveField*)ti->td;
    
    for (int step = 0; step < g_params.nt; step++) {
        SWG;                    // 等待组主信号
        update_wavefield(wf);   // 更新波场
        SSG;                    // 通知完成
    }
}
```

### 3. 组主线程函数

```c
static void group_main_thread() {
    for (int step = 0; step < g_params.nt; step++) {
        GWM;                    // 等待主线程
        GSS;                    // 通知工作线程
        update_wavefield(wf);   // 组主也参与计算
        GWS;                    // 等待工作线程
        swap_buffers(wf);       // 交换时间层
        GSM;                    // 通知主线程
    }
}
```

### 4. 主线程函数

```c
static void main_thread() {
    for (int step = 0; step < g_params.nt; step++) {
        MSG;                    // 通知所有组开始
        MWG;                    // 等待所有组完成
        
        // MPI 进程间通信
        MPI_Allreduce(&local_energy, &global_energy, ...);
        
        MSS;                    // 准备下一时间步
    }
}
```

## 性能优化建议

1. **负载均衡**：确保每个线程的计算量相近
2. **减少同步**：合并相邻的时间步计算（如果稳定性允许）
3. **向量化**：使用 SIMD 指令加速有限差分计算
4. **内存布局**：考虑使用 Structure of Arrays (SoA) 优化缓存

## 扩展功能

### 1. 吸收边界条件

在边界处添加吸收层，防止反射：

```c
// PML (Perfectly Matched Layer)
if (i < pml_width || i > nx - pml_width) {
    damping = exp(-sigma * depth);
    u_next *= damping;
}
```

### 2. 多种波源

```c
// 点源
void init_point_source(WaveField *wf, int cx, int cy) {
    wf->u_curr[cx * ny + cy] = 1.0;
}

// 线源
void init_line_source(WaveField *wf, int x) {
    for (int j = 0; j < ny; j++) {
        wf->u_curr[x * ny + j] = sin(j * 2 * M_PI / ny);
    }
}
```

### 3. 异质介质

```c
// 不同位置的波速不同
c[i][j] = c0 * (1 + 0.1 * sin(i * 0.1));
```

## 输出文件说明

| 文件 | 格式 | 内容 |
|------|------|------|
| `wave_step_*.txt` | 文本 | 中间切片的1D波形 |
| `wave_final.bin` | 二进制 | 最终2D波场 |
| `wave_slices.png` | 图像 | 所有切片对比 |
| `wave_animation.gif` | 动画 | 波传播过程 |
| `wave_3d.png` | 图像 | 3D表面图 |
