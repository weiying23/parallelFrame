#!/usr/bin/env python3
"""
波场数据可视化脚本

用法：
    python3 visualize_wave.py
    python3 visualize_wave.py --animate
    python3 visualize_wave.py --3d

功能：
    - 绘制波场快照
    - 创建动画
    - 3D 表面图
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import glob
import argparse

def load_binary_data(filename):
    """加载二进制波场数据"""
    with open(filename, 'rb') as f:
        dims = np.fromfile(f, dtype=np.int32, count=2)
        nx, ny = dims[0], dims[1]
        data = np.fromfile(f, dtype=np.float64, count=nx*ny)
        return data.reshape((nx, ny))

def load_slice_data(filename):
    """加载切片文本数据"""
    return np.loadtxt(filename)

def plot_wave_snapshot(filename, ax=None, title=None):
    """绘制单个波场快照"""
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 8))
    
    data = load_binary_data(filename)
    
    im = ax.imshow(data.T, origin='lower', cmap='RdBu_r', 
                   vmin=-1, vmax=1, interpolation='bilinear')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    if title:
        ax.set_title(title)
    plt.colorbar(im, ax=ax, label='Amplitude')
    
    return ax

def plot_all_snapshots():
    """绘制所有保存的切片"""
    files = sorted(glob.glob('wave_step_*.txt'))
    
    if not files:
        print("No wave data files found!")
        return
    
    n_files = len(files)
    n_cols = 3
    n_rows = (n_files + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 5*n_rows))
    axes = axes.flatten() if n_files > 1 else [axes]
    
    for i, fname in enumerate(files):
        data = load_slice_data(fname)
        step = int(fname.split('_')[2].split('.')[0])
        
        axes[i].plot(data[:, 0], data[:, 1])
        axes[i].set_title(f'Step {step}')
        axes[i].set_xlabel('X')
        axes[i].set_ylabel('Amplitude')
        axes[i].set_ylim(-1.5, 1.5)
        axes[i].grid(True)
    
    # 隐藏多余的子图
    for i in range(n_files, len(axes)):
        axes[i].axis('off')
    
    plt.tight_layout()
    plt.savefig('wave_slices.png', dpi=150)
    print("Saved: wave_slices.png")
    plt.show()

def animate_wave():
    """创建波传播动画"""
    files = sorted(glob.glob('wave_step_*.txt'))
    
    if not files:
        print("No wave data files found!")
        return
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # 初始化
    data = load_slice_data(files[0])
    line, = ax.plot(data[:, 0], data[:, 1], 'b-', linewidth=2)
    ax.set_xlim(0, data[-1, 0])
    ax.set_ylim(-1.5, 1.5)
    ax.set_xlabel('Position (m)')
    ax.set_ylabel('Wave Amplitude')
    ax.set_title('Wave Propagation')
    ax.grid(True)
    
    step_text = ax.text(0.02, 0.95, '', transform=ax.transAxes,
                       fontsize=12, verticalalignment='top',
                       bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    def init():
        line.set_ydata(np.zeros_like(data[:, 0]))
        step_text.set_text('')
        return line, step_text
    
    def update(frame):
        data = load_slice_data(files[frame])
        line.set_ydata(data[:, 1])
        step = int(files[frame].split('_')[2].split('.')[0])
        step_text.set_text(f'Step: {step}')
        return line, step_text
    
    anim = FuncAnimation(fig, update, frames=len(files),
                        init_func=init, blit=True, interval=200)
    
    anim.save('wave_animation.gif', writer='pillow', fps=5)
    print("Saved: wave_animation.gif")
    plt.show()

def plot_3d_surface(filename='wave_final.bin'):
    """绘制3D表面图"""
    try:
        from mpl_toolkits.mplot3d import Axes3D
    except:
        print("3D plotting not available")
        return
    
    data = load_binary_data(filename)
    nx, ny = data.shape
    
    x = np.linspace(0, 1, nx)
    y = np.linspace(0, 1, ny)
    X, Y = np.meshgrid(x, y)
    
    fig = plt.figure(figsize=(12, 5))
    
    # 3D 表面
    ax1 = fig.add_subplot(121, projection='3d')
    surf = ax1.plot_surface(X, Y, data.T, cmap='RdBu_r', 
                           vmin=-1, vmax=1, antialiased=True)
    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    ax1.set_zlabel('Amplitude')
    ax1.set_title('3D Wave Surface')
    fig.colorbar(surf, ax=ax1, shrink=0.5, aspect=5)
    
    # 2D 热图
    ax2 = fig.add_subplot(122)
    im = ax2.imshow(data.T, origin='lower', cmap='RdBu_r',
                   vmin=-1, vmax=1, extent=[0, 1, 0, 1])
    ax2.set_xlabel('X')
    ax2.set_ylabel('Y')
    ax2.set_title('Wave Amplitude Map')
    fig.colorbar(im, ax=ax2)
    
    plt.tight_layout()
    plt.savefig('wave_3d.png', dpi=150)
    print("Saved: wave_3d.png")
    plt.show()

def analyze_energy():
    """分析波场能量分布"""
    data = load_binary_data('wave_final.bin')
    
    # 计算能量密度
    energy = np.sum(data**2)
    max_amp = np.max(np.abs(data))
    
    print(f"\nWave Field Statistics:")
    print(f"  Total Energy: {energy:.6f}")
    print(f"  Max Amplitude: {max_amp:.6f}")
    print(f"  Mean Amplitude: {np.mean(np.abs(data)):.6f}")
    print(f"  RMS Amplitude: {np.sqrt(np.mean(data**2)):.6f}")
    
    # 绘制能量分布直方图
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    
    axes[0].hist(data.flatten(), bins=100, color='blue', alpha=0.7)
    axes[0].set_xlabel('Amplitude')
    axes[0].set_ylabel('Count')
    axes[0].set_title('Amplitude Distribution')
    axes[0].set_yscale('log')
    
    # 径向能量分布
    nx, ny = data.shape
    cx, cy = nx // 2, ny // 2
    max_r = min(cx, cy)
    radial_energy = []
    radii = range(0, max_r, max(1, max_r // 50))
    
    for r in radii:
        mask = np.zeros_like(data, dtype=bool)
        for i in range(nx):
            for j in range(ny):
                if abs(np.sqrt((i-cx)**2 + (j-cy)**2) - r) < 1:
                    mask[i, j] = True
        if np.any(mask):
            radial_energy.append(np.mean(data[mask]**2))
        else:
            radial_energy.append(0)
    
    axes[1].plot(radii, radial_energy, 'r-', linewidth=2)
    axes[1].set_xlabel('Radius (pixels)')
    axes[1].set_ylabel('Mean Energy Density')
    axes[1].set_title('Radial Energy Distribution')
    axes[1].grid(True)
    
    plt.tight_layout()
    plt.savefig('wave_energy_analysis.png', dpi=150)
    print("Saved: wave_energy_analysis.png")
    plt.show()

def main():
    parser = argparse.ArgumentParser(description='Wave Field Visualization')
    parser.add_argument('--animate', action='store_true', 
                       help='Create animation')
    parser.add_argument('--3d', action='store_true',
                       help='Plot 3D surface')
    parser.add_argument('--energy', action='store_true',
                       help='Analyze energy distribution')
    parser.add_argument('--all', action='store_true',
                       help='Generate all plots')
    
    args = parser.parse_args()
    
    if args.all:
        args.animate = args.energy = True
        plot_3d_surface()
    
    if args.animate:
        animate_wave()
    elif args.energy:
        analyze_energy()
    elif args._3d:
        plot_3d_surface()
    else:
        plot_all_snapshots()

if __name__ == '__main__':
    main()
