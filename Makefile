# Makefile for mythread examples

CC = mpicc
# CFLAGS = -Wall -g -O2 -I.
CFLAGS = -g -O2 -I.
LDFLAGS = -lpthread -lm

# 目标文件
TARGETS = wave_propagation wave_propagation_ghost

# 源文件
MYTHREAD_SRC = mythread/mythread.c

.PHONY: all clean

all: $(TARGETS)
 
# 简化示例
example_simple: example_simple.c $(MYTHREAD_SRC)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# 完整示例
example_mpi_mythread: example_mpi_mythread.c $(MYTHREAD_SRC)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# 波传播模拟（原始版本）
wave_propagation: wave_propagation.c $(MYTHREAD_SRC)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# 波传播模拟（带幽灵层版本 - 推荐）
wave_propagation_ghost: wave_propagation_ghost.c $(MYTHREAD_SRC)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# 波传播可视化版本
wave_visual: wave_visual.c $(MYTHREAD_SRC)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# 清理
clean:
	rm -f $(TARGETS) *.o *.txt *.bin *.png *.gif

# 运行示例
run_simple: example_simple
	mpirun -np 2 ./example_simple

run_full: example_mpi_mythread
	mpirun -np 4 ./example_mpi_mythread

run_wave: wave_propagation
	mpirun -np 4 ./wave_propagation

# 运行带幽灵层版本（推荐）
run_ghost: wave_propagation_ghost
	mpirun -np 4 ./wave_propagation_ghost

run_visual: wave_visual
	mpirun -np 1 ./wave_visual

# 可视化
visualize:
	python3 visualize_wave.py --all

# 调试模式编译
debug: CFLAGS += -DDEBUG -g
debug: all
