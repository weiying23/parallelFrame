# Makefile for mythread examples

CC = mpicc
CFLAGS = -Wall -g -O2 -I.. -I.
LDFLAGS = -lpthread -lm

# 目标文件
TARGETS = example_simple example_mpi_mythread

# 源文件
MYTHREAD_SRC = ../mythread/mythread.c

.PHONY: all clean

all: $(TARGETS)

# 简化示例
example_simple: example_simple.c $(MYTHREAD_SRC)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# 完整示例
example_mpi_mythread: example_mpi_mythread.c $(MYTHREAD_SRC)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# 清理
clean:
	rm -f $(TARGETS) *.o

# 运行示例
run_simple: example_simple
	mpirun -np 2 ./example_simple

run_full: example_mpi_mythread
	mpirun -np 4 ./example_mpi_mythread

# 调试模式编译
debug: CFLAGS += -DDEBUG -g
debug: all
