CC = mpicc
CFLAGS = -g -O2 -I.
LDFLAGS = -lpthread -lm

MYTHREAD_SRC = mythread/mythread_util.c \
               mythread/mythread_locv.c \
               mythread/mythread_timer.c \
               mythread/mythread_sync.c \
               mythread/mythread_pool.c \
               mythread/mythread_thread.c \
               mythread/mythread_fortran.c \
               mythread/mythread_decomp.c \
               mythread/mythread_field.c \
               mythread/mythread_halo.c \
               mythread/mythread_config.c
TARGETS = wave_propagation wave_propagation_ghost tests/taskpool_tests tests/test_bindcpu tests/test_field tests/test_config tests/test_decomp tests/test_locv tests/test_sync tests/test_timer tests/test_util

.PHONY: all clean run_simple run_full run_wave run_ghost run_visual run_taskpool_tests run_bindcpu_tests debug

all: $(TARGETS)

example_simple: example_simple.c $(MYTHREAD_SRC)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

example_mpi_mythread: example_mpi_mythread.c $(MYTHREAD_SRC)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

example_taskpool: example_taskpool.c $(MYTHREAD_SRC)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

tests/taskpool_tests: tests/taskpool_tests.c $(MYTHREAD_SRC)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

tests/test_bindcpu: tests/test_bindcpu.c $(MYTHREAD_SRC)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

tests/test_field: tests/test_field.c $(MYTHREAD_SRC)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

tests/test_config: tests/test_config.c $(MYTHREAD_SRC)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

tests/test_decomp: tests/test_decomp.c $(MYTHREAD_SRC)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

tests/test_locv: tests/test_locv.c $(MYTHREAD_SRC)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

tests/test_sync: tests/test_sync.c $(MYTHREAD_SRC)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

tests/test_timer: tests/test_timer.c $(MYTHREAD_SRC)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

tests/test_util: tests/test_util.c $(MYTHREAD_SRC)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

wave_propagation: wave_propagation.c $(MYTHREAD_SRC)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

wave_propagation_ghost: wave_propagation_ghost.c $(MYTHREAD_SRC)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

wave_visual: wave_visual.c $(MYTHREAD_SRC)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm -f $(TARGETS) *.o *.txt *.bin *.png *.gif mythread/*.o

run_simple: example_simple
	mpirun -np 2 ./example_simple

run_full: example_mpi_mythread
	mpirun -np 4 ./example_mpi_mythread

run_wave: wave_propagation
	mpirun -np 4 ./wave_propagation

run_ghost: wave_propagation_ghost
	mpirun -np 4 ./wave_propagation_ghost

run_visual: wave_visual
	mpirun -np 1 ./wave_visual

run_taskpool_tests: tests/taskpool_tests
	mpirun -np 1 ./tests/taskpool_tests grouped-basic
	mpirun -np 1 ./tests/taskpool_tests grouped-block-full
	mpirun -np 1 ./tests/taskpool_tests grouped-zero-task
	mpirun -np 1 ./tests/taskpool_tests grouped-isolation
	mpirun -np 1 ./tests/taskpool_tests grouped-invalid-order
	mpirun -np 1 ./tests/taskpool_tests grouped-try-full
	mpirun -np 1 ./tests/taskpool_tests ungrouped-basic
	mpirun -np 1 ./tests/taskpool_tests boundary

run_bindcpu_tests: tests/test_bindcpu
	./tests/test_bindcpu

run_field_tests: tests/test_field
	./tests/test_field

run_config_tests: tests/test_config
	./tests/test_config

run_decomp_tests: tests/test_decomp
	./tests/test_decomp

run_locv_tests: tests/test_locv
	./tests/test_locv

run_sync_tests: tests/test_sync
	./tests/test_sync

run_timer_tests: tests/test_timer
	./tests/test_timer

run_util_tests: tests/test_util
	./tests/test_util

visualize:
	python3 visualize_wave.py --all

debug: CFLAGS += -DDEBUG -g
debug: all
