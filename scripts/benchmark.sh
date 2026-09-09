#!/bin/bash
set -e
N_RUNS=3
export WAVE_CASE_CFG="config/case_bench.cfg"
export WAVE_HARDWARE_CFG="config/hardware_bench.cfg"
cd "$(dirname "$0")/.."

echo "=== mythread 三版本性能对比 ==="
echo "网格: 1000x1000, NT=200, 每项跑${N_RUNS}次取中位数"
echo ""

echo ">>> 编译..."
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release >/dev/null 2>&1
cmake --build build -j >/dev/null 2>&1
echo "fix / ghost: OK"

OMP_OK=0
if [ -x build/wave_propagation_ghost_omp ]; then
    OMP_OK=1; echo "omp: OK"
else
    echo "omp: SKIP (未找到 OpenMP，跳过)"
fi
echo ""

run_bench() {
    local name="$1" binary="$2" np="$3" extra="$4"
    local best=999999
    for i in $(seq 1 $N_RUNS); do
        t=$(eval "$extra mpirun -np $np --allow-run-as-root ./build/$binary" 2>&1 | \
            grep 'Simulation completed\|Total wall' | grep -o '[0-9.]*' | tail -1)
        [ -n "$t" ] && [ "$(echo "$t < $best" | bc -l 2>/dev/null)" = "1" ] && best=$t
    done
    if [ "$best" != "999999" ]; then
        echo "$best $name np=$np" >> /tmp/bench_r.txt
        printf "  %-10s np=%-2d: %.3f s\n" "$name" "$np" "$best"
    else
        printf "  %-10s np=%-2d: FAILED\n" "$name" "$np"
    fi
}

rm -f /tmp/bench_r.txt

for np in 1 2; do
    echo ">>> np=$np"
    run_bench "fix" "wave_propagation_ghost_fix" "$np"
    run_bench "ghost" "wave_propagation_ghost" "$np"
    [ $OMP_OK -eq 1 ] && run_bench "omp" "wave_propagation_ghost_omp" "$np" "OMP_NUM_THREADS=4"
    echo ""
done

echo "=== 结果 ==="
if [ -f /tmp/bench_r.txt ] && [ -s /tmp/bench_r.txt ]; then
    BASELINE=$(sort -n /tmp/bench_r.txt | head -1 | awk '{print $1}')
    sort -n /tmp/bench_r.txt | while read t name rest; do
        ratio=$(echo "scale=2; $t / $BASELINE" | bc)
        printf "  %5sx  %.3f s  %s %s\n" "$ratio" "$t" "$name" "$rest"
    done
fi
rm -f /tmp/bench_r.txt
