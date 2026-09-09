#!/bin/bash
# mythread 三模式多网格性能基准测试
# 用法: ./scripts/benchmark_full.sh [quick|full]

set -e
cd "$(dirname "$0")/.."
MODE="${1:-full}"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOGDIR="logs/bench_${TIMESTAMP}"
mkdir -p "$LOGDIR"
LOGFILE="${LOGDIR}/benchmark.log"
CSVFILE="${LOGDIR}/results.csv"
N_RUNS=5

# 颜色
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'
BLUE='\033[0;34m'; CYAN='\033[0;36m'; NC='\033[0m'; BOLD='\033[1m'

log()  { echo -e "$@" | tee -a "$LOGFILE"; }
log_n(){ echo -n "$@" | tee -a "$LOGFILE"; }

# ── 网格配置 ──
if [ "$MODE" = "quick" ]; then
    GRIDS=(
        "800x800x100"
        "2000x2000x100"
        "4000x4000x100"
    )
else
    GRIDS=(
        "1000x1000x100"
        "2000x2000x200"
        "4000x4000x200"
        "6000x6000x200"
        "8000x8000x200"
    )
fi

# ── 编译 ──
log "${BOLD}${CYAN}══════════════════════════════════════════════${NC}"
log "${BOLD}${CYAN}  mythread 多模式性能基准测试${NC}"
log "${BOLD}${CYAN}  时间: $(date)  模式: ${MODE}  重复: ${N_RUNS}次${NC}"
log "${BOLD}${CYAN}══════════════════════════════════════════════${NC}"
log ""

log "${YELLOW}>>> 编译所有版本...${NC}"
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release >/dev/null 2>&1
cmake --build build -j >/dev/null 2>&1
log "    fix / ghost 完成"

OMP_OK=0
if [ -x build/wave_propagation_ghost_omp ]; then
    OMP_OK=1; log "    omp  完成"
else
    log "    omp  跳过 (未找到 OpenMP)"
fi
log ""

# ── CSV 头 ──
echo "grid,binary,np,run,time_s,throughput_mpups,energy,l2_init" > "$CSVFILE"

# ── 辅助函数 ──
run_one() {
    local grid="$1" binary="$2" np="$3" run_id="$4" extra="$5"
    local nx=$(echo "$grid" | cut -d'x' -f1)
    local ny=$(echo "$grid" | cut -d'x' -f2)
    local nt=$(echo "$grid" | cut -d'x' -f3)

    # 生成临时配置文件
    cat > /tmp/bench_case.cfg << EOF
NX = $nx
NY = $ny
NT = $nt
A = 10.0
DT = 0.05
C0 = 0.1
USE_FIXED_DOMAIN = 0
HALO = 1
ENERGY_REPORT_INTERVAL = $((nt/2))
EOF

    local t e l2
    local output
    output=$(eval "$extra WAVE_CASE_CFG=/tmp/bench_case.cfg WAVE_HARDWARE_CFG=config/hardware.cfg mpirun -np $np --allow-run-as-root ./build/$binary" 2>&1)
    t=$(echo "$output" | grep 'Simulation completed\|Total wall' | grep -o '[0-9]*\.[0-9]*' | tail -1)
    e=$(echo "$output" | grep 'Step.*'"$nt" | grep -o 'E=[0-9.]*' | head -1 | cut -d'=' -f2)
    l2=$(echo "$output" | grep 'Initial:' | grep -o 'L2=[0-9.]*' | head -1 | cut -d'=' -f2)
    [ -z "$t" ] && t="0.000"
    [ -z "$e" ] && e="0"
    [ -z "$l2" ] && l2="0"

    local mps="0"
    if [ "$t" != "0.000" ]; then
        mps=$(echo "scale=1; $nx*$ny*$nt/$t/1000000" | bc 2>/dev/null || echo "0")
    fi

    echo "$grid,$binary,np=$np,$run_id,$t,$mps,$e,$l2" >> "$CSVFILE"
    printf "%s" "$t"
}

# ── 运行所有测试 ──
TMPFILE=$(mktemp /tmp/bench_results.XXXXXX)

for grid in "${GRIDS[@]}"; do
    nx=$(echo "$grid" | cut -d'x' -f1)
    ny=$(echo "$grid" | cut -d'x' -f2)
    nt=$(echo "$grid" | cut -d'x' -f3)
    gname="${nx}×${ny}×${nt}"
    points=$((nx * ny * nt))

    log ""
    log "${BOLD}${BLUE}┌──────────────────────────────────────────────┐${NC}"
    log "${BOLD}${BLUE}│  网格: ${gname}  (${points} 点·步)${NC}"
    log "${BOLD}${BLUE}└──────────────────────────────────────────────┘${NC}"

    # 存储每次运行的时间
    declare -a all_runs

    for np in 1; do
        for binary in wave_propagation_ghost_fix wave_propagation_ghost; do
            name="${binary##*_}"
            best=999999 total=0 worst=0
            unset all_runs

            log_n "  ${name}-8s np=${np}: "
            for run in $(seq 1 $N_RUNS); do
                t=$(run_one "$grid" "$binary" "$np" "$run" "")
                all_runs[$run]=$t
                log_n "${t}s "
                total=$(echo "$total + $t" | bc)
                [ "$(echo "$t < $best" | bc -l 2>/dev/null)" = "1" ] && best=$t
                [ "$(echo "$t > $worst" | bc -l 2>/dev/null)" = "1" ] && worst=$t
            done

            avg=$(echo "scale=3; $total / $N_RUNS" | bc)
            mps=$(echo "scale=1; $nx*$ny*$nt/$best/1000000" | bc)
            # 构建各次运行字符串
            runs_str=$(printf "%s " "${all_runs[@]}")
            log " →  min=${best}s  avg=${avg}s  max=${worst}s  [${runs_str}] ${mps} Mpts/s"

            echo "${grid}|${name}|np${np}|min=${best}|avg=${avg}|max=${worst}|runs=${runs_str}" >> "$TMPFILE"
        done

        # OMP 版本
        if [ "$OMP_OK" = "1" ]; then
            name="omp"; best=999999 total=0 worst=0
            unset all_runs

            ng=$(grep N_GROUPS config/hardware.cfg 2>/dev/null | awk '{print $3}' || echo 2)
            nw=$(grep N_WORKERS config/hardware.cfg 2>/dev/null | awk '{print $3}' || echo 3)
            nth=$((ng * nw))

            log_n "  ${name}-8s np=${np}: "
            for run in $(seq 1 $N_RUNS); do
                t=$(run_one "$grid" "wave_propagation_ghost_omp" "$np" "$run" "OMP_NUM_THREADS=$nth")
                all_runs[$run]=$t
                log_n "${t}s "
                total=$(echo "$total + $t" | bc)
                [ "$(echo "$t < $best" | bc -l 2>/dev/null)" = "1" ] && best=$t
                [ "$(echo "$t > $worst" | bc -l 2>/dev/null)" = "1" ] && worst=$t
            done

            avg=$(echo "scale=3; $total / $N_RUNS" | bc)
            mps=$(echo "scale=1; $nx*$ny*$nt/$best/1000000" | bc)
            runs_str=$(printf "%s " "${all_runs[@]}")
            log " →  min=${best}s  avg=${avg}s  max=${worst}s  [${runs_str}] ${mps} Mpts/s"

            echo "${grid}|omp|np${np}|min=${best}|avg=${avg}|max=${worst}|runs=${runs_str}" >> "$TMPFILE"
        fi
    done
done

# ═══════════════════════════════════════════════════════════════
#  汇总表格
# ═══════════════════════════════════════════════════════════════
log ""
log "${BOLD}${CYAN}╔══════════════════════════════════════════════════════════════╗${NC}"
log "${BOLD}${CYAN}║                    汇总结果表格                              ║${NC}"
log "${BOLD}${CYAN}╠══════════════════════════════════════════════════════════════╣${NC}"

# 辅助：从临时文件中读取某网格某版本的最佳时间
# 读取某网格某版本的 min/avg/max
get_field() {
    _g="$1"; _b="$2"; _f="$3"
    grep "^${_g}|${_b}|" "$TMPFILE" 2>/dev/null | head -1 | sed "s/.*${_f}=\([0-9.]*\).*/\1/"
}

# ── 表 1: 汇总 (avg ± spread) ──
printf "${BOLD}%-20s %18s %18s %18s${NC}\n" "网格" "fix (min/avg/max)" "ghost (min/avg/max)" "omp (min/avg/max)"
log "${BOLD}──────────────────────────────────────────────────────────────────────${NC}"

for grid in "${GRIDS[@]}"; do
    f_min=$(get_field "$grid" "fix" "min"); f_avg=$(get_field "$grid" "fix" "avg"); f_max=$(get_field "$grid" "fix" "max")
    g_min=$(get_field "$grid" "ghost" "min"); g_avg=$(get_field "$grid" "ghost" "avg"); g_max=$(get_field "$grid" "ghost" "max")
    o_min=$(get_field "$grid" "omp" "min"); o_avg=$(get_field "$grid" "omp" "avg"); o_max=$(get_field "$grid" "omp" "max")

    f_str="${f_min:-?}/${f_avg:-?}/${f_max:-?}"
    g_str="${g_min:-?}/${g_avg:-?}/${g_max:-?}"
    o_str="${o_min:-?}/${o_avg:-?}/${o_max:-?}"

    printf "%-20s ${GREEN}%18s${NC} ${YELLOW}%18s${NC} ${CYAN}%18s${NC}\n" "$grid" "$f_str" "$g_str" "$o_str"
done

log ""
log "${BOLD}  格式: min/avg/max (秒)${NC}"
log ""

# ── 表 2: 相对性能 ──
printf "${BOLD}%-20s %10s %10s %10s %8s %8s${NC}\n" \
    "网格" "fix(avg)" "ghost(avg)" "omp(avg)" "fix/omp" "ghost/fix"
log "${BOLD}──────────────────────────────────────────────────────────────────${NC}"

for grid in "${GRIDS[@]}"; do
    f_avg=$(get_field "$grid" "fix" "avg")
    g_avg=$(get_field "$grid" "ghost" "avg")
    o_avg=$(get_field "$grid" "omp" "avg")
    [ -z "$f_avg" ] && f_avg="N/A"
    [ -z "$g_avg" ] && g_avg="N/A"
    [ -z "$o_avg" ] && o_avg="N/A"

    r1="-" r2="-"
    if [ "$f_avg" != "N/A" ] && [ "$o_avg" != "N/A" ] && [ "$o_avg" != "0.000" ]; then
        r1=$(echo "scale=2; $f_avg / $o_avg" | bc 2>/dev/null || echo "-")
    fi
    if [ "$g_avg" != "N/A" ] && [ "$f_avg" != "N/A" ] && [ "$f_avg" != "0.000" ]; then
        r2=$(echo "scale=2; $g_avg / $f_avg" | bc 2>/dev/null || echo "-")
    fi

    printf "%-20s ${GREEN}%10s${NC} ${YELLOW}%10s${NC} ${CYAN}%10s${NC} %8s %8s\n" \
        "$grid" "$f_avg" "$g_avg" "$o_avg" "${r1}x" "${r2}x"
done

log "${BOLD}──────────────────────────────────────────────────────────────────${NC}"
log "  fix/omp: fix 相对 omp 的倍数 (>1 表示 omp 更快)"
log "  ghost/fix: ghost 相对 fix 的倍数 (>1 表示 fix 更快)"
log ""

# ═══════════════════════════════════════════════════════════════
log "${BOLD}${CYAN}╚══════════════════════════════════════════════════════════════╝${NC}"
log ""
log "${GREEN}日志目录: ${LOGDIR}${NC}"
log "${GREEN}详细 CSV: ${CSVFILE}${NC}"
log ""

# 输出 CSV 内容摘要
log "${CYAN}── CSV 数据预览 (前 10 行) ──${NC}"
head -10 "$CSVFILE" | column -t -s',' | tee -a "$LOGFILE"
log ""
log "完成。完整结果见 ${LOGDIR}/"
rm -f "$TMPFILE"
