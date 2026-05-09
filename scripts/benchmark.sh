#!/bin/bash
# 性能测试脚本 — 遍历参数组合，重复测试，汇总输出
# 用法: ./scripts/benchmark.sh [quick|full]
#   quick: 小网格快速对比（默认）
#   full:  完整参数矩阵

set -e
cd "$(dirname "$0")/.."

MODE="${1:-quick}"
REPEAT=3
LOG_DIR="logs/$(date +%Y%m%d_%H%M%S)"
SUMMARY="$LOG_DIR/summary.csv"
CASE_CFG="config/case.cfg"
HW_CFG="config/hardware.cfg"

mkdir -p "$LOG_DIR"

# ── 参数矩阵 ──
if [ "$MODE" = "quick" ]; then
  GRID_SIZES=("1000 1000 50")                    # NX NY NT
  GROUP_COUNTS=(1 2)                              # N_GROUPS
  WORKER_COUNTS=(1 2)                             # N_WORKERS
  DECOMPS=(0)                                     # 0=Y_ONLY, 1=XY_2D
  IMBALANCE=(0 50)                                # %
  PROGRAMS=("wave_propagation_ghost" "wave_propagation_ghost_fix")
elif [ "$MODE" = "full" ]; then
  GRID_SIZES=("500 500 100" "2000 2000 50")
  GROUP_COUNTS=(1 2 4)
  WORKER_COUNTS=(1 2 3)
  DECOMPS=(0 1)
  IMBALANCE=(0 30)
  PROGRAMS=("wave_propagation_ghost" "wave_propagation_ghost_fix")
else
  echo "Usage: $0 [quick|full]"
  exit 1
fi

# ── 备份原始配置 ──
cp "$CASE_CFG" "$LOG_DIR/case.cfg.bak"
cp "$HW_CFG"  "$LOG_DIR/hardware.cfg.bak"

# ── 写配置辅助函数 ──
write_case_cfg() {  # NX NY NT
  cat > "$CASE_CFG" <<EOF
NX = $1
NY = $2
NT = $3
A = 10.0
DT = 0.01
C0 = 0.1
USE_FIXED_DOMAIN = 0
HALO = 1
ENERGY_REPORT_INTERVAL = 10000
EOF
}

write_hw_cfg() {  # N_GROUPS N_WORKERS GROUP_DECOMP
  cat > "$HW_CFG" <<EOF
N_GROUPS = $1
N_WORKERS = $2
GROUP_DECOMP = $3
NCorePClu = 5
NCluPNode = 2
NCorePGrp = 4
ManageCoreId = 4
EOF
}

# ── CSV 表头 ──
echo "program,grid,groups,workers,decomp,imbalance,run,elapsed_s,throughput_Mps,energy_init,energy_final,l2_init,l2_final,max_amp,max_amp_pos" > "$SUMMARY"

# ── 主循环 ──
total_tests=0
for grid in "${GRID_SIZES[@]}"; do
  read nx ny nt <<< "$grid"
  for ng in "${GROUP_COUNTS[@]}"; do
    for nw in "${WORKER_COUNTS[@]}"; do
      for dc in "${DECOMPS[@]}"; do
        for imb in "${IMBALANCE[@]}"; do
          write_case_cfg "$nx" "$ny" "$nt"
          write_hw_cfg "$ng" "$nw" "$dc"

          for prog in "${PROGRAMS[@]}"; do
            for r in $(seq 1 $REPEAT); do
              total_tests=$((total_tests + 1))
              tag="${prog}_${nx}x${ny}_g${ng}w${nw}_d${dc}_i${imb}_r${r}"
              logfile="$LOG_DIR/${tag}.log"

              echo "[$total_tests] $tag"

              export WAVE_IMBALANCE_PCT="$imb"
              mpirun -np 1 "./$prog" > "$logfile" 2>&1 || {
                echo "  [FAIL] exit=$?" | tee -a "$logfile"
                continue
              }

              # ── 提取指标 ──
              elapsed=$(grep "Simulation completed" "$logfile" | awk '{for(i=1;i<=NF;i++) if($i~/^[0-9]/){print $i; exit}}')
              throughput=$(grep "Throughput" "$logfile" | awk '{for(i=1;i<=NF;i++) if($i~/^[0-9]/){print $i; exit}}')
              init_e=$(grep "Initial:.*E=" "$logfile" | head -1 | sed 's/.*E=\([0-9.]*\).*/\1/')
              final_e=$(grep "E=" "$logfile" | tail -1 | sed 's/.*E=\([0-9.]*\).*/\1/')
              init_l2=$(grep "Initial:.*L2=" "$logfile" | head -1 | sed 's/.*L2=\([0-9.]*\).*/\1/')
              final_l2=$(grep "L2=" "$logfile" | tail -1 | sed 's/.*L2=\([0-9.]*\).*/\1/')
              max_amp=$(grep "Initial:.*max|u|" "$logfile" | head -1 | sed 's/.*max|u|=\([0-9.]*\)@(\([0-9]*\),\([0-9]*\)).*/\1/')
              max_pos=$(grep "Initial:.*max|u|" "$logfile" | head -1 | sed 's/.*max|u|=[0-9.]*@(\([0-9]*\),\([0-9]*\))/\1,\2/')

              decomp_name="Y_ONLY"; [ "$dc" = "1" ] && decomp_name="XY_2D"
              echo "$prog,${nx}x${ny},$ng,$nw,$decomp_name,$imb,$r,$elapsed,$throughput,$init_e,$final_e,$init_l2,$final_l2,$max_amp,($max_pos)" >> "$SUMMARY"
            done
          done
        done
      done
    done
  done
done

# ── 恢复原始配置 ──
cp "$LOG_DIR/case.cfg.bak" "$CASE_CFG"
cp "$LOG_DIR/hardware.cfg.bak" "$HW_CFG"

# ── 生成汇总表 ──
echo ""
echo "============================================="
echo "  Benchmark complete: $total_tests tests"
echo "  Logs: $LOG_DIR"
echo "  Summary: $SUMMARY"
echo "============================================="
echo ""

# 按程序、网格、分组、worker数的均值汇总
awk -F, 'NR==1{print; next} {
  key=$1","$2","$3","$4","$5","$6
  cnt[key]++; et[key]+=$8; tp[key]+=$9
}
END{
  print "\n--- Averages ---"
  printf "%-25s %10s %6s %6s %7s %8s %8s\n", "program", "grid", "grp", "wrk", "decomp", "time(s)", "Mp/s"
  for(k in cnt){
    split(k,a,",")
    printf "%-25s %10s %6s %6s %7s %8.2f %8.1f\n", a[1],a[2],a[3],a[4],a[5], et[k]/cnt[k], tp[k]/cnt[k]
  }
}' "$SUMMARY"

echo ""
echo "Logs: $LOG_DIR"
echo "Summary: $SUMMARY"
