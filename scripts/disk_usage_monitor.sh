#!/usr/bin/env bash
# ============================================================================
# Disk Usage Monitor
#
# 파이프라인 실행 중 디스크 사용량을 모니터링합니다.
# 별도 터미널에서 실행: watch -n 30 bash scripts/disk_usage_monitor.sh
# ============================================================================

PROJECT_DIR="/media/jaeil/ff387598-88d2-44f2-b68e-cf799d26fb68/lncRNA-GBM-TMZ"
MOUNT="/media/jaeil/ff387598-88d2-44f2-b68e-cf799d26fb68"

echo "═══════════════════════════════════════════"
echo "  Disk Usage Monitor — $(date '+%Y-%m-%d %H:%M:%S')"
echo "═══════════════════════════════════════════"

# Overall disk
echo ""
df -h "${MOUNT}" | tail -1 | awk '{printf "  Total: %s | Used: %s | Avail: %s | %s\n", $2, $3, $4, $5}'

# Per-directory breakdown
echo ""
echo "  Directory breakdown:"
for dir in data/raw data/processed/trimmed data/processed/aligned data/processed/counts \
           data/reference results; do
    full_path="${PROJECT_DIR}/${dir}"
    if [[ -d "${full_path}" ]]; then
        size=$(du -sh "${full_path}" 2>/dev/null | cut -f1)
        n_files=$(find "${full_path}" -type f 2>/dev/null | wc -l)
        printf "    %-40s %8s  (%d files)\n" "${dir}/" "${size}" "${n_files}"
    fi
done

# Warning threshold
avail_gb=$(df --output=avail "${MOUNT}" | tail -1 | awk '{print int($1/1048576)}')
if [[ "${avail_gb}" -lt 100 ]]; then
    echo ""
    echo "  ⚠ WARNING: Less than 100 GB remaining!"
    echo "  Consider running with --keep-bam=false or cleaning up raw data."
fi

echo ""
