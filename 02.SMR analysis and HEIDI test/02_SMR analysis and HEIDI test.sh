#=====================================================================#
# SMR analysis and HEIDI test for multiple QTL datasets
# This script runs SMR for pQTL, eQTL, sQTL, and mQTL.
#=====================================================================#

SMR="smr-1.3.1"
BFILE="/path/to/reference/g1000_eur"
GWAS="/path/to/gwas/gwas.txt"
QTL_BESD_DIR="/path/to/data"
OUT_DIR="/path/to/smr_result"

declare -A qtls=(
    [pQTL]="pqtl_besd"
    [eQTL]="eqtl_besd"
    [sQTL]="sqtl_besd"
    [mQTL]="mqtl_besd"
)

for qtl in "${!qtls[@]}"; do
    echo "========================================="
    echo "Processing $qtl ..."
    echo "========================================="

    $SMR \
        --bfile "$BFILE" \
        --gwas-summary "$GWAS" \
        --beqtl-summary "$QTL_BESD_DIR/${qtls[$qtl]}" \
        --diff-freq-prop 0.2 \
        --thread-num "$THREADS" \
        --out "$OUT_DIR/${qtl}_smr"

    echo "$qtl finished."
    echo ""
done

echo "All SMR analyses completed."

# Post-processing notes (significance filtering)
# The SMR output file (*.smr) contains the following columns:
#   ProbeID, SMR_p, HEIDI_p
#
# To obtain significant associations after Bonferroni correction:
#   1. Determine the number of probes tested in each QTL dataset (N_probes).
#   2. Filter rows where:
#        SMR_p < (0.05 / N_probes)   (Bonferroni-corrected threshold)
#        HEIDI_p > 0.01 
