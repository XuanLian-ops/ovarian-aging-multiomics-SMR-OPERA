#=====================================================================#
# OPERA multi-omics integration
# This script runs OPERA to integrate multiple QTL datasets.
#=====================================================================#
BESD_FLIST="/path/to/mylist.txt"
GWAS="/path/to/mygwas.ma"
MBFILE="/path/to/mybdatalist.txt"
BFILE="/path/to/mydata"
OUT_PREFIX="/path/to/myopera"

for f in "$BESD_FLIST" "$GWAS" "$MBFILE"; do
    if [[ ! -f "$f" ]]; then
        echo "Error: File $f not found!"
        exit 1
    fi
done

if [[ ! -f "${BFILE}.bed" ]] || [[ ! -f "${BFILE}.bim" ]] || [[ ! -f "${BFILE}.fam" ]]; then
    echo "Error: Reference genotype files ($BFILE) missing .bed/.bim/.fam"
    exit 1
fi

out_dir=$(dirname "$OUT_PREFIX")
mkdir -p "$out_dir"

echo "========================================="
echo "Estimating prior parameters"
echo "========================================="

$OPERA \
    --besd-flist "$BESD_FLIST" \
    --gwas-summary "$GWAS" \
    --mbfile "$MBFILE" \
    --estimate-pi \
    --thread-num "$THREADS" \
    --out "${OUT_PREFIX}_prior"

if [[ $? -ne 0 ]]; then
    echo "Error: Prior estimation failed."
    exit 1
fi

# ------------------------------ #
# Multi-omics integration
# ------------------------------ #
echo "========================================="
echo "Running multi-omics integration"
echo "========================================="

$OPERA \
    --besd-flist "$BESD_FLIST" \
    --gwas-summary "$GWAS" \
    --bfile "$BFILE" \
    --sample-overlap \
    --estimate-pi \
    --thread-num "$THREADS" \
    --out "$OUT_PREFIX"

if [[ $? -ne 0 ]]; then
    echo "Error: Multi-omics integration failed."
    exit 1
fi

echo "========================================="
echo "OPERA analysis completed successfully."
echo "Results saved with prefix: $OUT_PREFIX"
echo "========================================="