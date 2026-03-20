#/bin/bash
#===================================================================#
#   GTEx eQTL Regional Plot Data Preparation (OPERA Analysis)
#   Gene: DHX58 (ENSG00000108771)
#   This script extracts eQTL data for the target region from
#   multiple eQTL datasets, merges them, and prepares input for SMR.
#===================================================================#

# ------------------------------ #
#   Set environment variables
# ------------------------------ #
SMR="/path/to/smr-1.3.1-linux-x86_64/smr-1.3.1"          # Update this path
pQTL="/path/to/pQTL_sig"                        # Update
eQTL="/path/to/eQTL_sig"
scriptQTL="/path/to/sQTL_sig"
gwas="/path/to/gwas_file.txt"                             # GWAS summary stats
reference="/path/to/ukbV3_eur_unrel_n10k_chr17"            # Reference PLINK files
genelist="/path/to/refseq_isoform_SMR.txt"                # Gene list for SMR
dir="/path/to/output_directory"                            # Working directory

cd $dir

# ------------------------------ #
#   Parameters for region
# ------------------------------ #
chr=17
center_bp=40134782
window=1000000
probe_id="ENSG00000108771"
start=$((center_bp - window))
end=$((center_bp + window))

# ------------------------------ #
#   1. Generate probe lists for each eQTL dataset
# ------------------------------ #
# Read .epi files to get probe IDs in the region
Rscript -e "
epi_pQTL  <- read.delim('${pQTL}.epi', header=FALSE)
epi_eQTL  <- read.table('${eQTL}.epi', header=FALSE)
epi_sQTL  <- read.delim('${sQTL}.epi', header=FALSE)

idx_pQTL   <- which(epi_pQTL$V1 == $chr & epi_pQTL$V4 >= $start & epi_pQTL$V4 <= $end)
idx_eQTL   <- which(epi_eQTL$V1 == $chr & epi_eQTL$V4 >= $start & epi_eQTL$V4 <= $end)
idx_sQTL   <- which(epi_sQTL$V1 == $chr & epi_sQTL$V4 >= $start & epi_sQTL$V4 <= $end)


write.table(epi_pQTL[idx_pQTL, 2], 'pQTL_${probe_id}.list', col.names=F, row.names=F, quote=F)
write.table(epi_eQTL[idx_pQTL, 2], 'eQTL_${probe_id}.list', col.names=F, row.names=F, quote=F)
write.table(epi_sQTL[idx_pQTL, 2], 'sQTL_${probe_id}.list', col.names=F, row.names=F, quote=F)
"

# ------------------------------ #
#   2. Extract BESD files for the selected probes
# ------------------------------ #
$SMR --beqtl-summary $pQTL --extract-probe pQTL_${probe_id}.list --chr $chr --make-besd --out pQTL_${probe_id}
$SMR --beqtl-summary $eQTL --extract-probe eQTL_${probe_id}.list --chr $chr --make-besd --out eQTL_${probe_id}
$SMR --beqtl-summary $sQTL --extract-probe sQTL_${probe_id}.list --chr $chr --make-besd --out sQTL_${probe_id}


# ------------------------------ #
#   3. Modify probe names in .epi files to avoid duplicates
# ------------------------------ #
Rscript -e "
epi_pQTL <- read.delim('pQTL_${probe_id}.epi', header=F)
epi_eQTL <- read.delim('eQTL_${probe_id}.epi', header=F)
epi_sQTL <- read.delim('sQTL_${probe_id}.epi', header=F)

epi_pQTL\$V2  <- paste0('pQTL_', epi_pQTL\$V2)
epi_eQTL\$V2  <- paste0('eQTL_', epi_eQTL\$V2)
epi_sQTL\$V2  <- paste0('sQTL_', epi_sQTL\$V2)

write.table(epi_pQTL, 'pQTL_${probe_id}.epi', row.names=F, col.names=F, sep='\t', quote=F)
write.table(epi_eQTL, 'eQTL_${probe_id}.epi', row.names=F, col.names=F, sep='\t', quote=F)
write.table(epi_sQTL, 'sQTL_${probe_id}.epi', row.names=F, col.names=F, sep='\t', quote=F)
"

# ------------------------------ #
#   4. Update gene positions in .esi file (if needed)
#    Example: copy positions from hypothalamus to eQTL
# ------------------------------ #
Rscript -e "
library(data.table)
esi1 <- fread('pQTL_${probe_id}.esi', header=F, data.table=F)
esi2 <- fread('eQTL_${probe_id}.esi', header=F, data.table=F)
idx <- match(esi1\$V2, esi2\$V2, nomatch=0)
esi1\$V4[which(idx!=0)] <- esi2\$V4[idx]
write.table(esi1, 'pQTL_${probe_id}.esi', row.names=F, col.names=F, sep='\t', quote=F)
"

# ------------------------------ #
#   5. Create a list of BESD files and merge them
# ------------------------------ #
echo "pQTL_${probe_id}" > merged_${probe_id}.list
echo "eQTL_${probe_id}" > merged_${probe_id}.list
echo "sQTL_${probe_id}" > merged_${probe_id}.list


$SMR --besd-flist merged_${probe_id}.list --make-besd --out merged_${probe_id}

# ------------------------------ #
#   6. Run SMR to generate plot data
# ------------------------------ #
$SMR \
    --beqtl-summary merged_${probe_id} \
    --gwas-summary $gwas \
    --bfile $reference \
    --plot \
    --diff-freq-prop 0.2 \
    --probe-wind 2000 \
    --gene-list $genelist \
    --out menopause_${probe_id}

# The SMR --plot option creates a .txt file named menopause_${probe_id}.txt
# which will be used in the R script.