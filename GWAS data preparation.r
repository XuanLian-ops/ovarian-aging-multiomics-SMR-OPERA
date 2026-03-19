#GWAS summary statistics for age at natural menopause (ANM) were obtained from public sources.
#### Processing steps:
#1. Remove SNPs with missing values  
#2. Filter SNPs with minor allele frequency (MAF) < 0.01  
#3. Format GWAS summary statistics with required columns: SNP, A1, A2, freq, beta, se, p

library(data.table)
library(dplyr)

input_file <- "raw_gwas.txt"
output_file <- "gwas.txt"

gwas <- fread(input_file)

gwas <- gwas %>%
  filter(
    !is.na(SNP),
    !is.na(A1),
    !is.na(A2),
    !is.na(freq),
    !is.na(beta),
    !is.na(se),
    !is.na(p)
  )

gwas <- gwas %>%
  filter(freq >= 0.01 & freq <= 0.99)

gwas <- gwas %>%
  mutate(
    A1 = toupper(A1),
    A2 = toupper(A2)
  )

gwas_formatted <- gwas %>%
  transmute(
    SNP  = SNP,
    A1   = A1,
    A2   = A2,
    freq = freq,
    b    = beta,
    se   = se,
    p    = p,
    n    = ifelse("n" %in% colnames(gwas), n, NA)
  )

#  Output File
fwrite(
  gwas_formatted,
  file = output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)