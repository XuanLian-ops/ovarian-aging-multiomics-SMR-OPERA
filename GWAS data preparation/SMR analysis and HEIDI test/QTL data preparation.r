#=====================================================================#
# Prepare QTL summary statistics (pQTL, eQTL, sQTL, mQTL) for SMR
# This script reads four CSV files (one per QTL type) and generates:
#   - Per-probe .esd files 
#   - One .flist file per QTL type 
#=====================================================================#

# ------------------------------ #
# 0. Configuration (EDIT THESE PATHS)
# ------------------------------ #
# Input CSV files
csv_files <- list(
  pQTL = "/path/to/pqtl_results.csv",
  eQTL = "/path/to/eqtl_results.csv",
  sQTL = "/path/to/sqtl_results.csv",
  mQTL = "/path/to/mqtl_results.csv"
)

# Output base directory
out_base <- "/path/to/output"


library(data.table) 
create_esd_files <- function(qtl_df, out_dir) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # Split by ProbeID
  probe_list <- split(qtl_df, qtl_df$ProbeID)
  
  # Prepare a list to collect probe info for .flist
  probe_info <- list()
  
  for (probe in names(probe_list)) {
    # Extract data for this probe and ensure correct column order
    esd_data <- probe_list[[probe]][, c("Chr", "SNP", "Bp", "A1", "A2", "Freq", "Beta", "SE", "P")]
    
    # Define output file path
    esd_file <- file.path(out_dir, paste0(probe, ".esd"))
    
    # Write .esd file (space-separated, no header)
    fwrite(esd_data, file = esd_file, sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    # Collect basic probe info (Chr and median Bp as approximate ProbeBp)
    probe_info[[probe]] <- data.frame(
      Chr = unique(esd_data$Chr)[1],   # assume all SNPs for a probe are on same chromosome
      ProbeID = probe,
      ProbeBp = median(esd_data$Bp),   # median SNP position as placeholder for probe location
      stringsAsFactors = FALSE
    )
  }
  
  # Combine probe info into a data.frame
  probe_info_df <- do.call(rbind, probe_info)
  rownames(probe_info_df) <- NULL
  return(probe_info_df)
}

create_flist <- function(probe_info, out_dir, anno = NULL) {
  # Merge with annotations if provided
  if (!is.null(anno)) {
    probe_info <- merge(probe_info, anno, by = "ProbeID", all.x = TRUE)
  } else {
    probe_info$Gene <- ""
    probe_info$Orientation <- "0"
  }
  
  # Build full path to .esd files
  probe_info$PathOfEsd <- file.path(normalizePath(out_dir), paste0(probe_info$ProbeID, ".esd"))
  
  # Create flist data with correct column order
  flist <- data.frame(
    Chr = probe_info$Chr,
    ProbeID = probe_info$ProbeID,
    GeneticDistance = 0,                # fixed to 0
    ProbeBp = probe_info$ProbeBp,
    Gene = probe_info$Gene,
    Orientation = probe_info$Orientation,
    PathOfEsd = probe_info$PathOfEsd,
    stringsAsFactors = FALSE
  )
  
  return(flist)
}


for (qtl_type in names(csv_files)) {
  cat("\nProcessing", qtl_type, "...\n")
  
  # Check input file exists
  if (!file.exists(csv_files[[qtl_type]])) {
    warning("File not found: ", csv_files[[qtl_type]], ". Skipping.")
    next
  }

  qtl_df <- fread(csv_files[[qtl_type]], data.table = FALSE)

  required_cols <- c("ProbeID", "Chr", "SNP", "Bp", "A1", "A2", "Freq", "Beta", "SE", "P")
  missing_cols <- setdiff(required_cols, colnames(qtl_df))
  if (length(missing_cols) > 0) {
    stop("Missing columns in ", qtl_type, " file: ", paste(missing_cols, collapse = ", "))
  }
  qtl_df$Chr <- gsub("chr", "", qtl_df$Chr, ignore.case = TRUE)
  qtl_df$Chr[qtl_df$Chr == "X"] <- 23
  qtl_df$Chr[qtl_df$Chr == "Y"] <- 24
  qtl_df$Chr <- as.integer(qtl_df$Chr)   # convert to integer
  # Create output subdirectory for this QTL type
  out_dir_qtl <- file.path(out_base, qtl_type)
  # Generate per-probe .esd files and collect probe info
  probe_info <- create_esd_files(qtl_df, out_dir_qtl)
  # Create .flist data frame
  flist_df <- create_flist(probe_info, out_dir_qtl, anno = probe_anno)
  # Write .flist file
  flist_file <- file.path(out_base, paste0(qtl_type, ".flist"))
  fwrite(flist_df, file = flist_file, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  cat("  Done. Generated", nrow(flist_df), "probes for", qtl_type, "\n")
  cat("  .flist saved to:", flist_file, "\n")
}
