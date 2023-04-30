
actg_tsv_to_gdata <- function(geno_path, scans_path,
  scans_col_map = 'scan_id', snps_path,
  snps_col_map = c('probe_id', 'chromosome', 'position'),
  na_encoding = '--') {

  require_defaults(snps_col_map, 'snps_col_map', actg_tsv_to_gdata)
  require_defaults(scans_col_map, 'scans_col_map', actg_tsv_to_gdata)

  # read all, convert snps and scans
  snps <- txt_snps_to_df(snps_path, snps_col_map)
  scans <- txt_scans_to_df(scans_path, scans_col_map)

  # dont use fread as matrix conversion takes forever
  geno <- scan(geno_path, 'character', sep = '\t', quiet = TRUE)
  geno <- matrix(geno, nrow = nrow(snps) + 1, byrow = TRUE)
  dimnames(geno) <- list(geno[, 1], geno[1, ])
  geno = geno[-1, -1]

  # order genotype following SNPs and scans
  geno <- geno[match(snps$probe_id, rownames(geno)), ]
  scan_order_idxs <- match(scans$scan_id, colnames(geno))
  scans <- scans[!is.na(scan_order_idxs), ]
  scans <- cbind(scanID = seq_along(scans[[1]]), scans)
  geno <- geno[, stats::na.omit(scan_order_idxs)]

  # convert genotype and write
  geno <- t(apply(geno, 1, actg_to_numeric, na_encoding))

  build_gwastools(geno, scans, snps)
}

require_defaults <- function(param, param_name, fun) {
  required_col <- as.character(formals(fun)[[param_name]])[-1]
  missing <- required_col[!(required_col %in% param)]
  if (!length(missing)) {
    paste0('Missing fields in ', param_name, ': ',
      paste0(missing, collapse = ', '))
  }
}

# Converts a character actg matrix to a numeric 0,1,2 matrix
actg_to_numeric <- function(snp, na_encoding = NA) {
  f_snp <- factor(snp, exclude = na_encoding, nmax = 3)
  i_snp <- as.integer(f_snp) - 1L
  if (nlevels(f_snp) == 2) {
    lvls <- levels(f_snp)
    dble_homo <- !grepl(substr(lvls[1], 1, 1), lvls[2])
    if (dble_homo) i_snp <- i_snp * 2L
  }

  i_snp
}

# read and convert to factors
txt_scans_to_df <- function(path, col_map) {
  data <- data.table::fread(path, sep = '\t', data.table = FALSE)
  data <- data[!is.na(col_map)]
  names(data) <- stats::na.omit(col_map)
  id_idx <- match('scan_id', names(data))
  data[-id_idx] <- lapply(data[-id_idx], factor)

  data
}

# read, convert to integers, and order
txt_snps_to_df <- function(path, col_map) {
  data <- data.table::fread(path, sep = '\t', data.table = FALSE)
  data <- data[!is.na(col_map)]
  names(data) <- stats::na.omit(col_map)
  data$chromosome <- as.integer(factor(data$chromosome,
      c(1:22, 'X', 'XY', 'Y', 'M'), nmax = 26, exclude = NULL))
  data$position <- as.integer(data$position)
  data <- data[order(data$chromosome, data$position), ]

  cbind(snpID = seq_along(data[[1]]), data, alleleA = NA, alleleB = NA)
}

build_gwastools <- function(geno, scans, snps) {
  stopifnot(inherits(geno, 'matrix'))

  rownames(scans) <- NULL
  rownames(snps) <- NULL

  if (!('snpID' %in% names(snps))) snps$snpID <- seq_len(nrow(snps))
  if (!('scanID' %in% names(scans))) scans$scanID <- seq_len(nrow(scans))
  geno <- MatrixGenotypeReader(geno, snps$snpID, snps$chromosome, snps$position,
    scans$scanID)
  snps <- SnpAnnotationDataFrame(snps)
  scans <- ScanAnnotationDataFrame(scans)

  GenotypeData(geno, snps, scans)
}

.gds_read_index <- function(field, gds) read.gdsn(index.gdsn(gds, field))

