#' Wrapper for snpgdsLDpruning to select Tag SNPs
#'
#' The tagged snp set is (by sliding window) representative and strongly not
#' redundant.
#'
#' @param gdata         A GenotypeData object
#' @param window_size   Max number of SNPs in LD window
#' @param window_length Max length in kb of the window
#' @param min_r2        Minimum r2 value to report
#' @param snps_idx      Indices of snps to use
#' @param scans_idx     Indices of scans to use
#' @param threads       The number of threads to use, currently ignored
#' @param quiet         Whether to be quiet
#' @inheritParams SNPRelate::snpgdsLDpruning
#' @param ...           Forwarded to SNPRelate::snpgdsLDpruning
#'
#' @return A list of SNP IDs stratified by chromosomes.
#' @export
snprelate_ld_select <- function(gdata,
  window_length = 500L,
  min_r2,
  window_size = NA,
  snps_idx = NULL,
  scans_idx = NULL,
  remove.monosnp = FALSE,
  autosome.only = FALSE,
  method = 'r',
  threads = 1,
  quiet = FALSE,
  ...
) {
  info <- if (quiet) function(...) {} else function(...) message(sprintf(...))

  info('managing input')
  gdsinfo <- request_snpgds_file(gdata, snps_idx, scans_idx)
  gds <- gdsinfo$snpgds
  if (gdsinfo$new_file) {
    on.exit({closefn.gds(gds); unlink(gds$filename)}, TRUE)
  }

  info('SNPRelate::snpgdsLDpruning')

  # there is currently a bug in snpgdsLDpruning with huge values
  bp <- min(1e9, window_length * 1000)

  suppressMessages(
    res <- SNPRelate::snpgdsLDpruning(gds, gdsinfo$scan_ids, gdsinfo$snp_ids,
      slide.max.bp = bp,
      slide.max.n = window_size,
      ld.threshold = sqrt(min_r2),
      remove.monosnp = remove.monosnp,
      autosome.only = autosome.only,
      method = method, num.thread = threads,
      verbose = !quiet)
  )

  res
}

#' Wrapper for snpgdsLDMat to compute r2
#'
#' @param gdata       A GenotypeData object
#' @param window_size Max number of SNPs in LD window, 0 for no window
#' @param min_r2      Minimum r2 value to report
#' @param snps_idx    Indices of snps to use
#' @param scans_idx   Indices of scans to use
#' @param threads     The number of threads to use
#' @param quiet       Whether to be quiet
#'
#' @return A data frame with columns SNP_A, SNP_B, R2 for r2 >= min_r2
#' @export
snprelate_ld <- function(gdata,
  window_size = 0,
  min_r2 = 0,
  snps_idx = NULL,
  scans_idx = NULL,
  threads = 1,
  quiet = FALSE
) {

  # lazy-loading, needs to be in search path
  Library <- library; Library('SNPRelate')

  info <- if (quiet) function(...) {} else function(...) message(sprintf(...))

  info('managing input')
  gdsinfo <- request_snpgds_file(gdata, snps_idx, scans_idx)
  gds <- gdsinfo$snpgds
  if (gdsinfo$new_file) {
    on.exit({closefn.gds(gds); unlink(gds$filename)}, TRUE)
  }

  info('SNPRelate::snpgdsLDMat')
  if (window_size <= 0) window_size <- length(snps_idx)
  suppressMessages(
    mat <- SNPRelate::snpgdsLDMat(gds, gdsinfo$scan_ids, gdsinfo$snp_ids,
      slide = window_size - 1, method = 'r', num.thread = threads,
      verbose = !quiet)
  )
  info('converting output')
  if (is.null(snps_idx)) snps_idx <- seq_len(nsnp(gdata))
  .convert_snpgdsLDMat_output(mat$LD^2, min_r2, getSnpID(gdata, snps_idx))
}

# the matrix output by snpgdsLDMat is n*m such that m[i, j] = LD(j, j+i)
.convert_snpgdsLDMat_output <- function(ld, min_r2, snp_ids) {

  nb_rows <- nrow(ld)
  nb_cols <- ncol(ld)

  goods <- which(ld >= min_r2, arr.ind = TRUE)
  # remove extra empty cells
  goods <- goods[rowSums(goods) <= nb_cols, , drop = FALSE]

  ind1 <- goods[, 2]
  ind2 <- rowSums(goods)

  ids1 <- as.character(snp_ids[ind1])
  ids2 <- as.character(snp_ids[ind2])
  r2 <- ld[goods]

  data.frame(SNP_A = ids1, SNP_B = ids2, R2 = r2, stringsAsFactors = FALSE)
}

#' Compute allele frequencie and snp missing rate
#'
#' Wrapper over SNPRelate::snpgdsSNPRateFreq
#'
#' @param gdata     A GenotypeData object
#' @param snps_idx  Vector of snps indices
#' @param scans_idx Vector of scans indices
#' @param quiet     Whether to be quiet
#'
#' @return A data frame of snps_idx, snps_ids, allele1, allele2, maf, missing
#'   where allele1 and allele2 are the rates of the alleles, and maf the minimum
#'   of the 2. Missing is the missing rate. N.B: the allele rates are computed
#'   on the non missing genotypes, i.e. their sum equals 1.
#' @export
snprelate_allele_frequencies <- function(gdata,
  snps_idx = NULL,
  scans_idx = NULL,
  quiet = FALSE
) {
  info <- if (quiet) function(...) {} else function(...) message(sprintf(...))

  # lazy-loading, needs to be in search path
  Library <- library; Library('SNPRelate')

  info('managing input')
  gi <- request_snpgds_file(gdata, snps_idx, scans_idx)
  gds <- gi$snpgds
  if (gi$new_file) {
    on.exit({closefn.gds(gds); unlink(gds$filename)}, TRUE)
  }

  # get rid of:Hint: it is suggested to call `snpgdsOpen' to open a GDS file
  # although we are actually using it (cf snp_gds_open)
  suppressMessages(
    res <- SNPRelate::snpgdsSNPRateFreq(gds, sample.id = gi$scan_ids,
      snp.id = gi$snp_ids)
  )
  if (is.null(snps_idx)) snps_idx <- seq_len(nsnp(gdata))
  snps_ids <- as.character(getSnpID(gdata, snps_idx))
  df <- data.frame(snps_idx = snps_idx, snps_ids = snps_ids,
    stringsAsFactors = FALSE)
  res_df <- as.data.frame(res)

  a1 <- res_df$AlleleFreq
  a2 <- 1 - a1

  df$allele1 <- a1
  df$allele2 <- a2
  df$maf <- res_df$MinorFreq
  df$missing <- res_df$MissingRate

  df
}

