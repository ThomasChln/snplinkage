
#' snprelate_qc
#'
#' Quality control using SNPRelate functions.
#'
#' @param gdata       Genotype data object
#' @param samples_nas NA threshold for samples, default 3 pct
#' @param ibs         Samples identity by state threshold, default 99 pct
#' @param keep_ids    Samples ids to keep even if IBS is higher than threshold.
#'                    Used for monozygotic twins.
#' @param snps_nas    NA threshold for SNPs, default 1 pct
#' @param maf         Minor allele frequency threshold, default 5 pct
#' @param tagsnp      TagSNP r2 correlation threshold, default 0.8
#' @param n_cores     Number of cores
#' @return List of gdata, Genotype data object, and df_qc, QC info data frame
#'
#' @examples
#' library(snplinkage)
#' gds_path <- save_hgdp_as_gds()
#' gdata <- load_gds_as_genotype_data(gds_path)
#' qc <- snprelate_qc(gdata, tagsnp = .99)
#'
#' @export
snprelate_qc <- function(gdata, samples_nas = 0.03, ibs = 0.99, keep_ids = NULL,
  snps_nas = 0.01, maf = 0.05, tagsnp = 0.8, n_cores = 1) {

  stopifnot(inherits(gdata, 'GenotypeData'))
  samples = getScanID(gdata)
  snps = getSnpID(gdata)
  df_qc <- .rbind_qc(NULL, 'Raw', NA, samples, snps)
  gds <- request_snpgds_file(gdata)$snpgds

  # samples NAs
  nas_rates <- suppressMessages(
    SNPRelate::snpgdsSampMissRate(gds, samples, snps))
  samples <- samples[nas_rates <= samples_nas]
  df_qc <- .rbind_qc(df_qc, 'Samples NAs', samples_nas, samples, snps)

  # samples IBS
  m_ibs <- suppressMessages(SNPRelate::snpgdsIBS(gds, samples, snps,
      num.thread = n_cores, verbose = FALSE))$ibs
  m_ibs[lower.tri(m_ibs, TRUE)] <- NA
  max_ibs <- apply(m_ibs[, -1], 2, max, na.rm = TRUE)
  samples_ibs <- samples[c(TRUE, max_ibs < ibs)]

  # keep ids
  if (!is.null(keep_ids) &&
    any(keep_ids %in% samples)) {
    samples_ibs <- union(samples_ibs, keep_ids[keep_ids %in% samples])
  }
  samples <- samples_ibs
  df_qc <- .rbind_qc(df_qc, 'Identity by state - Twins', ibs, samples, snps)

  # SNPs NAs
  snps <- unlist(suppressMessages(SNPRelate::snpgdsSelectSNP(gds,
        samples, snps, missing.rate = snps_nas, verbose = FALSE)))
  df_qc <- .rbind_qc(df_qc, 'SNPs NAs', snps_nas, samples, snps)

  # MAF
  snps <- unlist(suppressMessages(SNPRelate::snpgdsSelectSNP(gds,
       samples, snps, maf = maf, verbose = FALSE)))
  df_qc <- .rbind_qc(df_qc, 'MAF', maf, samples, snps)

  # TagSNP
  if (!is.na(tagsnp)) {
    snps <- unlist(suppressMessages(SNPRelate::snpgdsLDpruning(gds,
          samples, snps, ld.threshold = sqrt(tagsnp), method = 'r',
          verbose = FALSE)))
    df_qc <- .rbind_qc(df_qc, 'TagSNP', tagsnp, samples, snps)
  }
  samples = match(samples, getScanID(gdata))
  snps = match(snps, getSnpID(gdata))
  gdata <- genotype_data_subset(gdata, snps, samples)

  list(gdata = gdata, df_info = df_qc)
}

.rbind_qc <- function(df_qc, Step, Parameter, samples, snps) {
  df_qc_new_line <- data.frame(Step, Parameter,
    Samples = length(samples), SNPs = length(snps))
  rbind(df_qc, df_qc_new_line)
}

