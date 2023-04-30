
#' Compute Chi-squared p-values on a Genotype data object
#'
#' @param gdata           Genotype data object
#' @param snp_idxs        SNPs indexes
#' @param response_column Response column in gdata scans annotations data frame
#' @param response_value  Response value. The response vector will be a logical,
#'                        true if equal to the value, false otherwise.
#' @param threshold       Keep only associations greater than the threshold
#' @param ...             Passed to chisq_pvalues
#' @return SNPs annotation data frame, chi-squared p-values in column pvalues
#' @export
chisq_pvalues_gdata <- function(gdata, snp_idxs, response_column = 'region',
  response_value = 'Europe', threshold = 2, ...) {

  m_geno <- t(fetch_genotypes(gdata, snp_idxs))
  df_snp <- gdata@snpAnnot@data[snp_idxs, ]
  df_snp$pvalues <- chisq_pvalues(m_geno,
    gdata@scanAnnot@data[[response_column]] == response_value, ...)

  if (threshold > 0) df_snp %<>% subset(pvalues > threshold)

  df_snp
}

#' Compute Chi-squared p-values
#'
#' @param m_data           Data matrix of observations by variables
#' @param response         Response vector of length the number of observations
#' @param adjust_method    Multiple testing p-value adjustment method.
#'                         Passed to stats::p.adjust. 'fdr' by default.
#' @param mlog10_transform Logical, transform p-values by minus log10.
#'                         True by default.  
#' @param n_cores          Number of cores
#' @param ...              Passed to stats::chisq.test
#' @return Chi-squared p-values
#' @export
chisq_pvalues <- function(m_data, response, adjust_method = 'fdr',
  mlog10_transform = TRUE, n_cores = 1, ...) {

  pvals <- parallel_apply(m_data, .chisq_pvalues, n_cores, response, ...)
  if (!is.null(adjust_method)) pvals <- stats::p.adjust(pvals, adjust_method)
  if (mlog10_transform) pvals <- -log10(pvals)

  pvals
} 

.chisq_pvalues <- function(m_data, response, ...) {
  apply(m_data, 2, function(v_data) {
      catch_warnings(stats::chisq.test(v_data, response, ...)$p.value)$result
    })
}

#' Separate a matrix in a list of matrices of length the number of cores and
#' apply a function on the columns in parallel
#'
#' @param m_data    Data matrix
#' @param apply_fun Function to apply
#' @param n_cores   Number of cores
#' @param ...       Passed to apply_fun
#' @return apply_fun return
#' @export
parallel_apply <- function(m_data, apply_fun, n_cores = 1, ...) {

  if (n_cores == 1) return(apply_fun(m_data, ...))

  column_idxs <- floor(seq(1, ncol(m_data) + 1, length = n_cores + 1))
  df_idxs <- rbind.data.frame(utils::head(column_idxs, -1), column_idxs[-1] - 1)
  m_data <- lapply(df_idxs, function(idxs) m_data[, seq(idxs[1], idxs[2])])
  result <- parallel::mclapply(m_data, apply_fun, ..., mc.cores = n_cores)

  unlist(result, FALSE)
}
