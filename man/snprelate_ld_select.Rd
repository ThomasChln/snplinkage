% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/snprelate_ld.R
\name{snprelate_ld_select}
\alias{snprelate_ld_select}
\title{Wrapper for snpgdsLDpruning to select Tag SNPs}
\usage{
snprelate_ld_select(
  gdata,
  window_length = 500L,
  min_r2,
  window_size = NA,
  snps_idx = NULL,
  scans_idx = NULL,
  remove.monosnp = FALSE,
  autosome.only = FALSE,
  method = "r",
  threads = 1,
  quiet = FALSE,
  ...
)
}
\arguments{
\item{gdata}{A GenotypeData object}

\item{window_length}{Max length in kb of the window}

\item{min_r2}{Minimum r2 value to report}

\item{window_size}{Max number of SNPs in LD window}

\item{snps_idx}{Indices of snps to use}

\item{scans_idx}{Indices of scans to use}

\item{remove.monosnp}{if TRUE, remove monomorphic SNPs}

\item{autosome.only}{if \code{TRUE}, use autosomal SNPs only; if it is a
        numeric or character value, keep SNPs according to the specified
        chromosome}

\item{method}{"composite", "r", "dprime", "corr", see details}

\item{threads}{The number of threads to use, currently ignored}

\item{quiet}{Whether to be quiet}

\item{...}{Forwarded to SNPRelate::snpgdsLDpruning}
}
\value{
A list of SNP IDs stratified by chromosomes.
}
\description{
The tagged snp set is (by sliding window) representative and strongly not
redundant.
}
