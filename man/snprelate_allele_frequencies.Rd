% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/snprelate_ld.R
\name{snprelate_allele_frequencies}
\alias{snprelate_allele_frequencies}
\title{Compute allele frequencie and snp missing rate}
\usage{
snprelate_allele_frequencies(
  gdata,
  snps_idx = NULL,
  scans_idx = NULL,
  quiet = FALSE
)
}
\arguments{
\item{gdata}{A GenotypeData object}

\item{snps_idx}{Vector of snps indices}

\item{scans_idx}{Vector of scans indices}

\item{quiet}{Whether to be quiet}
}
\value{
A data frame of snps_idx, snps_ids, allele1, allele2, maf, missing
  where allele1 and allele2 are the rates of the alleles, and maf the minimum
  of the 2. Missing is the missing rate. N.B: the allele rates are computed
  on the non missing genotypes, i.e. their sum equals 1.
}
\description{
Wrapper over SNPRelate::snpgdsSNPRateFreq
}
