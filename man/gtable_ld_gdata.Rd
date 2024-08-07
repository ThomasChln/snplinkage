% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gtable_ld.R
\name{gtable_ld_gdata}
\alias{gtable_ld_gdata}
\title{Gtable of linkage disequilibrium and positions using a GenotypeData object}
\usage{
gtable_ld_gdata(
  gdata,
  snps_idx,
  maf = NULL,
  r2 = NULL,
  diamonds = length(snps_idx) < 40,
  window = 15,
  autotitle = TRUE,
  autotitle_bp = TRUE,
  double_title = FALSE,
  ...
)
}
\arguments{
\item{gdata}{GenotypeData object returned by load_gds_as_genotype_data}

\item{snps_idx}{SNPs indexes to select}

\item{maf}{Minor allele frequency threshold(s), see description}

\item{r2}{TagSNP r2 threshold(s), see description}

\item{diamonds}{Display the values as diamonds or as points
Default is TRUE for less than 40 SNPs.}

\item{window}{Window size for snprelate_ld.
Forced to the total number of SNPs if diamonds is FALSE}

\item{autotitle}{Set title to feature selection method(s), number of SNPs
and chromosome}

\item{autotitle_bp}{Set biplot title to feature selection method(s), number
of SNPs and chromosome}

\item{double_title}{Logical, if false (default) keep only biplot title}

\item{...}{Passed to gtable_ld}
}
\value{
gtable of ggplots
}
\description{
Compute linkage disequilibrium using snprelate_ld on a set of SNP indexes and
call gtable_ld.
Two parameters are available to compute and compare minor allele frequency
filtering and TagSNP selection by displaying two LD plots with their
positions in the center.
The maf and r2 parameters are used similarly and as follows:
   - compare baseline with MAF 5%
   gtable_ld(gdata, snps_idx, maf = 0.05)
   - compare baseline with TagSNP r2 = 0.8
   gtable_ld(gdata, snps_idx, r2 = 0.8)
   - compare 5% MAF with 5% MAF and r2 = 0.8
   gtable_ld(gdata, snps_idx, maf = c(0.05, 0.05), r2 = 0.8)
   - compare MAF 5%, r2 0.8 with MAF 10%, r2 = 0.6
   gtable_ld(gdata, snps_idx, maf = c(0.05, 0.1), r2 = c(0.8, 0.6))
}
\examples{
library(snplinkage)
gds_path <- save_hgdp_as_gds()
gdata <- load_gds_as_genotype_data(gds_path)
qc <- snprelate_qc(gdata, tagsnp = .99)

snp_idxs_1p13_large <- select_region_idxs(qc$gdata, chromosome = 1,
  position_min = 114e6, n_snps = 100)
plt <- gtable_ld_gdata(qc$gdata, snp_idxs_1p13_large)

}
