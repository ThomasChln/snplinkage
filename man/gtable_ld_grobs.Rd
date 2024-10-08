% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gtable_ld.R
\name{gtable_ld_grobs}
\alias{gtable_ld_grobs}
\title{Build gtable by combining ggplots}
\usage{
gtable_ld_grobs(plots, labels_colname, title)
}
\arguments{
\item{plots}{List of ggplots}

\item{labels_colname}{Does the SNP position plot contain labels}

\item{title}{Title text string}
}
\value{
gtable of ggplots
}
\description{
Build gtable by combining ggplots
}
\examples{

library(snplinkage)

# example rnaseq data frame, 20 variables of 20 patients
m_rna = matrix(runif(20 ^ 2), nrow = 20)

# pair-wise correlation matrix
m_ld = cor(m_rna) ^ 2

# keep only upper triangle and reshape to data frame
m_ld[lower.tri(m_ld, diag = TRUE)] = NA
df_ld = reshape2::melt(m_ld) |> na.omit()

# rename for SNPLinkage
names(df_ld) = c('SNP_A', 'SNP_B', 'R2')

# visualize with ggplot_ld
gg_ld = ggplot_ld(df_ld)
# let's imagine the 20 variables came from 3 physically close regions
positions = c(runif(7, 10e5, 15e5), runif(6, 25e5, 30e5),
              runif(7, 45e5, 50e5)) |> sort()

# build the dataframe
df_snp_pos = data.frame(position = positions)
df_snp_pos$label = c(rep('HLA-A', 7), rep('HLA-B', 6), rep('HLA-C', 7))
gg_snp_pos = ggplot_snp_pos(df_snp_pos, labels_colname = 'label')

l_ggs = list(snp_pos = gg_snp_pos, ld = gg_ld)
gt_ld = gtable_ld_grobs(l_ggs, labels_colname = TRUE,
                        title = 'RNASeq correlations')
grid::grid.draw(gt_ld)

}
