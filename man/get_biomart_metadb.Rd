% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/biomart_annots.R
\name{get_biomart_metadb}
\alias{get_biomart_metadb}
\title{get_biomart_metadb}
\usage{
get_biomart_metadb(
  filepath = extdata_filepath("bmart_meta.rds"),
  host = "https://grch37.ensembl.org"
)
}
\arguments{
\item{filepath}{Path to save the biomaRt objects}

\item{host}{BiomaRt Ensembl host, by default https://grch37.ensembl.org}
}
\value{
List of slots snpmart and ensembl as detailed above
}
\description{
To query gene names of SNPs, it is necessary to retrieve two objects using
biomaRt::useMart. First, the object required to map SNP rs identifiers to
ENSEMBL identifiers. Second, the object required to map ENSEMBL identifiers
to common gene names. The function returns a list of two slots named
snpmart and ensembl corresponding to each one, respectively.
Once obtained it is saved to a local file.
}
