#' select_region_idxs
#'
#' Select SNP indexes corresponding to a specific genomic region.
#'
#' @param gdata        Genotype Data object
#' @param chromosome   Chromosome to select
#' @param position_min Minimum base pair position to select
#' @param position_max Maximum base pair position to select
#' @param n_snps       Maximum number of SNPs to return
#' @param offset       Number of SNPs to offset
#' @return SNP indexes of Genotype Data object
#' @export
select_region_idxs <- function(gdata, chromosome, position_min = -Inf,
  position_max = Inf, n_snps = 0, offset = 0) {

  df_snp <- gdata_snps_annots(gdata)
  snp_idxs <- which(df_snp$chromosome == chromosome &
      df_snp$position >= position_min & df_snp$position <= position_max)

  if (offset != 0) snp_idxs %<>% utils::tail(-offset)
  if (n_snps != 0) snp_idxs %<>% utils::head(n_snps)

  snp_idxs
}

#' get_biomart_metadb
#'
#' To query gene names of SNPs, it is necessary to retrieve two objects using
#' biomaRt::useMart. First, the object required to map SNP rs identifiers to
#' ENSEMBL identifiers. Second, the object required to map ENSEMBL identifiers
#' to common gene names. The function returns a list of two slots named
#' snpmart and ensembl corresponding to each one, respectively.
#' Once obtained it is saved to a local file.
#'
#' @param filepath Path to save the biomaRt objects
#' @param host     BiomaRt Ensembl host, by default https://grch37.ensembl.org
#' @return List of slots snpmart and ensembl as detailed above
#' @export
get_biomart_metadb <- function(filepath = extdata_filepath('bmart_meta.rds'),
  host = "https://grch37.ensembl.org") {
  if (!file.exists(filepath)) {
    snpmart <- biomaRt::useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp",
      host = host)
    ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl",
      host = host)
 
    biomart_metadb <- list(snpmart = snpmart, ensembl = ensembl)
    save(biomart_metadb, file = filepath)
    return(biomart_metadb)
  }

  get(load(filepath))
}

biomart_gene_annots <- function(rsids, biomart) {

  df_snps <- biomaRt::getBM(c('refsnp_id', 'ensembl_gene_stable_id'), 
    'snp_filter', rsids, biomart$snpmart) %>%
    subset(ensembl_gene_stable_id != '')

  df_genes <- biomaRt::getBM(c('ensembl_gene_id', 'hgnc_symbol'), 
    'ensembl_gene_id', df_snps$ensembl_gene_stable_id, biomart$ensembl)

  df_annots <- merge(df_snps, df_genes, by.x = 'ensembl_gene_stable_id',
    by.y = 'ensembl_gene_id')

  df_annots <- subset(df_annots[-1], !duplicated(refsnp_id)) %>%
    merge(rsids, by = 1, all = TRUE)

  df_annots[match(rsids, df_annots$refsnp_id), ] %$%
    ifelse(is.na(hgnc_symbol) | hgnc_symbol == '', refsnp_id, hgnc_symbol)
}


#' gdata_add_gene_annots
#'
#' Add biomaRt gene annotations to Genotype Data object.
#'
#' @param gdata          Genotype Data object
#' @param snp_idxs       SNP indexes
#' @param rsids_colname  Column of SNP annotation data frame with rs identifiers
#' @param biomart_metadb List with slots snpmart and ensembl, corresponding to 
#'                       the biomart databases to query for SNP identifiers
#'                       and gene names, respectively. See get_biomart_metadb
#'                       function.
#' @return Genotype Data object
#' @export
gdata_add_gene_annots <- function(gdata, snp_idxs, rsids_colname = 'probe_id',
  biomart_metadb = get_biomart_metadb()) {

  if (is.null(gdata@snpAnnot)) {
    stop("Can't add gene annotations if the gdata SNPs data frame is NULL.")
  }
  df_snps <- gdata@snpAnnot@data

  rsids <- grep('^rs', df_snps[[rsids_colname]][snp_idxs], value = TRUE)

  annots <- biomart_gene_annots(rsids, biomart_metadb)
  if (!'gene' %in% names(df_snps)) df_snps$gene <- df_snps[[rsids_colname]]
  df_snps$gene[match(rsids, df_snps[[rsids_colname]])] <- annots

  gdata@snpAnnot@data <- df_snps 
  gdata
}

#' gdata_add_gene_annots_hladr_example
#'
#' Add HLA-DR gene annotations to Genotype Data object.
#' Convenience function for the vignette to avoid querying biomaRt on build.
#'
#' @param gdata        Genotype Data object
#' @param hla_dr_idxs  HLA-DR indexes in the example Genotype data object
#' @return Genotype Data object
#' @export
gdata_add_gene_annots_hladr_example <- function(gdata, hla_dr_idxs) {

  annots <- c("HLA-DRA",
    "HLA-DRA", "HLA-DRA", "rs2395182", "rs3129890", "rs13209234",
    "rs9268832", "rs6903608", "rs2395185", "HLA-DRB1", "HLA-DRB1",  
    "HLA-DRB1", "HLA-DRB1", "rs9271366", "rs9271568", "rs17533090",
    "rs3129763", "HLA-DQA1", "HLA-DQA1", "HLA-DQB1")

  df_snps <- gdata@snpAnnot@data
  df_snps$gene <- df_snps$probe_id
  df_snps$gene[hla_dr_idxs] <- annots
  gdata@snpAnnot@data <- df_snps 

  gdata
}

#' gdata_add_gene_annots_aim_example
#'
#' Add ancestry informative markers gene annotations to Genotype Data object.
#' Convenience function for the vignette to avoid querying biomaRt on build.
#'
#' @param gdata    Genotype Data object
#' @param aim_idxs AIM indexes in the example Genotype data object
#' @return Genotype Data object
#' @export
gdata_add_gene_annots_aim_example <- function(gdata, aim_idxs) {

  annots <- c("GABBR1", "rs11752362",
    "rs3132630", "rs6926530", "rs2233956", "CCHCR1",
    "CCHCR1", "CCHCR1", "HLA-B", "HLA-B",
    "HLA-B", "HLA-B", "HLA-B",
    "rs4413654", "BAG6", "CSNK2B",
    "LY6G5B", "LY6G6D", "SKIV2L", "rs9404942")

  df_snps <- gdata@snpAnnot@data
  df_snps$gene <- df_snps$probe_id
  df_snps$gene[as.numeric(aim_idxs)] <- annots
  gdata@snpAnnot@data <- df_snps 

  gdata
}
