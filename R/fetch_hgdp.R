
#' save_hgdp_as_gds
#'
#' Save the HGDP SNP data text file as a Genomic Data Structure file
#'
#' @param paths    Paths of the zip, txt, and gds files
#' @param outpath  Output GDS file path
#' @param ...      Passed to save_genotype_data_as_gds
#' @return Path of the saved gds file
#' @export
save_hgdp_as_gds <- function(paths = hgdp_filepaths(), outpath = tempfile(),
  ...) {

  unzip_paths <- unzip_hgdp(paths)
  actg_gdata <- actg_tsv_to_gdata(unzip_paths$geno, paths$scans,
    c('scan_id', 'gender', 'population', 'geographic_origin', 'region'),
    unzip_paths$snps)
  lapply(unzip_paths, file.remove)
  save_genotype_data_as_gds(actg_gdata, outpath, quiet = TRUE, ...)

  outpath
}

extdata_filepath = function(filepath) {
  file.path(system.file('extdata', package = 'snplinkage'), filepath)
}

unzip_hgdp = function(paths) {
  utils::unzip(paths$zip, c(paths$geno, paths$snps), junkpaths = TRUE) %>%
    as.list %>% stats::setNames(c('geno', 'snps'))
}

hgdp_filepaths = function() {
  list(zip = extdata_filepath('hgdp.zip'), scans = extdata_filepath('hgdp.txt'),
    geno = 'hgdp/HGDP_FinalReport_Forward.txt', snps = 'hgdp/HGDP_Map.txt')
}

