
# fetch and keep only a subset of the hgdp file

fetch_hgdp <- function(paths = c('hgdp.zip', 'hgdp.txt')) {
  urls <- c('http://www.hagsc.org/hgdp/data/hgdp.zip', 
    'ftp://ftp.cephb.fr/hgdp_v3/hgdp-ceph-unrelated.out')
  for (i in 1:2) {
    if (!file.exists(paths[i])) utils::download.file(urls[i], paths[i])
  }
}

reduce_hgdp <- function(paths = hgdp_filepaths(),
  scan_selection = 'Europe|North Africa|Middle Est',
  snp_selection = 'SADs') {

  force(paths)
  setup_temp_dir()
  unzip_paths <- unzip_hgdp(paths)

  geno <- data.table::fread(unzip_paths$geno, data.table = FALSE)
  geno %<>% reduce_hgdp_scans(paths$scans, scan_selection)
  geno %<>% reduce_hgdp_snps(unzip_paths$snps, snp_selection)

  utils::write.table(geno, unzip_paths$geno, row.names = FALSE,
    col.names = FALSE, sep = '\t', quote = FALSE)

  dir.create('hgdp')
  file.rename(unzip_paths$geno, paths$geno) 
  file.rename(unzip_paths$snps, paths$snps) 
  utils::zip(paths$zip, c(paths$geno, paths$snps))
}

reduce_hgdp_scans <- function(geno, filepath, scan_selection) {

  scan_annots <- data.table::fread(filepath, data.table = FALSE)

  if (!is.null(scan_selection)) {
    scan_annots <- scan_annots[grep(scan_selection, scan_annots$Region), ]
  }

  scan_idxs <- match(scan_annots$hgdp_id, geno[1, ])
  scan_annots <- scan_annots[!is.na(scan_idxs), ]

  utils::write.table(scan_annots, filepath, row.names = FALSE,
    sep = '\t', quote = FALSE)

  geno[c(1, stats::na.omit(scan_idxs))]
}

reduce_hgdp_snps <- function(geno, filepath, snp_selection) {

  snp_annots <- data.table::fread(filepath, data.table = FALSE)

  if (!is.null(snp_selection)) {
    if (snp_selection == 'SADs') {
      # 5,000 regularly spaced; 2,656 selected from SADs related regions

      snp_sample <- seq(1, nrow(snp_annots), length.out = 5e3)
      snp_annots <- subset(snp_annots,
        (V1 %in% V1[snp_sample]) |
          (V2 == 1 & V3 > 113e6 & V3 < 115e6) |
          (V2 == 6 & V3 > 29e6 & V3 < 33e6) |
          (V2 == 8 & V3 > 11e6 & V3 < 12e6) |
          (V2 == 11 & V3 > 49e6 & V3 < 56e6))
    }
  }

  snp_idxs <- match(snp_annots$V1, geno$V1)
  snp_annots <- snp_annots[!is.na(snp_idxs), ]

  utils::write.table(snp_annots, filepath, row.names = FALSE,
    col.names = FALSE, sep = '\t', quote = FALSE)

  geno[c(1, stats::na.omit(snp_idxs)), ]
}

fetch_hgdp()
reduce_hgdp()


