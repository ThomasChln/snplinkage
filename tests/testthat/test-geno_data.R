context('genotype data')
snplinkage:::setup_temp_dir()

test_actg_tsv_to_gdata  <- function() {


  paths <- system.file('extdata', paste0('hgdp.', c('zip', 'txt')),
    package = 'snplinkage')
  zippaths <- file.path('hgdp',
    c('HGDP_FinalReport_Forward.txt', 'HGDP_Map.txt'))
  txts_paths <- unzip(paths[1], zippaths, junkpaths = TRUE)
  actg_gdata <- snplinkage:::actg_tsv_to_gdata(txts_paths[1], paths[2],
    c('scan_id', 'gender', 'population', 'geographic_origin', 'region'),
    txts_paths[2])

  expect_is(actg_gdata, 'GenotypeData')

  # Open all
  l_files <- lapply(c(txts_paths[1], paths[2], txts_paths[2]),
    data.table::fread, sep = '\t', data.table = FALSE)
  file.remove(txts_paths)
  geno_test <- l_files[[1]][-1]
  geno_test <- setNames(geno_test[-1, ], geno_test[1, ])
  scans_test <- l_files[[2]]
  snps_test <- l_files[[3]]

  # Check SNPs order
  snps <- actg_gdata@snpAnnot@data$probe_id
  # Order following chromosome and position
  snps_test[2] <- as.integer(factor(snps_test[[2]],
      c(1:22, 'X', 'XY', 'Y', 'M'), nmax = 26, exclude = NULL))
  snps_order <- order(snps_test[[2]], snps_test[[3]])
  snps_test <- snps_test[[1]][snps_order]
  expect_true(all(snps == snps_test))

  # Check scans order
  scans <- actg_gdata@scanAnnot@data$scan_id
  # Discard missings
  scans_test_merge <- match(scans_test[[1]], names(geno_test))
  scans_test <- scans_test[[1]][!is.na(scans_test_merge)]
  expect_true(all(scans == scans_test))

  geno <- snplinkage:::fetch_genotypes(actg_gdata)
  geno[is.na(geno)] <- 3L
  geno_test <- geno_test[snps_order, ]
  geno_test <- geno_test[, na.omit(scans_test_merge)]
  colnames(geno_test) = scans_test
  geno_test <- t(apply(geno_test, 1, snplinkage:::actg_to_numeric, '--'))
  geno_test[is.na(geno_test)] <- 3L
  expect_true(all(geno == geno_test))
}
test_that('Convert hgdp text files to gds', test_actg_tsv_to_gdata())

