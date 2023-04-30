context('SNPRelate LD')

gds <- snplinkage:::save_hgdp_as_gds()
hgdp_gdata <- snplinkage:::load_gds_as_genotype_data(gds)
GDATA <- hgdp_gdata
on.exit(GWASTools::close(GDATA))

.snprelate_allele_frequencies <- function() {
  freqs <- snprelate_allele_frequencies(GDATA)

  expect_equal(dim(freqs), c(nsnp(GDATA), 6))

  with(freqs, expect_true(all(allele1 == maf | allele2 == maf)))

  sums <- with(freqs, allele1 + allele2)
  delta <- max(abs(sums - 1))
  expect_equal(delta, 0)

  freqs2 <- snprelate_allele_frequencies(GDATA, 7:102)
  expect_equivalent(freqs2, freqs[7:102, ])

  # scans_idx
  freqs3 <- snprelate_allele_frequencies(GDATA, 1:10, 1)
  expect_true(all(freqs3$missing == 0))
  expect_true(all(freqs3$maf == 0 | freqs3$maf == 0.5))
}
test_that("snprelate_allele_frequencies", .snprelate_allele_frequencies())



..convert_snpgdsLDMat_output <- function() {
  convert <- snplinkage:::.convert_snpgdsLDMat_output
  # dummy ld mat
  n <- 1000
  k <- 5
  m <- matrix(seq_len(n * k), nrow = k, byrow = FALSE)
  all_ind <- which(!is.na(m), arr.ind = TRUE)
  bads <- which(all_ind[, 1] + all_ind[, 2] > n, arr.ind = TRUE)
  m[all_ind[bads]] <- NaN

  ids <- 1:n
  system.time(df <- convert(m, 0, ids))
  .get <- function(m, df, row) {
    i <- as.integer(df[row, 'SNP_A'])
    j <- as.integer(df[row, 'SNP_B'])
    m[j - i, i]
  }

  row <- 117
  expect_identical(df[row, 'R2'], .get(m, df, row))
  row <-   499
  expect_identical(df[row, 'R2'], .get(m, df, row))

  ids <- paste0('SNP_' , 1:n)
  system.time(df <- convert(m, 4991, ids))
  truth <- data.frame(SNP_A = "SNP_999", SNP_B = "SNP_1000", R2 = 4991,
    stringsAsFactors = FALSE)
  expect_identical(df, truth)

}
test_that(".convert_snpgdsLDMat_output", ..convert_snpgdsLDMat_output())

.snprelate_ld <- function() {
  ld <- snprelate_ld(GDATA, snps_idx = 1:100, window_size = 2, threads = 1,
    quiet = TRUE)

  expect_equal(dim(ld), c(99, 3))
  expect_equal(names(ld), c('SNP_A', 'SNP_B', 'R2'))
}
test_that("snprelate_ld", .snprelate_ld())

.snprelate_ld_select <- function() {
  setup_temp_dir(TRUE)

  gdata <- hgdp_gdata

  # extreme setting => select one per chromosome !
  ws <- 20000
  set.seed(1)
  res <- snprelate_ld_select(gdata, window_size = ws, window_length = 1e9,
    min_r2 = 0, quiet = TRUE, remove.monosnp = TRUE)
  # N.B: it may happen, for some seeds, that several snps per chromosome are
  # selected if their r^2 is exactly 0
  expect_equal(length(res), 24)

  # window length
  set.seed(1)
  res_wl2 <- snprelate_ld_select(gdata, window_size = 2, window_length = 1e9,
    min_r2 = 0.1, quiet = TRUE)
  set.seed(1)
  res_wl3 <- snprelate_ld_select(gdata, window_size = 3, window_length = 1e9,
    min_r2 = 0, quiet = TRUE)
  # increasing window length offer more opportunity to discover redundant snps
  # so less tag snps
  expect_gt(length(unlist(res_wl2)), length(unlist(res_wl3)))

  wl <- 1e6L
  MIN_R2 <- 0.1
  chr_by_id <- getChromosome(gdata)
  ws <- max(table(chr_by_id))

  set.seed(1)
  tagged_by_chr <- snprelate_ld_select(gdata,
    window_size = ws, window_length = wl, min_r2 = MIN_R2,
    quiet = TRUE)
  chrs <- names(tagged_by_chr)
  for (chr in chrs) {
    tagged <- as.character(tagged_by_chr[[chr]])
    chr_number <- as.integer(sub('chr', '', chr))
    all_in_chr <- which(chr_by_id == chr_number)
    pruned <- setdiff(as.character(getSnpID(gdata, all_in_chr)), tagged)
    ld <- snprelate_ld(gdata, snps_idx = all_in_chr,
      window_size = length(all_in_chr), min_r2 = MIN_R2, quiet = TRUE)

    expect_equal(names(ld), c('SNP_A', 'SNP_B', 'R2'))
  }

}
test_that('snprelate_ld_select', .snprelate_ld_select())

