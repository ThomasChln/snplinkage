context('ggplots LD')


pdf(NULL)

gds <- save_hgdp_as_gds()
hgdp_gdata <- load_gds_as_genotype_data(gds)

on.exit(GWASTools::close(hgdp_gdata))
df_snp = hgdp_gdata@snpAnnot@data

df_ld_diams <- snprelate_ld(hgdp_gdata, window = 3, snps_idx = 1:5, quiet = TRUE)
df_ld_pts <- snprelate_ld(hgdp_gdata, 0, snps_idx = 1:50, quiet = TRUE)

test_ggplot_ld <- function() {

  # diamonds
  ggplt <- ggplot_ld(df_ld_diams)
  expect_is(ggplot2::ggplotGrob(ggplt), 'gtable')

  # points
  ggplt <- ggplot_ld(df_ld_pts)
  expect_is(ggplot2::ggplotGrob(ggplt), 'gtable')

  # reverse reindex, for biplots
  ggplt <- ggplot_ld(df_ld_pts, reverse = TRUE, reindex = FALSE)
  expect_is(ggplot2::ggplotGrob(ggplt), 'gtable')
}
test_that('ggplot_ld', test_ggplot_ld())

test_ggplot_snp_pos <- function() {

  # minimal call
  ggplt <- ggplot_snp_pos(df_snp[1:10, ])
  expect_is(ggplot2::ggplotGrob(ggplt), 'gtable')

  # add upper subset
  ggplt <- ggplot_snp_pos(df_snp[1:10, ], 1:5)
  expect_is(ggplot2::ggplotGrob(ggplt), 'gtable')

  # add SNP labels
  ggplt <- ggplot_snp_pos(df_snp, labels_colname = 'snpID')
  expect_is(ggplot2::ggplotGrob(ggplt), 'gtable')
}
test_that('ggplot_snp_pos', test_ggplot_snp_pos())

test_gtable_ld <- function() {

  # with biplot
  plt <- gtable_ld(df_ld_pts, df_snp, 1:40)
  expect_is(plt, 'gtable')

  # with labels
  plt <- gtable_ld(df_ld_diams, df_snp, labels_colname = 'snpID')
  expect_is(plt, 'gtable')
}
test_that('gtable_ld', test_gtable_ld())

test_gtable_ld_gdata <- function() {

  plt <- gtable_ld_gdata(hgdp_gdata, 1:8)
  expect_is(plt, 'gtable')

  # compare MAF 10% r2 .05 with MAF 20% r2 .025
  plt <- gtable_ld_gdata(hgdp_gdata, 1:8, c(.2, .1), c(.025, .05))
  expect_is(plt, 'gtable')
}
test_that('gtable_ld_gdata', test_gtable_ld_gdata())

snp_idxs_mhc <- select_region_idxs(hgdp_gdata,
  chromosome = 6, position_min = 31e6, position_max = 32e6)
df_snp_associations <- chisq_pvalues_gdata(hgdp_gdata, snp_idxs_mhc)
df_top_aim <- subset(df_snp_associations, rank(-pvalues) <= 30)

test_ggplot_associations <- function() {

  # large, points
  ggplt <- ggplot_associations(df_snp_associations)
  expect_is(ggplot2::ggplotGrob(ggplt), 'gtable')

  # small, linked area 
  ggplt <- ggplot_associations(df_top_aim, byindex = TRUE)
  expect_is(ggplot2::ggplotGrob(ggplt), 'gtable')
}
test_that('ggplot_associations', test_ggplot_associations())

  
test_gtable_ld_associations_gdata <- function() {
  # small, linked area and diamonds
  plt <- gtable_ld_associations_gdata(df_top_aim, hgdp_gdata)
  expect_is(plt, 'gtable')

  # large, points
  plt <- gtable_ld_associations_gdata(df_snp_associations, hgdp_gdata)
  expect_is(plt, 'gtable')
}
test_that('gtable_ld_associations_gdata', test_gtable_ld_associations_gdata())

