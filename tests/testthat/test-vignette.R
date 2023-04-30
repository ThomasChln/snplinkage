context('vignette')

test_that('vignette', {
    gds_path <- save_hgdp_as_gds()
    gdata <- load_gds_as_genotype_data(gds_path)
    qc <- snprelate_qc(gdata, tagsnp = .99)
    tex_table = capture.output(print_qc_as_tex_table(qc))
    expect_true(length(tex_table) > 0)

    snp_idxs_hladr <- select_region_idxs(qc$gdata,
      chromosome = 6, position_min = 32.5e6, n_snps = 20, offset = 9)
 
    # qc$gdata <- gdata_add_gene_annots(qc$gdata, snp_idxs_hladr)
    qc$gdata <- gdata_add_gene_annots_hladr_example(qc$gdata, snp_idxs_hladr)
 
    plt <- gtable_ld_gdata(qc$gdata, snp_idxs_hladr, labels_colname = 'gene')
    expect_is(plt, 'gtable')
 
    snp_idxs_mhc <- select_region_idxs(qc$gdata,
      chromosome = 6, position_min = 29e6, position_max = 33e6)
    df_assocs <- chisq_pvalues_gdata(qc$gdata, snp_idxs_mhc)
 
    df_top_aim <- subset(df_assocs, rank(-pvalues, ties.method = 'first') <= 20)
 
    #qc$gdata <- gdata_add_gene_annots(qc$gdata, rownames(df_top_aim))
    qc$gdata <- gdata_add_gene_annots_aim_example(qc$gdata, rownames(df_top_aim))
 
    plt <- gtable_ld_associations_gdata(df_top_aim, qc$gdata,
      labels_colname = 'gene')
    expect_is(plt, 'gtable')
 
    plt <- gtable_ld_associations_gdata(df_assocs, qc$gdata,
      labels_colname = 'gene')
    expect_is(plt, 'gtable')
  })
