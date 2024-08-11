
#' Gtable of linkage disequilibrium and associations using a GenotypeData object
#'
#' Compute linkage disequilibrium using snprelate_ld on the set of SNPs in the
#' associations data frame and call gtable_ld_associations.
#' Creates a gtable of a linkage disequilibrium, chromosomic positions, and
#' association scores ggplots.
#'
#' @param df_assocs      SNP annotation data frame with columns chromosome,
#'                       position, and as specified by parameters pvalue_colname
#'                       and optionally labels_colname.
#' @param gdata          GenotypeData object, as returned by
#'                       load_gds_as_genotype_data
#' @param pvalue_colname Column name of df_snp with association values
#' @param labels_colname Optional column name of df_snp with labels.
#'                       Set NULL to remove labels.
#' @param diamonds       Should the values be displayed as diamonds or points ?
#'                       Default is TRUE for up to 40 SNPs.
#' @param window         Window size for snprelate_ld.
#'                       Forced to the total number of SNPs if diamonds is FALSE
#' @param ...            Passed to gtable_ld_associations
#' @return gtable
#'
#' @examples
#' library(snplinkage)
#' gds_path <- save_hgdp_as_gds()
#' gdata <- load_gds_as_genotype_data(gds_path)
#' qc <- snprelate_qc(gdata, tagsnp = .99)
#'
#' snp_idxs_mhc <- select_region_idxs(qc$gdata,
#'   chromosome = 6, position_min = 29e6, position_max = 33e6)
#' df_assocs <- chisq_pvalues_gdata(qc$gdata, snp_idxs_mhc)
#'
#' df_top_aim <- subset(df_assocs, rank(-pvalues, ties.method = 'first') <= 20)
#'
#' #qc$gdata <- gdata_add_gene_annots(qc$gdata, rownames(df_top_aim))
#' qc$gdata <- gdata_add_gene_annots_aim_example(qc$gdata, rownames(df_top_aim))
#'
#' plt <- gtable_ld_associations_gdata(df_top_aim, qc$gdata,
#'   labels_colname = 'gene')
#'
#' @export
gtable_ld_associations_gdata <- function(df_assocs, gdata,
  pvalue_colname = 'pvalues', labels_colname = 'probe_id',
  diamonds = nrow(df_assocs) <= 40, window = 15, ...) {

  if (!diamonds) window = 0
  snp_idxs <- match(df_assocs$snpID, gdata@snpAnnot@data$snpID)
  df_ld <- snprelate_ld(gdata, window_size = window, snps_idx = snp_idxs,
    quiet = TRUE)

  df_assocs[[labels_colname]] <- gdata@snpAnnot@data[[labels_colname]][snp_idxs]

  gtable_ld_associations(df_assocs, df_ld, pvalue_colname = pvalue_colname,
    labels_colname = labels_colname, ...)
}

#' Gtable of linkage disequilibrium and associations
#'
#' Creates a gtable of a linkage disequilibrium, chromosomic positions, and
#' association scores ggplots.
#'
#' @param df_assocs      SNP annotation data frame with columns chromosome,
#'                       position, and as specified by parameters pvalue_colname
#'                       and optionally labels_colname.
#' @param df_ld          Data frame with columns SNP_A, SNP_B, and R2, as
#'                       returned by the snprelate_ld function.
#' @param pvalue_colname Column name of df_snp with association values
#' @param labels_colname Optional column name of df_snp with labels.
#'                       Set NULL to remove labels.
#' @param n_labels       Number of labels of most associated SNPs to display.
#' @param diamonds       Should the values be displayed as diamonds or points ?
#'                       Default is TRUE for up to 40 SNPs.
#' @param linked_area    Add a linked area to associations points.
#'                       Default same as diamonds.
#' @param point_size     Point size for ggplot_ld, ignored if diamonds is TRUE.
#' @param colors         Colors of SNPs
#' @param ...            Passed to ggplot_associations
#' @return gtable
#' @export
gtable_ld_associations <- function(df_assocs, df_ld, pvalue_colname = 'pvalues',
  labels_colname = 'probe_id', n_labels = 5,
  diamonds = nrow(df_assocs) <= 40, linked_area = diamonds,
  point_size = 150 / nrow(df_assocs),
  colors = snp_position_colors(nrow(df_assocs)), ...) {

  gg_pval <- ggplot_associations(df_snp = df_assocs, n_labels = n_labels,
    linked_area = linked_area, labels_colname = labels_colname,
    pvalue_colname = pvalue_colname,
    colors = if (linked_area) colors else 'black', ...)

  if (diamonds) {
    gg_pval <- gg_pval + theme(axis.text.x = element_blank())
  }

  title <- gg_pval$labels$x %>% gsub(' (Mbp)', '', ., fixed = TRUE) %>%
    paste('-', nrow(df_assocs), 'SNPs')
  gg_pval <- gg_pval + labs(title = title, x = NULL)

  gg_ld <- ggplot_ld(df_ld, point_size = point_size)
  gg_pos <- ggplot_snp_pos(df_assocs, diamonds, if (diamonds) labels_colname,
    colors = colors)

  list(pval = gg_pval, pos = gg_pos, ld = gg_ld) %>%
    gtable_ld_associations_combine(diamonds)
}

#' Build gtable by combining ggplots 
#'
#' @param ggplots  List of ggplots
#' @param diamonds Does the LD visualization use diamond-type layout
#' @return gtable of ggplots
#'
#' @examples
#'
#' library(snplinkage)
#'
#' # example rnaseq data frame, 20 variables of 20 patients
#' m_rna = matrix(runif(20 ^ 2), nrow = 20)
#'
#' # pair-wise correlation matrix
#' m_ld = cor(m_rna) ^ 2
#'
#' # keep only upper triangle and reshape to data frame
#' m_ld[lower.tri(m_ld, diag = TRUE)] = NA
#' df_ld = reshape2::melt(m_ld) |> na.omit()
#'
#' # rename for SNPLinkage
#' names(df_ld) = c('SNP_A', 'SNP_B', 'R2')
#'
#' # visualize with ggplot_ld
#' gg_ld = ggplot_ld(df_ld)
#'
#' # let's imagine the 20 variables came from 3 physically close regions
#' positions = c(runif(7, 10e5, 15e5), runif(6, 25e5, 30e5),
#'               runif(7, 45e5, 50e5)) |> sort()
#'
#' # build the dataframe
#' df_snp_pos = data.frame(position = positions)
#' df_snp_pos$label = c(rep('HLA-A', 7), rep('HLA-B', 6), rep('HLA-C', 7))
#'
#' gg_pos_biplot = ggplot_snp_pos(df_snp_pos, labels_colname = 'label',
#'                                upper_subset = TRUE)
#'
#' # let's assume HLA-B is more associated with the outcome than the other genes
#' pvalues = c(runif(7, 1e-3, 1e-2), runif(6, 1e-8, 1e-6), runif(7, 1e-3, 1e-2))
#' log10_pvals = -log10(pvalues)
#'
#' # we can reuse the df_snp_pos object
#' df_snp_pos$pvalues = log10_pvals
#' 
#' # add the chromosome column
#' df_snp_pos$chromosome = 6
#'
#' gg_assocs = ggplot_associations(df_snp_pos, labels_colname = 'label',
#'                                 linked_area = TRUE, nudge = c(0, 0.5),
#'                                 n_labels = 12)
#'
#' l_ggs = list(pos = gg_pos_biplot, ld = gg_ld, pval = gg_assocs)
#' gt_ld = gtable_ld_associations_combine(l_ggs, diamonds = TRUE)
#' grid::grid.draw(gt_ld)
#'
#' @export
gtable_ld_associations_combine = function(ggplots, diamonds) {

  plots <- lapply(ggplots, ggplotGrob)
  plot_pos = if (getRversion() >= "4.3.3") {
    list(height = 12, width = 7)
  } else {
    list(height = 10, width = 5)
  }

  position_size <- grid::unit(if (diamonds) 1 else .5, 'null')
  plots$pval %<>% gtable::gtable_add_rows(position_size, plot_pos$height)

  plots$pval %<>% gtable::gtable_add_grob(
    gtable::gtable_filter(plots$pos, 'panel'), plot_pos$height + 1, plot_pos$width)

  ld_size <- grid::unit(1.8, 'null')
  plots$pval %<>% gtable::gtable_add_rows(ld_size, plot_pos$height + 1)

  gtable::gtable_add_grob(plots$pval,
    gtable::gtable_filter(plots$ld, 'panel|guide-box'), plot_pos$height + 2, plot_pos$width)
}

#' Ggplot associations
#'
#' Get SNPs associations ggplot, either as points or as a linked area.
#' Optionally add labels to most associated points using ggrepel.
#'
#' @param df_snp         SNP annotation data frame with columns chromosome,
#'                       position, and as specified by parameters pvalue_colname
#'                       and optionally labels_colname.
#' @param pvalue_colname Column name of df_snp with association values
#' @param labels_colname Optional column name of df_snp with labels.
#'                       Set to NULL to remove.
#' @param n_labels       Number of labels of most associated points to display.
#' @param nudge          Nudge parameter passed to ggrepel::geom_label_repel.
#' @param linked_area    Add a linked area to associations points, default FALSE
#' @param byindex        Display by SNP index or chromosomic position (default)
#' @param colors         Colors of SNPs
#' @return ggplot
#' @export
ggplot_associations <- function(df_snp, pvalue_colname = 'pvalues',
  labels_colname = 'probe_id', n_labels = 10, nudge = c(0, 1),
  linked_area = FALSE, byindex = linked_area,
  colors = if (linked_area) snp_position_colors(nrow(df_snp)) else 'black') {
  
  pvals <- df_snp[[pvalue_colname]]

  assoc_labels = if (!is.null(labels_colname)) {
    varnames <- df_snp[[labels_colname]]
    varnames[rank(-pvals) > n_labels] <- NA
    varnames
  }

  df_assocs <- df_snp %$%
    cbind.data.frame(chromosome, y = pvals, vnames = assoc_labels,
      x = if (byindex) seq_len(nrow(df_snp)) else position / 1e6)

  yrange <- range(df_assocs$y)
  yrange[1] <- yrange[1] - yrange[1] / 10
  yrange[2] <- yrange[2] + yrange[2] / 10

  gg <- ggplot(df_assocs, aes(x, y)) +
    facet_wrap(~ chromosome, nrow = 1, scales = 'free_x') +
    coord_cartesian(if (linked_area) c(.5, nrow(df_snp) + .5),
      yrange, expand = !linked_area)

  if (linked_area) {
    gg <- gg +
      geom_area(fill = grDevices::grey(.5), alpha = .5, color = 'black') +
      geom_segment(aes(xend = x), yend = 0, color = colors) +
      scale_x_continuous(breaks = df_assocs$x)
  }

  if (!is.null(labels_colname)) {
    gg <- gg + ggrepel::geom_label_repel(label = df_assocs$vnames, na.rm = TRUE,
        label.size = .1, label.padding = unit(.1, 'lines'), color = colors,
        nudge_x = nudge[1], nudge_y = nudge[2])
  }

  xlabel_gg <- if (length(unique(df_assocs$chromosome)) == 1) {
    xlabel <- paste('Chromosome', unique(df_assocs$chromosome),
      if (!byindex) '(Mbp)')
    list(labs(x = xlabel),
      theme(strip.background = element_blank(), strip.text = element_blank()))
  } else labs(x = 'Chromosomes (Mbp)')

  gg + geom_point(color = colors) + cowplot::theme_cowplot() +
    labs(y = '-log10(p-values)') + xlabel_gg +
      theme(axis.ticks.x = element_line(color = colors, lineend = 'square'))
}

