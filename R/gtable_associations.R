
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

gtable_ld_associations_combine = function(ggplots, diamonds) {

  plots <- lapply(ggplots, ggplotGrob)
  insert_idx <- 10

  position_size <- grid::unit(if (diamonds) 1 else .5, 'null')
  plots$pval %<>% gtable::gtable_add_rows(position_size, insert_idx)

  plots$pval %<>% gtable::gtable_add_grob(
    gtable::gtable_filter(plots$pos, 'panel'), insert_idx + 1, 5)

  ld_size <- grid::unit(1.8, 'null')
  plots$pval %<>% gtable::gtable_add_rows(ld_size, insert_idx + 1)

  gtable::gtable_add_grob(plots$pval,
    gtable::gtable_filter(plots$ld, 'panel|guide-box'), insert_idx + 2, 5)
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

  if (!is.null(labels_colname)) {
    varnames <- df_snp[[labels_colname]]
    varnames[rank(-pvals) > n_labels] <- NA
  }

  df_assocs <- df_snp %$%
    cbind.data.frame(chromosome, y = pvals, vnames = varnames,
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
      if (!linked_area) '(Mbp)')
    list(labs(x = xlabel),
      theme(strip.background = element_blank(), strip.text = element_blank()))
  } else labs(x = 'Chromosomes (Mbp)')

  gg + geom_point(color = colors) + cowplot::theme_cowplot() +
    labs(y = '-log10(p-values)') + xlabel_gg +
      theme(axis.ticks.x = element_line(color = colors, lineend = 'square'))
}

