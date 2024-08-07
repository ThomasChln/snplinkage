#' Gtable of linkage disequilibrium and positions using a GenotypeData object
#'
#' Compute linkage disequilibrium using snprelate_ld on a set of SNP indexes and
#' call gtable_ld.
#' Two parameters are available to compute and compare minor allele frequency
#' filtering and TagSNP selection by displaying two LD plots with their
#' positions in the center.
#' The maf and r2 parameters are used similarly and as follows:
#'    - compare baseline with MAF 5%
#'    gtable_ld(gdata, snps_idx, maf = 0.05)
#'    - compare baseline with TagSNP r2 = 0.8
#'    gtable_ld(gdata, snps_idx, r2 = 0.8)
#'    - compare 5% MAF with 5% MAF and r2 = 0.8
#'    gtable_ld(gdata, snps_idx, maf = c(0.05, 0.05), r2 = 0.8)
#'    - compare MAF 5%, r2 0.8 with MAF 10%, r2 = 0.6
#'    gtable_ld(gdata, snps_idx, maf = c(0.05, 0.1), r2 = c(0.8, 0.6))
#'
#' @param gdata        GenotypeData object returned by load_gds_as_genotype_data
#' @param snps_idx     SNPs indexes to select
#' @param maf          Minor allele frequency threshold(s), see description
#' @param r2           TagSNP r2 threshold(s), see description
#' @param diamonds     Display the values as diamonds or as points
#'                     Default is TRUE for less than 40 SNPs.
#' @param window       Window size for snprelate_ld.
#'                     Forced to the total number of SNPs if diamonds is FALSE
#' @param autotitle    Set title to feature selection method(s), number of SNPs
#'                     and chromosome
#' @param autotitle_bp Set biplot title to feature selection method(s), number
#'                     of SNPs and chromosome
#' @param double_title Logical, if false (default) keep only biplot title
#' @param ...         Passed to gtable_ld
#' @return gtable of ggplots
#'
#' @examples
#' library(snplinkage)
#' gds_path <- save_hgdp_as_gds()
#' gdata <- load_gds_as_genotype_data(gds_path)
#' qc <- snprelate_qc(gdata, tagsnp = .99)
#'
#' snp_idxs_1p13_large <- select_region_idxs(qc$gdata, chromosome = 1,
#'   position_min = 114e6, n_snps = 100)
#' plt <- gtable_ld_gdata(qc$gdata, snp_idxs_1p13_large)
#'
#' @export
gtable_ld_gdata <- function(gdata, snps_idx, maf = NULL, r2 = NULL,
  diamonds = length(snps_idx) < 40, window = 15, autotitle = TRUE,
  autotitle_bp = TRUE, double_title = FALSE, ...) {

  # used to verify that maf and r2 are null or numerics of length 1 or 2
  # check_fx_args(maf = 'N1:2', r2 = 'N1:2')
  check_fx_args(snps_idx = '!I+', window = '!N1')
  die_unless(inherits(gdata, 'GenotypeData'), 'gdata must be a GenotypeData')
  chromosome <- unique(getChromosome(gdata, snps_idx))
  die_unless(length(chromosome) == 1,
    'Only one chromosome must be provided by the snps_idx parameter')
  if (!diamonds) window = 0

  # Manage MAF and TagSNP options
  # maf for both plots
  titles <- rep(paste('Chromosome', chromosome, '-'), 2) %>% as.list %>%
    stats::setNames(c('ld', 'biplot'))
  titles$ld <- paste(titles$ld, length(snps_idx), 'SNPs')
  df_snp <- gdata_snps_annots(gdata)

  if (length(maf) == 2) {
    titles$ld <- paste0(titles$ld, ' MAF ', maf[2] * 100, '% -')
    mafs <- snprelate_allele_frequencies(gdata, snps_idx, quiet = TRUE)
    snps_idx <- snps_idx[mafs$maf > maf[2]]
  }
  # tagsnp for both
  if (length(r2) == 2) {
    titles$ld <- paste(titles$ld, 'TagSNP', r2[2], '-')
    tsnp_ids <- snprelate_ld_select(gdata, min_r2 = r2[2],
      snps_idx = snps_idx, quiet = TRUE)
    snps_idx <- match(tsnp_ids[[1]], df_snp$snpID)
  }
  # maf on second plot
  maf_subset <- if (length(maf) > 0) {
      titles$biplot <- paste0(titles$biplot, ' MAF ', maf[1] * 100, '% -')
      mafs <- snprelate_allele_frequencies(gdata, snps_idx, quiet = TRUE)
      mafs$maf > maf[1]
    } else {
      TRUE
    }
  # tagsnp on second plot
  tsnp_subset <- if (length(r2) > 0) {
      titles$biplot <- paste(titles$biplot, 'TagSNP', r2[1], '-')
      tsnp_ids <- snprelate_ld_select(gdata, min_r2 = r2[1],
        snps_idx = snps_idx[maf_subset], quiet = TRUE)
      tsnp_idxs <- match(tsnp_ids[[1]], df_snp$snpID)
      maf_subset <- TRUE
      snps_idx %in% tsnp_idxs
    } else {
      TRUE
    }
  biplot_subset <- maf_subset & tsnp_subset
  biplot_subset <- if (!all(biplot_subset)) {
    titles$biplot <- paste(titles$biplot, sum(biplot_subset), 'SNPs')
    if (!double_title) autotitle <- FALSE
    which(biplot_subset)
  }

  df_ld <- snprelate_ld(gdata, window, snps_idx = snps_idx, quiet = TRUE)

  gtable_ld(df_ld, df_snp, biplot_subset,
    title = if (autotitle) titles$ld else '',
    title_biplot =  if (autotitle_bp) titles$biplot else '', ...)
}

#' Gtable of linkage disequilibrium and chromosomic positions
#'
#' Creates a gtable of linkage disequilibrium and chromosomic positions
#' ggplots. A biplot_subset parameter is available to add a second linkage
#' disequibrium ggplot to visualize the effect of a SNP selection.
#'
#' @param df_ld          Data frame returned by snprelate_ld 
#' @param df_snp         SNP annotations with columns snpID and position
#' @param biplot_subset  SNP indexes of the subset for the second ld plot
#' @param labels_colname Column name of df_snp to use as SNP labels 
#' @param diamonds       Display the values as diamonds or as points
#'                       Default is TRUE for less than 40 SNPs.
#' @param point_size     Size for geom_point. Ignored if diamonds is TRUE.
#' @param title          Plot title
#' @param title_biplot   Optional biplot title
#' @param ...            Passed to ggplot_ld
#' @return gtable of ggplots
#'
#' @examples
#' library(snplinkage)
#' gds_path <- save_hgdp_as_gds()
#' gdata <- load_gds_as_genotype_data(gds_path)
#' qc <- snprelate_qc(gdata, tagsnp = .99)
#' snp_idxs_8p23 <- select_region_idxs(qc$gdata, chromosome = 8,
#'   position_min = 11e6, position_max = 12e6)
#'
#' df_ld <- snprelate_ld(qc$gdata, snps_idx = snp_idxs_8p23, quiet = TRUE)
#' plt <- gtable_ld(df_ld, df_snp = gdata_snps_annots(qc$gdata))
#'
#' @export
gtable_ld <- function(df_ld, df_snp, biplot_subset = NULL,
  labels_colname = NULL, diamonds = length(unique(df_ld$SNP_A)) < 40,
  point_size = ifelse(is.null(biplot_subset), 120, 80) / sqrt(nrow(df_ld)),
  title = '', title_biplot = '', ...) {

  unique_ids <- unique(unlist(df_ld[1:2]))
  df_snp <- df_snp[match(unique_ids, df_snp$snpID), ]

  df_ld[1:2] <- lapply(df_ld[1:2],
    function(ids) as.numeric(factor(ids, unique_ids)))
  
  plots <- list(snp_pos = ggplot_snp_pos(df_snp, biplot_subset, labels_colname),
    ld = ggplot_ld(df_ld, ..., point_size = point_size) + labs(title = title))

  # Biplot, for diamonds remove subset, for points set to NA
  if (!is.null(biplot_subset)) {
    snp_subset <- with(df_ld,
      SNP_A %in% biplot_subset & SNP_B %in% biplot_subset)
    if (diamonds) df_ld <- df_ld[snp_subset, ] else df_ld$R2[!snp_subset] = NA

    plots$biplot <- ggplot_ld(df_ld, diamonds, ..., reverse = TRUE,
      reindex = FALSE, point_size = point_size) +
      labs(title = title_biplot)
  }

  gtable_ld_grobs(plots, labels_colname, title)
}

# build gtable by combining grobs
gtable_ld_grobs <- function(plots, labels_colname, title) {

  ## gtable probably had an update
  # in r 4.3.0 plot_pos is 7, 5 out of 12
  # after, 9, 7, out of 16
  plot_pos = if (getRversion() >= "4.3.3") {
    list(height = 9, width = 7, bottom = 17)
  } else {
    list(height = 7, width = 5, bottom = 13)
  }

  plots <- lapply(plots, ggplotGrob)
  is_biplot <- 'biplot' %in% names(plots)
  ld_relative_size <- if (!is.null(labels_colname) || is_biplot) 3 else 7

  plots$snp_pos <- gtable::gtable_add_rows(plots$snp_pos,
    grid::unit(ld_relative_size, 'null'), plot_pos$height)
  plots$snp_pos <- gtable::gtable_add_grob(plots$snp_pos,
    gtable::gtable_filter(plots$ld, 'panel|guide-box'),
    plot_pos$height + 1, plot_pos$width)

  # set ld title to bottom if biplot
  if (title != '') {
    title_pos <- if (is_biplot) plot_pos$bottom else 2
    plots$snp_pos <- gtable::gtable_add_rows(plots$snp_pos,
      grid::unit(0.2, 'null'), title_pos)
    plots$snp_pos <- gtable::gtable_add_grob(plots$snp_pos,
      gtable::gtable_filter(plots$ld, 'title'), title_pos + 1, plot_pos$width)
  }

  if (is_biplot) {
    plots$snp_pos <- gtable::gtable_add_rows(plots$snp_pos,
      grid::unit(3, 'null'), 2)
    plots$snp_pos <- gtable::gtable_add_grob(plots$snp_pos,
      gtable::gtable_filter(plots$biplot, 'panel|guide-box|title'), 3, plot_pos$width)
  }

  plots$snp_pos
}

#' Ggplot SNPs position
#'
#' Get SNPs position ggplot with mappings to combine with other ggplots.
#' Optionally add labels and an upper subset.
#'
#' @param df_snp         SNP annotation data frame with a column named position
#'                       and, if specified, one named as the labels_colname
#'                       parameter.
#' @param upper_subset   Subset of df_snp for the positions on the upper side
#' @param labels_colname Optional column name of df_snp to use as SNP labels.
#' @param colors         Colors for each SNP
#' @return ggplot
#' @export
ggplot_snp_pos <- function(df_snp, upper_subset = NULL, labels_colname = NULL,
  colors = snp_position_colors(nrow(df_snp))) {

  check_fx_args(df_snp = '!d+')
  required_names <- c('position', labels_colname)
  die_unless(required_names %in% names(df_snp),
    paste('df_snp must have columns',
      paste0(required_names, collapse = ', ')))

  n_snp <- nrow(df_snp)
  range_pos <- range(df_snp$position)

  # Scale positions
  pos_size <- 0.15 * n_snp
  xrange <- c(1 + pos_size, n_snp - pos_size)
  df_snp$position <- (df_snp$position - range_pos[1]) * diff(xrange)
  df_snp$position <- df_snp$position / diff(range_pos) + xrange[1]
  df_snp$index <- seq(xrange[1], xrange[2], length = n_snp)
  df_snp$seq <- seq_len(n_snp)

  # main plot
  map <- aes(.data$position, xend = .data$index)
  yrange <- c(-2, 1.5)
  xrange[1:2] <- c(xrange[1] - 0.25, xrange[2] + 0.25)
  plot <- ggplot(df_snp, map) +
    geom_segment(aes(xend = .data$position), y = 0, yend = 1,
      color = colors) +
    geom_segment(aes(xend = .data$seq), y = -0.025, yend = yrange[1],
      color = colors, lineend = 'round')

  # position annotations
  annots_params <- list('text', y = .5)
  xpos <- c(xrange[1] / 2, xrange[2] + (n_snp - xrange[2]) / 2)
  range_pos <- paste0(.format_eng_range(range_pos), 'bp')
  annots <- lapply(1:2, function(idx) {
      annots_params[c('x', 'label')] <- list(xpos[idx], range_pos[idx])
      do.call(annotate, annots_params)
    })
  plot <- plot + annots

  # names and y scale
  if (!is.null(labels_colname)) {
    yrange[1] <- -7
    plot <- plot + list(geom_text(x = df_snp$seq - 0.25, y = -5, color = colors,
        label = df_snp[[labels_colname]], size = 2.5, angle = 90),
      geom_segment(x = df_snp$seq, xend = df_snp$seq,
        y = yrange[1], yend = -2, color = colors, lineend = 'round'))
  }

  # optional second position map
  if (!is.null(upper_subset)) {
    colors <- colors[upper_subset]
    df_snp_subset <- df_snp[upper_subset, ]
    yrange[2] <- 3
    plot <- plot + geom_segment(aes(xend = .data$seq), df_snp_subset,
      y = 1.025, yend = yrange[2] - 0.05, color = colors, lineend = 'round')
  }

  # borders and theme
  plot +
    annotate('segment', x = rep(xrange[1], 2), xend = rep(xrange[2], 2),
      y = 0:1, yend = 0:1, lineend = 'square') +
    annotate('segment', x = xrange[1:2], xend = xrange[1:2],
      y = rep(0, 2), yend = rep(1, 2)) +
    coord_cartesian(xlim = c(0.5, n_snp + 0.5),
      ylim = c(yrange[1], yrange[2]),
      expand = FALSE) +
    theme(panel.background = element_blank(), text = element_blank(),
      axis.ticks = element_blank(), panel.grid = element_blank())
}

snp_position_colors <- function(n_snps,
  colors = c('#D8B70A', '#02401B', '#A2A475', '#972D15', '#81A88D')) {
  colors[seq_len(n_snps) %% length(colors) + 1]
}

.format_eng_range <- function(nums) {
  bases <- floor(log(nums, 1e3))
  base_strs <- bases
  base_strs[base_strs > 2] <- 2
  base_strs <- list(NULL, 'k', 'M')[base_strs + 1]
  nums <- nums / 1e3 ^ bases

  # keep 3 digits
  digit_base <- 10 ^ (3 - nchar(floor(nums)))
  nums <- nums * digit_base
  nums <- c(floor(nums[1]), ceiling(nums[2])) / digit_base
  nums <- paste(nums, base_strs)
}
