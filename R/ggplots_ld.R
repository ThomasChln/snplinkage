
#' Ggplot linkage disequilibrium
#'
#' Display SNP r2 correlations using points or diamonds with text.
#'
#' @param df_ld       Data frame with columns SNP_A, SNP_B, and R2.
#'                    As returned by the snprelate_ld function.
#' @param diamonds    Should the values be displayed as diamonds or points ?
#'                    Default is TRUE for less than 40 SNPs.
#' @param reverse     Reverse the display (horizontal symmetry)
#' @param reindex     If FALSE, SNPs are positionned following their IDs
#' @param point_size  Size for geom_point. Ignored if diamonds is TRUE.
#' @return ggplot
#' @export
ggplot_ld <- function(df_ld, diamonds = length(unique(df_ld$SNP_A)) < 40,
  point_size = 120 / sqrt(nrow(df_ld)),
  reverse = FALSE, reindex = TRUE) {

  check_fx_args(df_ld = '!d+', diamonds = '!B1')
  required_names <- c('SNP_A', 'SNP_B', 'R2')
  names_idxs <- match(required_names, names(df_ld))
  missing <- is.na(names_idxs)
  die_if(any(missing), paste('Df_ld must have columns:',
      paste0(required_names[missing], collapse = ', ')))
  df_ld <- as.data.frame(lapply(df_ld[required_names], as.numeric))

  if (reindex) {
    uniq_ids <- unique(unlist(df_ld[1:2]))
    df_ld[1:2] <- lapply(df_ld[1:2],
      function(ids) as.numeric(factor(ids, uniq_ids)))
  }
  args <- list(df_ld, reverse)
  fun <- paste0('.ggplot_ld_', if (diamonds) {
      'diamonds'
    } else {
      args[[3]] <- point_size
      'points'
    })
  ggplt <- do.call(fun, args)

  ggplt + theme(panel.background = element_blank(),
    panel.grid = element_blank(), axis.text = element_blank(),
    axis.title = element_blank(), axis.ticks = element_blank())
}

diamond_grid <- function(df_grid_snp, r2, window) {
  grid_snp = df_grid_snp
  grid_snp <- grid_snp[apply(grid_snp, 1, diff) <= window, ]
  grid_snp[, 1] <- grid_snp[, 1] - grid_snp[, 2]
  grid_snp[, 2] <- grid_snp[, 2] + .5 * grid_snp[, 1]
  grid_snp[, 1] <- grid_snp[, 1] - .5 * grid_snp[, 1]

  df_ld <- cbind.data.frame(grid_snp, r2)[c(2, 1, 3)]

  stats::setNames(df_ld, c('x', 'y', 'R2'))
}

.ggplot_ld_diamonds <- function(df_ld, reverse) {
  all_snps <- unlist(df_ld[1:2])
  seq_snp <- unique(all_snps)
  n_snp <- max(seq_snp)
  window <- max(apply(df_ld[1:2], 1, diff))
  df_ld = diamond_grid(df_ld[1:2], df_ld$R2, window)

  # Plot LD values
  inv_r2 <- pmin(1, pmax(0, 1 - df_ld$R2))
  df_ld$color <- grDevices::rgb(1, inv_r2, inv_r2)
  ggplt <- ggplot(df_ld, aes(.data$x, .data$y)) +
    diamond_annots(df_ld) +
    geom_text(label = round(df_ld$R2 * 100), size = 2.5)

  # Plot SNPs indexes
  df_snp <- data.frame(x = seq_snp, y = 0, color = 'gray')
  ggplt <- ggplt + diamond_annots(df_snp) +
    annotate('text', df_snp$x, df_snp$y,
      label = df_snp$x, size = 2.5, fontface = 'bold')

  ywindow <- c(- window / 2 - 0.5, 0.5)
  if (reverse) ywindow %<>% rev

  ggplt +
    coord_cartesian(xlim = c(0.5, n_snp + 0.5), ylim = ywindow, expand = FALSE)
}

.ggplot_ld_points <- function(df_ld, reverse, size) {

  # add snp indexes
  uniq_snps <- unique(unlist(df_ld[1:2]))
  df_snps <- data.frame(SNP_A = uniq_snps, SNP_B = uniq_snps)
  df_snps <- rbind(df_ld[if (reverse) 1:2 else 2:1], df_snps)

  # rotate points
  m_ld <- matrix(unlist(df_snps), 2, byrow = TRUE)
  theta <- -pi / 4
  rot_mat <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2)
  m_ld <- rot_mat %*% m_ld

  # range values
  ranges <- apply(m_ld, 1, range)
  df_snps <- lapply(list(x = 1, y = 2),
    function(idx) ((m_ld[idx, ] - ranges[1, idx]) / diff(ranges[, idx]) * length(uniq_snps)) + 0.5)

  # separate indexes and correlations
  seq_rows <- seq_len(nrow(df_ld))
  df_ld[1:2] <- lapply(df_snps, '[', seq_rows)
  df_snps <- lapply(df_snps, '[', -seq_rows)

  map <- aes(.data$SNP_A, .data$SNP_B, color = .data$R2)
  xframe <- c(0, length(uniq_snps) + 1)
  yframe <- if (reverse) rev(xframe) else xframe

  ggplot(df_ld, map) + geom_point(shape = 18, size = size) +
    annotate('point', x = df_snps$x, y = df_snps$y, shape = 18,
      color = grDevices::grey(.6), size = size * 0.9) +
    annotate('segment', x = xframe + c(-0.3, 0.3), xend = xframe[2] / 2,
      y = yframe[2], yend = yframe[1] + 0.6, lineend = 'round') +
    scale_colour_gradient(limits = 0:1, breaks = c(0, 0.5, 1), low = 'white',
      high = '#972D15', na.value = 'grey',
      guide = guide_colorbar(barheight = 5)) +
    coord_cartesian(xlim = xframe, ylim = xframe + c(0, 0.6), expand = FALSE) +
    labs(color = 'R2') +
    theme(legend.justification = c(1, reverse), legend.position = c(1, reverse))
}


#' Get diamond ggplot layer.
#'
#' Diamond ggplot layer for ggplot_ld
#'
#' @param data Data frame of 3 columns defining the diamonds
#' @param x Name of the column for horizontal positions
#' @param y Name of the column for vertical positions
#' @param color Name of the column for color values
#' @param size Radius of the diamonds
#' @return gglayers
#' @export
diamond_annots <- function(data, x = 'x', y = 'y', color = 'color', size = .5) {

  check_fx_args(x = '!C1', y = '!C1', color = '!C1')

  die_unless(all(c(x, y, color) %in% names(data)),
    paste('Data does not have columns',
      paste(c(x, y, color), collapse = ',')))

  data[c(x, y)] <- lapply(data[c(x, y)], as.numeric)
  bounds <- apply(data[c(x, y)], 1, diamond_bounds, x, y, size)

  color <- data[[color]]
  diamonds <- lapply(seq_along(bounds), function(bounds_idx) {
      bounds <- bounds[[bounds_idx]]
      annotate('polygon', bounds$X, bounds$Y, color = 'black',
        fill = color[bounds_idx])
    })

  diamonds
}


diamond_bounds <- function(diamond_coord, x, y, size) {
  bounds <- sapply(diamond_coord[c(x, y)],
    function(dim) sapply(c('-', '+'), Reduce, c(dim, size)))
  bounds <- data.frame('X' = c(bounds[, 1], rep(diamond_coord[[x]], 2)),
    'Y' = c(rep(diamond_coord[[y]], 2), bounds[, 2]))
  bounds[2:3, ] <- bounds[3:2, ]
  bounds
}

