
#' Load GDS as Genotype Data
#'
#' Open a connection to a snpgds file (cf. SNPRelate package) as a Genotype
#' Data object.
#'
#' @param gds_file        Path of snpgds file
#' @param read_snp_annot  Read the SNPs' annotations
#' @param read_scan_annot Read the scans' annotations
#' @return Genotype Data object
#'
#' @examples
#' library(snplinkage)
#' gds_path <- save_hgdp_as_gds()
#' gdata <- load_gds_as_genotype_data(gds_path)
#'
#' @export
load_gds_as_genotype_data <- function(gds_file, read_snp_annot = TRUE,
  read_scan_annot = TRUE) {

  gds <- snp_gds_open(gds_file, TRUE, TRUE, TRUE)

  .read_node <- function(name) {
    if (is.null(index.gdsn(gds, name, silent = TRUE)))
      return(NULL)
    .gds_unserialize_object(gds, name)
  }

  snp_annot <- if (read_snp_annot)
      .read_node('snp_annot.serialized')
    else
      NULL

  scan_annot <- if (read_scan_annot)
      .read_node('scan_annot.serialized')
    else
      NULL

  reader <- gds_genotype_reader(gds)
  GenotypeData(reader, snp_annot, scan_annot)
}

#' gdata_snp_annots
#'
#' Get SNPs annotations from a Genotype Data object or a subset.
#'
#' @param gdata   Genotype Data object
#' @param snp_ids SNP identifiers to subset
#' @return SNP annotation data frame
#' @export
gdata_snps_annots = function(gdata, snp_ids = NULL) {
  df_snp = if (is.null(gdata@snpAnnot)) {
    data.frame(snpID = GWASTools::getSnpID(gdata),
      probe_id = GWASTools::getSnpID(gdata),
      chromosome = GWASTools::getChromosome(gdata),
      position = GWASTools::getPosition(gdata))
  } else gdata@snpAnnot@data

  if (!is.null(snp_ids)) df_snp = df_snp[match(snp_ids, df_snp$snpID), ]

  df_snp
}

#' gdata_scan_annots
#'
#' Get scans annotations from a Genotype Data object or a subset.
#'
#' @param gdata    Genotype Data object
#' @param scan_ids Scan identifiers to subset
#' @return Scans annotations data frame
#' @export
gdata_scans_annots = function(gdata, scan_ids) {
  if (is.null(gdata@scanAnnot)) {
    scan_ids = GWASTools::getScanID(gdata)
    if (!is.integer(scan_ids)) scan_ids = seq_along(scan_ids)
    data.frame(scanID = scan_ids)
  } else gdata@scanAnnot@data
}

genotype_data_subset <- function(gdata, snps_idxs, scans_idxs) {

  qced_geno <- fetch_genotypes(gdata)[snps_idxs, scans_idxs]
  df_scans = gdata_scans_annots(gdata)

  df_snps = gdata_snps_annots(gdata)
  build_gwastools(qced_geno, df_scans[scans_idxs, , drop = FALSE],
    df_snps[snps_idxs, ])
}


fetch_genotypes <- function(
  gdata,
  snps_idx = NULL,
  scans_idx = NULL,
  char = FALSE,
  snps_first = is_snp_first_dim(gdata)) {

  stopifnot(methods::is(gdata, 'GenotypeData'))

  genos <- NULL
  if (is.null(snps_idx) && is.null(scans_idx)) {
    genos <- getGenotype(gdata)
    # getGenotype may drop
    if (is.null(dim(genos))) {
      dim(genos) <- c(nsnp(gdata), nscan(gdata))
    }
  } else {
    if (is.null(snps_idx)) snps_idx <- seq_len(nsnp(gdata))
    if (is.null(scans_idx)) scans_idx <- seq_len(nscan(gdata))
    genos <- fetch_genotypes_by_idx(gdata, snps_idx, scans_idx, snps_first)
  }

  if (char) {
    if (is.null(snps_idx)) snps_idx <- seq_len(nsnp(gdata))
    if (is.null(scans_idx)) scans_idx <- seq_len(nscan(gdata))
    a1 <- fetch_allele1(gdata, snps_idx)
    a2 <- fetch_allele2(gdata, snps_idx)
    genos <- genotypeToCharacter(genos, a1, a2)
  }

  genos
}

.is_block <- function(v) {
  (utils::tail(v, 1) == v[1] + length(v) - 1L) && !is.unsorted(v)
}


.split_sorted_ints_by_blocks <- function(ints) {
  # try to be smart by testing obvious use cases first
  if (utils::tail(ints, 1) == ints[1] + length(ints) - 1) {
    return(list(c(ints[1], length(ints))))
  }

  grps <- cumsum(c(0, diff(ints) > 1))
  blks <- split(ints, grps)
  lapply(blks, function(x) c(x[1], length(x)))
}

fetch_genotypes_by_idx <- function(
  gdata,
  snps_idx = seq_len(nsnp(gdata)),
  scans_idx = seq_len(nscan(gdata)),
  snps_first = is_snp_first_dim(gdata)
) {

  snps_order <- NULL
  is_snps_block <- .is_block(snps_idx)
  snps_blocks <- if (is_snps_block) { # be smart; is it already a block ?
      list(c(snps_idx[1], length(snps_idx)))
      } else {
          snps_order <- order(snps_idx)
          .split_sorted_ints_by_blocks(snps_idx[snps_order])
      }

  scans_order <- NULL
  is_scans_block <- .is_block(scans_idx)
  scans_blocks <- if (is_scans_block) { # be smart; is it already a block ?
          list(c(scans_idx[1], length(scans_idx)))
      } else {
          scans_order <- order(scans_idx)
          .split_sorted_ints_by_blocks(scans_idx[scans_order])
      }

  genos <- fetch_genotypes_by_blocks(gdata, snps_blocks, scans_blocks, snps_first)

  if (!is_snps_block) { # need snp reordering
     if (!is_scans_block) { # need scan reordering
       genos <- genos[order(snps_order), order(scans_order), drop = FALSE]
     } else {
       genos <- genos[order(snps_order), , drop = FALSE]
     }
  } else {
    if (!is_scans_block) { # need scan reordering
      genos <- genos[, order(scans_order), drop = FALSE]
    }
  }

  genos
}


fetch_genotypes_by_blocks <- function(
  gdata,
  snps_blocks = list(c(1, nsnp(gdata))),
  scans_blocks = list(c(1, nscan(gdata))),
  snps_first = is_snp_first_dim(gdata)) {

  ### special case: one snp block - one scan block
  if (length(snps_blocks) == 1 && length(scans_blocks) == 1) {
    snp = snps_blocks[[1]]
    scan = scans_blocks[[1]]
    m <- getGenotype(gdata, snp = snp, scan = scan)
    if (is.null(dim(m))) dim(m) <- c(snp[2], scan[2])
    return(m)
  }

  nb_snps <- sum(sapply(snps_blocks, '[[', 2))
  nb_scans <- sum(sapply(scans_blocks, '[[', 2))

  # preallocate results matrix
  genotypes <- matrix(0L, nb_snps, nb_scans)
  if (snps_first) {
    col <- 1L
    for (scan in scans_blocks) {
      y <- seq.int(from = col, length.out = scan[2])
      col <- col + scan[2]
      row <- 1L
      for (snp in snps_blocks) {
        m <- getGenotype(gdata, snp = snp, scan = scan)
        if (is.null(dim(m))) dim(m) <- c(snp[2], scan[2])
        x <- seq.int(from = row, length.out = snp[2])
        row <- row + snp[2]
        genotypes[x, y] <- m
      }
    }
  } else {
    row <- 1L
    for (snp in snps_blocks) {
      x <- seq.int(from = row, length.out = snp[2])
      row <- row + snp[2]
      col <- 1L
      for (scan in scans_blocks) {
        m <- getGenotype(gdata, snp = snp, scan = scan)
        if (is.null(dim(m))) {
          dim(m) <- c(snp[2], scan[2])
        }
        y <- seq.int(from = col, length.out = scan[2])
        col <- col + scan[2]
        genotypes[x, y] <- m
      }
    }
  }

  genotypes
}

save_genotype_data_as_gds <- function(
  gdata,
  filename,
  save_annotations = TRUE,
  compress = "ZIP.max",
  quiet = FALSE,
  chunk_size = 1000,
  ...
) {
# code adapted from GWASTools::convertNcdfGds

  info <- if (quiet) function(...) {} else function(...) message(sprintf(...))

  ### snp gds

  tt <- system.time({
    allele1 <- fetch_allele1(gdata) %//% 'A'
    allele2 <- fetch_allele2(gdata) %//% 'B'

    alleles <- paste(allele1, allele2, sep = '/')
    alleles[is.na(allele1) | is.na(allele2)] <- ''
  }, gcFirst = FALSE)
  info('fetched alleles in %.2fs', tt[3])


  tt <- system.time({
    snp_infos <- data.frame(
      snpID = getSnpID(gdata),
      chromosome = getChromosome(gdata),
      position = getPosition(gdata),
      alleles = alleles
      , stringsAsFactors = FALSE)
  }, gcFirst = FALSE)
  info('fetched snp_infos in %.2fs', tt[3])

  info('saving genotypes...')
  tt <- system.time({
    it <- fetch_genotypes_by_chunk_iterator(gdata, chunk_size = chunk_size, ...)
    write_snp_gds(filename, it,
      sample_ids = getScanID(gdata), snp_infos = snp_infos, quiet = quiet)
  }, gcFirst = FALSE)
  info('took %s', tt[3])

  ### annotations
  if (save_annotations) {
    tt <- system.time({
      gds <- snp_gds_open(filename, FALSE)
      .gds_serialize_object(gds, 'snp_annot.serialized', get_snp_annot(gdata))
      .gds_serialize_object(gds, 'scan_annot.serialized', get_scan_annot(gdata))
      closefn.gds(gds)
    }, gcFirst = FALSE)
    info('saved annotations in %.2fs', tt[3])
  }

  invisible()
}


.gds_serialize_object <- function(gds, name, object, compress = "ZIP.max") {
  ser <- serialize(object, NULL)
  add.gdsn(gds, name, as.integer(ser), storage = "uint8", compress = compress)
}


.gds_unserialize_object <- function(gds, name) {
  ser <- read.gdsn(index.gdsn(gds, name))
  unserialize(as.raw(ser))
}



write_genotypes_to_gds <- function(
  gds,
  geno_iterator,
  node_name = 'genotype',
  quiet = FALSE)
{
  info <- if (quiet) function(...) {} else function(...) message(sprintf(...))

  .convert_chunk <- function(x) {
    x[is.na(x)] <- 3L
    t(x)
  }

  nb_chunks <- geno_iterator$nb

  ### save the first chunk
  chunk <- .convert_chunk(geno_iterator$get(1))
  node <- add.gdsn(gds, node_name, chunk, storage = "bit2")

  current <- 0L
  for (i in seq.int(from = 2, length.out = nb_chunks - 1)) {
    perc <- ceiling(i / nb_chunks * 100)
    if (perc > current) {
      info('writing genotypes %i%%', perc)
      current <- perc
    }

    chunk <- .convert_chunk(geno_iterator$get(i))
    append.gdsn(node, chunk)
  }

  put.attr.gdsn(node, "sample.order")

  node
}


# workaround for current bug when forking in gdsfmt
.create_gds <- function(filename) {
  gds <- createfn.gds(filename)
  closefn.gds(gds)
  openfn.gds(filename, readonly = FALSE, allow.duplicate = FALSE,
    allow.fork = TRUE)
}

write_snp_gds <- function(filename, geno_iterator, sample_ids, snp_infos,
  compress = "ZIP.max", quiet = FALSE, ...) {

  stopifnot(!anyDuplicated(sample_ids))
  cols <- c('snpID', 'chromosome', 'position', 'alleles')
  stopifnot(all(cols %in% colnames(snp_infos)), !anyDuplicated(snp_infos$snpID))

  info <- if (quiet) function(...) {} else function(...) message(sprintf(...))

  gfile <- .create_gds(filename)
  on.exit(closefn.gds(gfile), TRUE)

  .add <- function(name, val) {
    info('writing %s...', name)
    add.gdsn(gfile, name, val, compress = compress, closezip = TRUE)
  }

  .add("sample.id", sample_ids)
  .add("snp.id", snp_infos$snpID)
  .add("snp.position", snp_infos$position)
  .add("snp.chromosome", snp_infos$chromosome)
  .add("snp.allele", snp_infos$alleles)

  write_genotypes_to_gds(gfile, geno_iterator, quiet = quiet, ...)
}




# temp wrapper for snpgdsOpen while it is not mainstream
snp_gds_open <- function(...) {
  SNPRelate <- NULL; Library <- library; Library(SNPRelate)
  gds <- SNPRelate::snpgdsOpen(...)

  # hack: some GWASTools functions check for gds.class
  class(gds) <- "gds.class"

  gds
}



request_snpgds_file <- function(
    gdata,
    snps_idx = NULL,
    scans_idx = NULL
) {
  gds <- fetch_gds(gdata)
  new_file <- NULL
  if (is.null(gds)) {
    # no gds, must create one
    gdsfn <- tempfile(fileext = '.gds')

    # subset if needed
    gdata2 <- if (is.null(snps_idx) && is.null(scans_idx)) {
        gdata
      } else {
        if (is.null(snps_idx)) snps_idx = seq_len(nsnp(gdata))
        if (is.null(scans_idx)) scans_idx = seq_len(nscan(gdata))
        genotype_data_subset(gdata, snps_idx, scans_idx)
      }

    save_genotype_data_as_gds(gdata2, gdsfn, quiet = TRUE)
    gds <- snp_gds_open(gdsfn, TRUE, TRUE, TRUE)
    new_file <- TRUE
    # gdata already subsetted now
    snps_idx <- scans_idx <- NULL
    snp_ids <- scan_ids <- NULL

  } else {
    gdsfn <- gds$filename
    new_file <- FALSE
    snp_ids <- if (is.null(snps_idx)) NULL else {
              getSnpID(gdata, snps_idx)
    }
    scan_ids <- if (is.null(scans_idx)) NULL else {
              getScanID(gdata, scans_idx)
          }

  }

  list(snpgds = gds, snps_idx = snps_idx, scans_idx = scans_idx,
      new_file = new_file, snp_ids = snp_ids, scan_ids = scan_ids)
}


gds_genotype_reader <- function(gds, ...) {
  if (is.character(gds)) {
    gds <- snp_gds_open(gds, TRUE, TRUE, TRUE)
  }

  genotypeDim <- if (is_snp_first_dim(gds)) "snp,scan" else "scan,snp"
  gdsreader <- methods::new('GdsReader', filename = gds$filename, handler = gds)

  methods::new('GdsGenotypeReader', gdsreader, genotypeDim = genotypeDim, ...)
}


fetch_genotypes_by_chunk_iterator <- function(
  gdata,
  chunk_size = 100L,
  nb_chunks = NULL) {

  nb_snps <- nsnp(gdata)

  if (is.null(nb_chunks)) {
    nb_chunks <- as.integer(ceiling(nb_snps / chunk_size))
  } else {
    chunk_size <- as.integer(ceiling(nb_snps / nb_chunks))
  }

  .get_chunk <- function(n) {
    die_unless(n > 0 && n <= nb_chunks, 'bad chunk index %i', n)
    i <- (n - 1) * chunk_size + 1L
    j <- min(nb_snps, i + chunk_size - 1L)
    idx <- seq.int(i, j)
    fetch_genotypes_by_idx(gdata, idx)
  }

  list(nb = nb_chunks, chunk_size = chunk_size, get = .get_chunk)
}



read_snp_gds_alleles <- function(gds) {
  if (is.null(index.gdsn(gds, "snp.allele", silent = TRUE)))
    return(NULL)

  alleles <- read.gdsn(index.gdsn(gds, "snp.allele"))
  alleles[!nzchar(alleles)] <- NA

  alleles
}



