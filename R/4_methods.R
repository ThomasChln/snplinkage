

#' Fetch GDS (default)
#'
#' @param obj Default object
#' @param ... Not passed
#' @return NULL
#' @export
fetch_gds.default <- function(obj, ...) {
  NULL
}

#' Fetch GDS (GdsGenotypeReader)
#'
#' @param obj GdsGenotypeReader object
#' @param ... Not passed
#' @return S4 slot 'handler' of obj
#' @export
fetch_gds.GdsGenotypeReader <- function(obj, ...) {
  obj@handler
}

#' Fetch GDS (GenotypeData)
#'
#' @param obj GenotypeData object
#' @param ... Not passed
#' @return fetch_gds output on S4 slot 'data' of obj
#' @export
fetch_gds.GenotypeData <- function(obj, ...) {
  fetch_gds(obj@data)
}


#' Fetch GDS (GenotypeDataSubset)
#'
#' @param obj GenotypeDataSubset object
#' @param ... Not passed
#' @return NULL
#' @export
fetch_gds.GenotypeDataSubset <- function(obj, ...) {
  NULL
}


#' Is SNP first dimension (default)
#'
#' @param obj Default object
#' @param ... Not passed
#' @return NA
#' @export
is_snp_first_dim.default <- function(obj, ...) {
  NA
}

#' Is SNP first dimension (GDS object)
#'
#' @param obj GDS object
#' @param ... Not passed
#' @return Logical, TRUE if SNP is first dimension
#' @export
is_snp_first_dim.gds.class <- function(obj, ...) {
  snpfirstdim <- NA
  rd <- names(get.attr.gdsn(index.gdsn(obj, "genotype")))
	if ("snp.order" %in% rd) snpfirstdim <- TRUE
	if ("sample.order" %in% rd) snpfirstdim <- FALSE

  snpfirstdim
}


#' Is SNP first dimension (GenotypeData object)
#'
#' @param obj Genotype data object
#' @param ... Not passed
#' @return is_snp_first_dim output on S4 slot 'data'
#' @export
is_snp_first_dim.GenotypeData <- function(obj, ...) {
  is_snp_first_dim(obj@data)
}



#' Is SNP first dimension (GdsGenotypeReader object)
#'
#' @param obj GdsGenotypeReader object
#' @param ... Not passed
#' @return is_snp_first_dim output on S4 slot 'handler'
#' @export
is_snp_first_dim.GdsGenotypeReader <- function(obj, ...) {
  is_snp_first_dim(obj@handler)
}



#' Is SNP first dimension (MatrixGenotypeReader object)
#'
#' @param obj MatrixGenotypeReader object
#' @param ... Not passed
#' @return TRUE
#' @export
is_snp_first_dim.MatrixGenotypeReader <- function(obj, ...) {
  TRUE
}


#' Is SNP first dimension (NcdfGenotypeReader object)
#'
#' @param obj NcdfGenotypeReader object
#' @param ... Not passed
#' @return TRUE
#' @export
is_snp_first_dim.NcdfGenotypeReader <- function(obj, ...) {
  TRUE
}








#' Fetch allele 1 (default object)
#'
#' @param obj Default object
#' @param snps_idx SNPs indexes
#' @return NULL 
#' @export
fetch_allele1.default <- function(obj, snps_idx) { NULL }


#' Fetch allele 2 (default object)
#'
#' @param obj Default object
#' @param snps_idx SNPs indexes
#' @return NULL 
#' @export
fetch_allele2.default <- function(obj, snps_idx) { NULL }


#' Fetch allele 1 (GenotypeData object)
#'
#' @param obj GenotypeData object
#' @param ... Passed to getAlleleA
#' @return Allele 1 
#' @export
fetch_allele1.GenotypeData <- function(obj, ...) {
  alleles <- getAlleleA(obj, ...)
  if (!is.null(alleles)) return(alleles)

  fetch_allele1(obj@data, ...)
}


#' Fetch allele 2 (GenotypeData object)
#'
#' @param obj GenotypeData object
#' @param ... Passed to getAlleleB
#' @return Allele 2 
#' @export
fetch_allele2.GenotypeData <- function(obj, ...) {
  alleles <- getAlleleB(obj, ...)
  if (!is.null(alleles)) return(alleles)

  fetch_allele2(obj@data, ...)
}


fetch_snpgds_alleles <- function(gdsobj, snps_idx) {
  node <- index.gdsn(gdsobj, "snp.allele", silent = TRUE)
  if (is.null(node)) return(NULL)

  if (missing(snps_idx)) {
    return(read.gdsn(node))
  }

  sel <- rep(FALSE, objdesp.gdsn(node)$dim)
  sel[snps_idx] <- TRUE
  alleles <- readex.gdsn(node, sel)
  alleles[order(snps_idx)]
}


#' Fetch allele 1 (GdsGenotypeReader object)
#'
#' @param obj GenotypeData object
#' @param snps_idx SNPs indexes
#' @return Allele 1 
#' @export
fetch_allele1.GdsGenotypeReader <- function(obj, snps_idx) {
  alleles <- fetch_snpgds_alleles(obj@handler, snps_idx)
  if (is.null(alleles)) return(NULL)
  substr(alleles, 1, 1)
}


#' Fetch allele 2 (GdsGenotypeReader object)
#'
#' @param obj GenotypeData object
#' @param snps_idx SNPs indexes
#' @return Allele 2 
#' @export
fetch_allele2.GdsGenotypeReader <- function(obj, snps_idx) {
  alleles <- fetch_snpgds_alleles(obj@handler, snps_idx)
  if (is.null(alleles)) return(NULL)
  substr(alleles, 3, 3)
}




#' Get scans annotations (GenotypeData object)
#'
#' @param obj GenotypeData object
#' @param ... Not passed
#' @return Data frame
#' @export
get_scan_annot.GenotypeData <- function(obj, ...) {
  obj@scanAnnot
}


#' Get SNPs annotations (GenotypeData object)
#'
#' @param obj GenotypeData object
#' @param ... Not passed
#' @return Data frame
#' @export
get_snp_annot.GenotypeData <- function(obj, ...) {
  obj@snpAnnot
}


#' Get scans annotations (GenotypeDataSubset object)
#'
#' @param obj GenotypeDataSubset object
#' @param ... Not passed
#' @return Data frame
#' @export
get_scan_annot.GenotypeDataSubset <- function(obj, ...) {
  obj@scanAnnot[obj@scans_idx, ]
}


#' Get SNPs annotations (GenotypeDataSubset object)
#'
#' @param obj GenotypeDataSubset object
#' @param ... Not passed
#' @return Data frame
#' @export
get_snp_annot.GenotypeDataSubset <- function(obj, ...) {
  obj@snpAnnot[obj@snps_idx, ]
}

