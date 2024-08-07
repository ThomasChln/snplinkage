# SNPLinkage

[![CRAN version](https://www.r-pkg.org/badges/version/snplinkage)](https://cran.r-project.org/package=snplinkage)
[![CRAN total downloads](https://cranlogs.r-pkg.org/badges/grand-total/snplinkage)](https://cran.r-project.org/package=snplinkage)
[![CRAN monthly downloads](https://cranlogs.r-pkg.org/badges/snplinkage)](https://cran.r-project.org/package=snplinkage)
[![Coverage status](https://codecov.io/github/ThomasChln/snplinkage/branch/main/graph/badge.svg?token=DFWQHUXPNE)](https://codecov.io/github/thomaschln/snplinkage)

This R package provides linkage disequilibrium visualizations by displaying correlation matrices annotated with chromosomic positions and gene names. Two types of displays are provided to focus on small or large regions, and both can be extended to combine associations results or investigate feature selection methods.

## Summary

  Physically close SNPs exhibit correlation structures due to biological, hereditary and evolutionary factors. When they are more correlated than randomly expected they are said to be in linkage disequilibrium and groups of SNPs in strong linkage disequilibrium are called haplotypes. Genetic association studies often investigate correlation structures to explore specific regions associated with a trait or disease, and genome-wide association studies usually remove highly correlated SNPs to increase statistical power using TagSNP selection.

  This package provides linkage disequilibrium visualizations that combine the correlation matrix of SNPs with their chromosomic positions and gene names in order to investigate correlation structures of specific regions. Two types of displays are available to focus on small or large regions, and both can be extended to combine association results or to investigate feature selection methods. The correlations are computed using the SNPRelate package and the plots are customizable ggplot2 and gtable objects and are annotated using the biomaRt package.

## Installation

Stable CRAN version, in R:

```r
  install.packages('snplinkage')
```

Development version, using the devtools package in R:

```r
  devtools::install_git('https://gitlab.com/thomaschln/snplinkage.git')
```

## Usage

Open Genotype Data (GWASTools package) file and perform quality control and feature selection.

```r
  library('snplinkage')
  gds_path <- save_hgdp_as_gds()
  gdata <- load_gds_as_genotype_data(gds_path)
  qc <- snprelate_qc(gdata, tagsnp = .99)
```

Select a small region (20 SNPs) and visualize

```r
  snp_idxs_1p13 <- select_region_idxs(qc$gdata,
    chromosome = 1, position_min = 114.4e6, n_snps = 20, offset = 12)

  plt <- gtable_ld_gdata(qc$gdata, snp_idxs_1p13, labels_colname = 'probe_id')
  grid::grid.draw(plt)
```

Select a large region (300 SNPs) and visualize

```r
  snp_idxs_8p23 <- select_region_idxs(qc$gdata, chromosome = 8,
    position_min = 11e6, position_max = 12e6)

  df_ld <- snprelate_ld(qc$gdata, snps_idx = snp_idxs_8p23, quiet = TRUE)
  plt <- gtable_ld(df_ld, df_snp = gdata_snps_annots(qc$gdata))
  grid::grid.draw(plt)
```

Add gene annotations using biomaRt

```r
  snp_idxs_hladr <- select_region_idxs(qc$gdata,
    chromosome = 6, position_min = 32.5e6, n_snps = 20, offset = 9)

  # qc$gdata <- gdata_add_gene_annots(qc$gdata, snp_idxs_hladr)
  qc$gdata <- gdata_add_gene_annots_hladr_example(qc$gdata, snp_idxs_hladr)

  plt <- gtable_ld_gdata(qc$gdata, snp_idxs_hladr, labels_colname = 'gene')
  grid::grid.draw(plt)
```

Combine with association studies results

* Perform basic geographical association study and fetch gene annotations

```r
  snp_idxs_mhc <- select_region_idxs(qc$gdata,
    chromosome = 6, position_min = 29e6, position_max = 33e6)
  df_assocs <- chisq_pvalues_gdata(qc$gdata, snp_idxs_mhc)

  df_top_aim <- subset(df_assocs, rank(-pvalues, ties.method = 'first') <= 20)

  #qc$gdata <- gdata_add_gene_annots(qc$gdata, rownames(df_top_aim))
  qc$gdata <- gdata_add_gene_annots_aim_example(qc$gdata, rownames(df_top_aim))
```

* Visualize only the 20 most associated SNPs

```r
  plt <- gtable_ld_associations_gdata(df_top_aim, qc$gdata,
    labels_colname = 'gene')
  grid::grid.draw(plt)
```

* Visualize the whole region

```r
  plt <- gtable_ld_associations_gdata(df_assocs, qc$gdata,
    labels_colname = 'gene')
  grid::grid.draw(plt)
```

Visualize TagSNP feature selection biplot

* Using a small region

```r
  plt <- gtable_ld_gdata(qc$gdata, snp_idxs_1p13, r2 = 0.8)
  grid::grid.draw(plt)
```

* Using a large region

```r
  plt <- gtable_ld_gdata(qc$gdata, snp_idxs_1p13_large, r2 = 0.8)
  grid::grid.draw(plt)
```

See the [vignette](https://cran.r-project.org/web/packages/snplinkage/vignettes/snplinkage.pdf) for results and further details.

## Further reading

The code in this package was used to perform large-scale visualizations of up to 500 SNPs of three genomic regions (*HLA*, *1p13.2*, *8p23*) in a population of 1,200 systemic autoimmune diseases patients and healthy controls, sampled by the European research project PreciseSADs, in my [Ph.D. thesis](https://archive-ouverte.unige.ch/unige:161795).

## Acknowledgements

The package uses code first published in the [SNPClust package](https://github.com/ThomasChln/snpclust), which was co-authored with Karl Forner, Alessandro Di Cara and Jérôme Wojcik.

## License

This package is free and open source software, licensed under GPL-3.
