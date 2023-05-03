from rocker/verse:4.3.0
run apt-get update && \
  apt-get install -y --no-install-recommends texlive texlive-latex-recommended texlive-fonts-extra qpdf && \
  R -e "install.packages(c('pkgbuild', 'roxygen2', 'testthat'))"
run R -e "install.packages(c('cowplot', 'gdsfmt', 'ggplot2', 'ggrepel', 'gtable', 'GWASTools', 'biomaRt', 'data.table', 'dplyr', 'knitr', 'plyr', 'SNPRelate'))"
add ./DESCRIPTION /snplinkage/DESCRIPTION
run R -e "devtools::install_deps('snplinkage', dependencies = TRUE)"
add ./ /snplinkage
run R -e "devtools::install('snplinkage', dependencies = TRUE)"
