from rocker/verse:4.4.0
run apt-get update && \
  apt-get install -y --no-install-recommends texlive texlive-latex-recommended texlive-fonts-extra qpdf && \
  R -e "install.packages(c('pkgbuild', 'roxygen2', 'testthat'))"

run R -e 'install.packages("BiocManager")'
run R -e 'BiocManager::install("S4Vectors")'
run R -e 'BiocManager::install(c("GWASTools", "gdsfmt", "SNPRelate", "biomaRt"))'

run R -e "install.packages(c('cowplot', 'ggplot2', 'ggrepel', 'gtable', 'data.table', 'dplyr', 'knitr', 'plyr'))"
add ./DESCRIPTION /snplinkage/DESCRIPTION
run R -e "devtools::install_deps('snplinkage', dependencies = TRUE)"
add ./ /snplinkage
run R -e "devtools::install('snplinkage', dependencies = TRUE)"
