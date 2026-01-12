pkgs <- c('Seurat', 'Matrix', 'dplyr', 'ggplot2', 'qvalue', 'GEOquery', 'testthat', 'tidyr')
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    cat(sprintf("Installing %s...\n", p))
    install.packages(p, repos='https://cloud.r-project.org')
  } else {
    cat(sprintf("%s: OK\n", p))
  }
}
