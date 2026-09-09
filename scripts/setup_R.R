dir.create('env/R_library',recursive=TRUE,showWarnings=FALSE)
.libPaths(c(normalizePath('env/R_library'),.libPaths()))
options(timeout=300,repos=c(CRAN='https://cloud.r-project.org'))
if(!requireNamespace('BiocManager',quietly=TRUE)) install.packages('BiocManager',lib='env/R_library')
for(p in c('limma','statmod','jsonlite')) {
  if(!requireNamespace(p,quietly=TRUE)) BiocManager::install(p,lib='env/R_library',ask=FALSE,update=FALSE)
}
