.onLoad <- function(libname, pkgname) {
  # set global variables in order to avoid CHECK notes
  utils::globalVariables("type")
  utils::globalVariables("seqnames")
  utils::globalVariables("queryHits")
  utils::globalVariables("width")
  utils::globalVariables("strand")
  utils::globalVariables("label")
  utils::globalVariables("end")
  utils::globalVariables("gene_name")
  utils::globalVariables("Type")

  invisible()
}
