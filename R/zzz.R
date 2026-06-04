.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "saguiSED is an optional SED-fitting extension. ",
    "Use check_bagpipes() before running the Bagpipes backend."
  )
}
