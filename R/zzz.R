.onAttach <- function(libname, pkgname) {
  version <- utils::packageDescription(pkgname, fields = "Version")
  msg <- paste0("====================
", "This is ", "version ", version, " of ", pkgname, "
SCPA GitHub page: https://github.com/jackbibby1/SCPA
For SCPA tutorials: https://jackbibby1.github.io/SCPA/
====================
")

  packageStartupMessage(msg)

}
