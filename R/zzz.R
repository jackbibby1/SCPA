.onAttach <- function(libname, pkgname) {
  version <- utils::packageDescription(pkgname, fields = "Version")
  msg <- paste0("====================
", "This is ", "version ", version, " of ", pkgname, "
For SCPA tutorials and latest version: https://jackbibby1.github.io/SCPA/
For the SCPA GitHub page: https://github.com/jackbibby1/SCPA
If you use SCPA, please cite: Bibby JA. et al. Cell Rep. 2022
====================
")

  packageStartupMessage(msg)

}
