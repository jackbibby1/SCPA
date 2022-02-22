.onAttach <- function(libname, pkgname) {
  version <- packageDescription(pkgname, fields = "Version")
  msg <- paste0("====================
", "This is ", "version ", version, " of ", pkgname, "
SCPA GitHub page: https://github.com/jackbibby1/SCPA
For SCPA tutorials: https://jackbibby1.github.io/SCPA/
If you use this, please cite: https://www.biorxiv.org/content/10.1101/2022.02.07.478807v1
====================
")

  packageStartupMessage(msg)

}
