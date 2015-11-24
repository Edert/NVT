
###==========================
### NVT
###==========================
##' Class "NVTdata"
##'
##' A class with all required NVTdata
##'
##'
##' @name NVTdata class
##' @rdname NVT-class
##' @aliases NVT-class NVT
##' @docType class
##' @section Objects from the Class: Objects can be created by calls of the
##' form \code{new("NVT", ...)}.
##' @author Thomas Eder
##' @keywords classes
##' @examples
##'
##' showClass("NVTdata")
##'
##' @exportClass NVTdata
setClass(
  Class="NVTdata",
  representation=representation(
    method="character",
    exp1="data.frame",
    exp2="data.frame",
    hklist="character")
)
