
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
##' form \code{new("NVTdata", ...)}.
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
    norm_method="character",
    norm_method_name="character",
    is_norm="logical",
    exp1="data.frame",
    exp2="data.frame",
    hklist="character",
    norm1="data.frame",
    norm2="data.frame",
    length="data.frame"
    ),
  prototype=prototype(
    norm_method="N",
    norm_method_name="Not",
    is_norm=FALSE,
    exp1=data.frame(),
    exp2=data.frame(),
    hklist=character(),
    norm1=data.frame(),
    norm2=data.frame(),
    length=data.frame()
  )
)
