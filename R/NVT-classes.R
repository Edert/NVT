
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

#'Show method, to retrieve the normalized expression values
#'
#'@docType methods
#'@name show
#'@rdname show
#'@aliases show,NVTdata-method
#'@param object A previously initialized and normalized NVTobject
#'@return Returns the normalized expression values of both samples
#'@examples
#'library("NVT")
#'data(myexp1)
#'data(myexp2)
#'data(mylen)
#'mylist1<-c("ENSG00000111640","ENSG00000163631","ENSG00000075624","ENSG00000172053",
#'"ENSG00000170950","ENSG00000165704","ENSG00000196839","ENSG00000168938","ENSG00000177700")
#'
#'mynvt <- NVTinit(mylist1,myexp1,myexp2,"N")
#'mynorm <- NVTnormalize(mynvt)
#'
#'show(mynorm)
#'@export
setMethod(f = "show", signature = "NVTdata",
          definition = function(object) {
            if(check_norm_list(object@norm1)  && check_norm_list(object@norm2)){
              myout <- cbind(object@norm1,object@norm2)
              return(myout)
            }else{
              stop("Not a valid NVTdata object with normalized values!")
            }
          }
)
