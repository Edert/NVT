#global methods availible
method_v <- c("N","TC","Med","TMM","UQ","Q","RPKM","DEQ","ER","TPM","G")

#initialize and load input
init <- function(hkgene_list, exp_list1, exp_list2, method){

  if (missing('hkgene_list')){
    stop("No housekeeping-gene list specified!")
  }
  if (missing('exp_list1')){
    stop("No expression list specified!")
  }
  if (missing('exp_list2')){
    stop("No second expression list specified!")
  }
  if (missing('method')){
    stop("No method specified!")
  }

  if(check_expression_list(exp_list1) && check_expression_list(exp_list2) && check_hkgene_list(hkgene_list) && check_method(method)) {
    NVTData = return(new("NVTdata",exp1=exp_list1,exp2=exp_list2,hklist=hkgene_list,method=method))
    return(NVTData)
  }
}

check_hkgene_list <- function(hkgene_list) {
  if( is.character(hkgene_list) ){
    return(TRUE)
  }else{
    stop("Housekeeping-gene list is no list!")
  }
}

check_expression_list <- function(exp_list) {
  if( is.data.frame(exp_list) ){
    return(TRUE)
  }else{
    stop("Expression list is no data frame!")
  }
}

check_method <- function(method) {
  if( method %in% method_v ){
    return(TRUE)
  }else{
    stop("Unknown method specified!")
  }
}
