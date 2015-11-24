#global methods availible
method_v <- c("N","TC","Med","TMM","UQ","Q","RPKM","DEQ","ER","TPM","G")

#initialize and load input
NVTinit <- function(hkgene_list, exp_list1, exp_list2, method){

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

#normalize
NVTnormalize <- function(NVTdata) {

  if(check_expression_list(NVTdata@exp1) && check_expression_list(NVTdata@exp2) && check_hkgene_list(NVTdata@hklist) && check_method(NVTdata@method)){
    switch(NVTdata@method,
           N={
             #NVTdata@norm1 <- NVTdata@exp1
             #NVTdata@norm2 <- NVTdata@exp2
           },
           TC={

           },
           Med={

           },
           TMM={

           },
           UQ={

           },
           Q={

           },
           RPKM={

           },
           DEQ={

           },
           ER={

           },
           TPM={

           },
           G={

           },
           {
             NVTdata@norm1 <- NVTdata@exp1
             NVTdata@norm2 <- NVTdata@exp2
           }
    )

    return(NVTdata)
  }else{
    stop("Not a valid NVTdata object!")
  }
}

#plot data
NVTplot <- function(NVTdata) {
  if(check_expression_list(NVTdata@exp1) && check_expression_list(NVTdata@exp2)
     && check_hkgene_list(NVTdata@hklist) && check_method(NVTdata@method)
     && check_expression_list(NVTdata@norm1)  && check_expression_list(NVTdata@norm2)){

    plot(log(NVTdata@norm1[,1]),log(NVTdata@norm2[,1]),main=paste("MA-plot", names(mynorm@norm1),"vs.",names(mynorm@norm2)),
         xlab=paste("log( normalized expression",names(mynorm@norm1),")"),ylab=paste("log( normalized expression",names(mynorm@norm2),")")
         ,asp=1,pch=20,col="grey")

    m1 <- log(NVTdata@norm1[NVTdata@hklist,])
    m2 <- log(NVTdata@norm2[NVTdata@hklist,])

    points(m1,m2,col="blue",pch=19)

    #fm <- lm(m2[,1] ~ m1[,1])
    #abline(fm, col = "red")

  }else{
    stop("Not a valid NVTdata object with normalized values!")
  }
}


#calculate pearson correclation of housekeeping genes
NVTpearson <- function(NVTdata) {
  if(check_expression_list(NVTdata@exp1) && check_expression_list(NVTdata@exp2)
     && check_hkgene_list(NVTdata@hklist) && check_method(NVTdata@method)
     && check_expression_list(NVTdata@norm1)  && check_expression_list(NVTdata@norm2)){


  }else{
    stop("Not a valid NVTdata object with normalized values!")
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
