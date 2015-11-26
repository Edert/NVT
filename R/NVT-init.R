#global methods availible
method_v <- c("N","TC","Med","TMM","UQ","UQ2","Q","RPKM","RPKM2","RPM","DEQ","TPM","G")

#initialize and load input
NVTinit <- function(hkgene_list, exp_list1, exp_list2, method, length){

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
  #length missing
  if (missing('length')){
    if(check_expression_list(exp_list1) && check_expression_list(exp_list2) && check_hkgene_list(hkgene_list) && check_method(method)) {
      NVTData = return(new("NVTdata",exp1=exp_list1,exp2=exp_list2,hklist=hkgene_list,method=method))
      return(NVTData)
    }
  }

  #all here add length
  if(check_expression_list(exp_list1) && check_expression_list(exp_list2) && check_hkgene_list(hkgene_list) && check_method(method)) {
    NVTData = return(new("NVTdata",exp1=exp_list1,exp2=exp_list2,hklist=hkgene_list,method=method,length=length))
    return(NVTData)
  }
}

#normalize
NVTnormalize <- function(NVTdataobj) {

  if(check_expression_list(NVTdataobj@exp1) && check_expression_list(NVTdataobj@exp2) && check_hkgene_list(NVTdataobj@hklist) && check_method(NVTdataobj@method)){
    switch(NVTdataobj@method,
           N={
             print ("No normalization!")

             NVTdataobj@norm1 <- NVTdataobj@exp1
             NVTdataobj@norm2 <- NVTdataobj@exp2
           },
           TC={
             print ("Total count normalization!")

             ci1 <- NVTdataobj@exp1[,1]
             N1 <- sum(ci1)
             ci2 <- NVTdataobj@exp2[,1]
             N2 <- sum(ci2)
             m=mean(c(N1,N2))
             NVTdataobj@norm1 <-  as.data.frame(ci1/(N1+N2)*m)
             NVTdataobj@norm2 <-  as.data.frame(ci2/(N1+N2)*m)
             colnames(NVTdataobj@norm1)=names(NVTdataobj@exp1)
             colnames(NVTdataobj@norm2)=names(NVTdataobj@exp2)
             rownames(NVTdataobj@norm1)=rownames(NVTdataobj@exp1)
             rownames(NVTdataobj@norm2)=rownames(NVTdataobj@exp2)
           },
           Med={
             print ("Median normalization!")

             ci1 <- NVTdataobj@exp1[,1]
             N1 <- median(ci1)
             ci2 <- NVTdataobj@exp2[,1]
             N2 <- median(ci2)
             m=mean(c(N1,N2))
             NVTdataobj@norm1 <-  as.data.frame(ci1/(N1+N2)*m)
             NVTdataobj@norm2 <-  as.data.frame(ci2/(N1+N2)*m)
             colnames(NVTdataobj@norm1)=names(NVTdataobj@exp1)
             colnames(NVTdataobj@norm2)=names(NVTdataobj@exp2)
             rownames(NVTdataobj@norm1)=rownames(NVTdataobj@exp1)
             rownames(NVTdataobj@norm2)=rownames(NVTdataobj@exp2)
           },
           TMM={
             print ("Trimmed Mean of M-values normalization!")

             library("NOISeq")
             mymatrix <- as.matrix(cbind(NVTdataobj@exp1,NVTdataobj@exp2))

             mynormmatrix <- tmm(mymatrix)
             NVTdataobj@norm1 <- as.data.frame(mynormmatrix[,1])
             NVTdataobj@norm2 <- as.data.frame(mynormmatrix[,2])
             colnames(NVTdataobj@norm1)=names(NVTdataobj@exp1)
             colnames(NVTdataobj@norm2)=names(NVTdataobj@exp2)
           },
           UQ={
             print ("Upper Quartile normalization!")

             ci1 <- NVTdataobj@exp1[,1]
             N1 <- quantile(ci1, 0.75)
             ci2 <- NVTdataobj@exp2[,1]
             N2 <- quantile(ci2, 0.75)
             m=mean(c(N1,N2))
             NVTdataobj@norm1 <- as.data.frame(ci1/(N1+N2)*m)
             NVTdataobj@norm2 <- as.data.frame(ci2/(N1+N2)*m)
             colnames(NVTdataobj@norm1)=names(NVTdataobj@exp1)
             colnames(NVTdataobj@norm2)=names(NVTdataobj@exp2)
             rownames(NVTdataobj@norm1)=rownames(NVTdataobj@exp1)
             rownames(NVTdataobj@norm2)=rownames(NVTdataobj@exp2)
           },
           UQ2={
             print ("Upper Quartile normalization!")

             library("NOISeq")
             mymatrix <- as.matrix(cbind(NVTdataobj@exp1,NVTdataobj@exp2))

             mynormmatrix <- uqua(mymatrix)
             NVTdataobj@norm1 <- as.data.frame(mynormmatrix[,1])
             NVTdataobj@norm2 <- as.data.frame(mynormmatrix[,2])
             colnames(NVTdataobj@norm1)=names(NVTdataobj@exp1)
             colnames(NVTdataobj@norm2)=names(NVTdataobj@exp2)
           },
           Q={
             print ("Quantile normalization!")
             library("limma")

             mymatrix <- as.matrix(cbind(NVTdataobj@exp1,NVTdataobj@exp2))
             mynormmatrix <- normalizeBetweenArrays(mymatrix,"quantile")
             NVTdataobj@norm1 <- as.data.frame(mynormmatrix[,1])
             NVTdataobj@norm2 <- as.data.frame(mynormmatrix[,2])
             colnames(NVTdataobj@norm1)=names(mynormmatrix[1,])[1]
             colnames(NVTdataobj@norm2)=names(mynormmatrix[1,])[2]
           },
           RPKM={
             print ("RPKM normalization!")

             if(nrow(NVTdataobj@length)==0){
               stop("No RPKM normalization possible without gene length")
             }

             li <- NVTdataobj@length

             ci1 <- NVTdataobj@exp1[,1]
             N1 <- sum(ci1)
             NVTdataobj@norm1 <- ci1/(N1/10^9)/(li)
             #NVTdataobj@norm1 <- 10^9*ci1/li*N1

             ci2 <- NVTdataobj@exp2[,1]
             N2 <- sum(ci2)
             NVTdataobj@norm2 <- ci2/(N1/10^9)/(li)
             #NVTdataobj@norm2 <- 10^9*ci2/li*N2
             colnames(NVTdataobj@norm1)=names(NVTdataobj@exp1)
             colnames(NVTdataobj@norm2)=names(NVTdataobj@exp2)
           },
           RPM={
             print ("RPM normalization!")

             if(nrow(NVTdataobj@length)==0){
               stop("No RPKM normalization possible without gene length")
             }

             li <- NVTdataobj@length

             ci1 <- NVTdataobj@exp1[,1]
             N1 <- sum(ci1)
             NVTdataobj@norm1 <- ci1/(N1/10^9)

             ci2 <- NVTdataobj@exp2[,1]
             N2 <- sum(ci2)
             NVTdataobj@norm1 <- ci2/(N1/10^9)

             colnames(NVTdataobj@norm1)=names(NVTdataobj@exp1)
             colnames(NVTdataobj@norm2)=names(NVTdataobj@exp2)
           },
           RPKM2={
             print ("RPKM normalization!")

             if(nrow(NVTdataobj@length)==0){
               stop("No RPKM normalization possible without gene length")
             }

             library("NOISeq")
             mymatrix <- as.matrix(cbind(NVTdataobj@exp1,NVTdataobj@exp2))

             mynormmatrix <- rpkm(mymatrix,long=NVTdataobj@length[,1])
             NVTdataobj@norm1 <- as.data.frame(mynormmatrix[,1])
             NVTdataobj@norm2 <- as.data.frame(mynormmatrix[,2])
             colnames(NVTdataobj@norm1)=names(NVTdataobj@exp1)
             colnames(NVTdataobj@norm2)=names(NVTdataobj@exp2)
           },
           DEQ={
             print ("DESeq normalization!")

             library("DESeq")
             condition = factor( c( "untreated", "treated"))
             mymatrix <- as.matrix(cbind(NVTdataobj@exp1,NVTdataobj@exp2))
             cds = newCountDataSet( mymatrix, condition )
             cds = estimateSizeFactors( cds )
             NVTdataobj@norm1 <- as.data.frame(counts( cds, normalized=TRUE )[,1])
             NVTdataobj@norm2 <- as.data.frame(counts( cds, normalized=TRUE )[,2])
             colnames(NVTdataobj@norm1)=names(NVTdataobj@exp1)
             colnames(NVTdataobj@norm2)=names(NVTdataobj@exp2)
           },
           TPM={
             print ("TPM normalization!")

             if(nrow(NVTdataobj@length)==0){
               stop("No TPM normalization possible without gene length")
             }

             li <- NVTdataobj@length
             ml <- mean (NVTdataobj@length[,1])

             ci1 <- NVTdataobj@exp1[,1]
             N1 <- sum(ci1)
             #RPKM1 <- 10^9*ci1/li*N1
             #NVTdataobj@norm1 <- ml*RPKM1/sum(RPKM1)*10^3
             RPK1 <- ci1/(li/1000)
             NVTdataobj@norm1 <- RPK1/(sum(RPK1)/10^6)

             ci2 <- NVTdataobj@exp2[,1]
             N2 <- sum(ci2)
             #RPKM2 <- 10^9*ci2/li*N2
             #NVTdataobj@norm2 <- ml*RPKM2/sum(RPKM2)*10^3
             RPK2 <- ci2/(li/1000)
             NVTdataobj@norm2 <- RPK2/(sum(RPK2)/10^6)

             colnames(NVTdataobj@norm1)=names(NVTdataobj@exp1)
             colnames(NVTdataobj@norm2)=names(NVTdataobj@exp2)
           },
           G={
             print ("Normalization be given gene-set!")
              print(NVTdataobj@hklist)
             gn1 <- mean(NVTdataobj@exp1[NVTdataobj@hklist,])
             gn2 <- mean(NVTdataobj@exp2[NVTdataobj@hklist,])
             NVTdataobj@norm1 <- NVTdataobj@exp1/gn1
             NVTdataobj@norm2 <- NVTdataobj@exp2/gn2

           },
           {
             print ("No normalization!")

             NVTdataobj@norm1 <- NVTdataobj@exp1
             NVTdataobj@norm2 <- NVTdataobj@exp2
           }
    )

    return(NVTdataobj)
  }else{
    stop("Not a valid NVTdata object!")
  }
}

#plot data
NVTplot <- function(NVTdataobj) {
  if(check_expression_list(NVTdataobj@exp1) && check_expression_list(NVTdataobj@exp2)
     && check_hkgene_list(NVTdataobj@hklist) && check_method(NVTdataobj@method)
     && check_expression_list(NVTdataobj@norm1)  && check_expression_list(NVTdataobj@norm2)){

    plot(log(NVTdataobj@norm1[,1]),log(NVTdataobj@norm2[,1]),main=paste("MA-plot", names(NVTdataobj@norm1),"vs.",names(NVTdataobj@norm2)),
         xlab=paste("log( normalized expression",names(NVTdataobj@norm1),")"),ylab=paste("log( normalized expression",names(NVTdataobj@norm2),")")
         ,asp=1,pch=20,col="grey")

    m1 <- log(NVTdataobj@norm1[NVTdataobj@hklist,])
    m2 <- log(NVTdataobj@norm2[NVTdataobj@hklist,])

    points(m1,m2,col="blue",pch=19)

    #fm <- lm(m2[,1] ~ m1[,1])
    #abline(fm, col = "red")

  }else{
    stop("Not a valid NVTdata object with normalized values!")
  }
}

#calculate pearson correclation of housekeeping genes
NVTpearson <- function(NVTdataobj) {
  if(check_expression_list(NVTdataobj@exp1) && check_expression_list(NVTdataobj@exp2)
     && check_hkgene_list(NVTdataobj@hklist) && check_method(NVTdataobj@method)
     && check_expression_list(NVTdataobj@norm1)  && check_expression_list(NVTdataobj@norm2)){

    m1 <- log(NVTdataobj@norm1[NVTdataobj@hklist,])
    m2 <- log(NVTdataobj@norm2[NVTdataobj@hklist,])

    pearson <- cor(m1,m2,method="pearson")
    return(pearson)

  }else{
    stop("Not a valid NVTdata object with normalized values!")
  }
}

#checks
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
