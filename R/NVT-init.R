#global methods availible
method_v <- c("N","TC","Med","TMM","UQ","UQ2","Q","RPKM","RPM","DEQ","TPM","G")
#"RPKM2",

#'initialize and load input dato into NVTobject
#'
#'@param hkgene_list A list of housekeeping-genes
#'@param exp_list1 The first data frame of expression values per gene
#'@param exp_list2 Second data frame of expression values per gene
#'@param method The normalization method to us [N,TC,Med,TMM,UQ,UQ2,Q,RPKM,RPM,DEQ,TPM,G]
#'@param length A data frame of length per gene
#'@return A NVTobject ready for normalization
#'@examples
#'myexp1 = read.table("exp1.txt")
#'myexp2 = read.table("exp2.txt")
#'mylen = read.table("length.txt")
#'mylist1=c('gene1','gene2','gene3','gene10')
#'
#'mynvt <- NVTinit(mylist1,myexp1,myexp2,"N")
#'mynvt2 <- NVTinit(mylist1,myexp1,myexp2,"RPKM",mylen)
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
      return(new("NVTdata",exp1=exp_list1,exp2=exp_list2,hklist=hkgene_list,method=method))
    }
  }

  #all here add length
  if(check_expression_list(exp_list1) && check_expression_list(exp_list2) && check_hkgene_list(hkgene_list) && check_method(method)) {
    return(new("NVTdata",exp1=exp_list1,exp2=exp_list2,hklist=hkgene_list,method=method,length=length))
  }
}

#'Normalize a NVTobject,
#'the normalization method has been set already in the initialization step
#'
#'@param NVTdataobj A previously initialized NVTobject
#'@return A normalized NVTobject
#'@examples
#'myexp1 = read.table("data/exp1.txt")
#'myexp2 = read.table("data/exp2.txt")
#'mylen = read.table("data/length.txt")
#'mylist1=c('gene1','gene2','gene3','gene10')
#'
#'mynvt <- NVTinit(mylist1,myexp1,myexp2,"N")
#'
#'mynorm <- NVTnormalize(mynvt)
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
             NVTdataobj@norm1 <-  as.data.frame(ci1/N1*m)
             NVTdataobj@norm2 <-  as.data.frame(ci2/N2*m)
           },
           Med={
             print ("Median normalization!")

             ci1 <- NVTdataobj@exp1[,1]
             N1 <- median(ci1)
             ci2 <- NVTdataobj@exp2[,1]
             N2 <- median(ci2)
             m=mean(c(N1,N2))
             NVTdataobj@norm1 <-  as.data.frame(ci1/N1*m)
             NVTdataobj@norm2 <-  as.data.frame(ci2/N2*m)

           },
           TMM={
             print ("Trimmed Mean of M-values normalization!")

             if (requireNamespace("NOISeq", quietly = TRUE)) {
               #library("NOISeq")
               mymatrix <- as.matrix(cbind(NVTdataobj@exp1,NVTdataobj@exp2))

               mynormmatrix <- tmm(mymatrix)
               NVTdataobj@norm1 <- as.data.frame(mynormmatrix[,1])
               NVTdataobj@norm2 <- as.data.frame(mynormmatrix[,2])
             } else {

               print ("NOISeq package not available, data not normalized!")
               NVTdataobj@norm1 <- NVTdataobj@exp1
               NVTdataobj@norm2 <- NVTdataobj@exp2
             }

           },
           UQ={
             print ("Upper Quartile normalization!")

             ci1 <- NVTdataobj@exp1[,1]
             ci2 <- NVTdataobj@exp2[,1]

             ci <- cbind(ci1,ci2)
             #row_sub <- apply(ci, 1, function(row) all(row !=0 )) #remove all 0 rows
             #res <- ci[row_sub,]

             #N1 <- quantile(res[,1], 0.75)
             #N2 <- quantile(res[,2], 0.75)
             N1 <- quantile(ci1, 0.75)
             N2 <- quantile(ci2, 0.75)
             m=mean(N1,N2)

             NVTdataobj@norm1 <- as.data.frame( ci1/N1*m )
             NVTdataobj@norm2 <- as.data.frame( ci2/N2*m )
           },
           UQ2={
             print ("Upper Quartile normalization (from NOISeq)!")

             if (requireNamespace("NOISeq", quietly = TRUE)) {
               #library("NOISeq")
               mymatrix <- as.matrix(cbind(NVTdataobj@exp1,NVTdataobj@exp2))

               mynormmatrix <- uqua(mymatrix, long = 1000, lc = 0, k = 0)
               NVTdataobj@norm1 <- as.data.frame(mynormmatrix[,1])
               NVTdataobj@norm2 <- as.data.frame(mynormmatrix[,2])
             } else {

               print ("NOISeq package not available, data not normalized!")
               NVTdataobj@norm1 <- NVTdataobj@exp1
               NVTdataobj@norm2 <- NVTdataobj@exp2
             }

           },
           Q={
             print ("Quantile normalization!")

             if (requireNamespace("limma", quietly = TRUE)) {
               #library("limma")

               mymatrix <- as.matrix(cbind(NVTdataobj@exp1,NVTdataobj@exp2))
               mynormmatrix <- normalizeBetweenArrays(mymatrix,"quantile")
               NVTdataobj@norm1 <- as.data.frame(mynormmatrix[,1])
               NVTdataobj@norm2 <- as.data.frame(mynormmatrix[,2])

             } else {
               print ("limma package not available, data not normalized!")
               NVTdataobj@norm1 <- NVTdataobj@exp1
               NVTdataobj@norm2 <- NVTdataobj@exp2
             }

           },
           RPKM={
             print ("RPKM normalization!")

             if(nrow(NVTdataobj@length)==0){
               stop("No RPKM normalization possible without gene length")
             }

             li <- NVTdataobj@length

             ci1 <- NVTdataobj@exp1[,1]
             N1 <- sum(ci1)
             NVTdataobj@norm1 <- as.data.frame(ci1/(N1/10^9)/(li))
             #NVTdataobj@norm1 <- 10^9*ci1/li*N1

             ci2 <- NVTdataobj@exp2[,1]
             N2 <- sum(ci2)
             NVTdataobj@norm2 <- as.data.frame(ci2/(N1/10^9)/(li))
             #NVTdataobj@norm2 <- 10^9*ci2/li*N2
           },
           RPM={
             print ("RPM normalization!")

             if(nrow(NVTdataobj@length)==0){
               stop("No RPKM normalization possible without gene length")
             }

             li <- NVTdataobj@length

             ci1 <- NVTdataobj@exp1[,1]
             N1 <- sum(ci1)
             NVTdataobj@norm1 <- as.data.frame(ci1/(N1/10^9))

             ci2 <- NVTdataobj@exp2[,1]
             N2 <- sum(ci2)
             NVTdataobj@norm2 <- as.data.frame(ci2/(N1/10^9))
           },
           #RPKM2={
           #  print ("RPKM normalization!")

           #  if(nrow(NVTdataobj@length)==0){
           #    stop("No RPKM normalization possible without gene length")
           #  }

           #  library("NOISeq")
           #  mymatrix <- as.matrix(cbind(NVTdataobj@exp1,NVTdataobj@exp2))

           #  mynormmatrix <- rpkm(mymatrix,long=NVTdataobj@length[,1])
           #  NVTdataobj@norm1 <- as.data.frame(mynormmatrix[,1])
           #  NVTdataobj@norm2 <- as.data.frame(mynormmatrix[,2])
           #},
           DEQ={
             print ("DESeq normalization!")

             if (requireNamespace("DESeq", quietly = TRUE)) {
               #library("DESeq")
               condition = factor( c( "untreated", "treated"))

               if(is.integer(NVTdataobj@exp1) && is.integer(NVTdataobj@exp1) ){
                 mymatrix <- as.matrix(cbind(NVTdataobj@exp1,NVTdataobj@exp2))
               }else{
                 print ("Input counts are not integer, converting them!")
                 mymatrix <- as.matrix(cbind(as.integer(NVTdataobj@exp1[,1]),as.integer(NVTdataobj@exp2[,1])))
               }

               cds = newCountDataSet( mymatrix, condition )
               cds = estimateSizeFactors( cds )
               NVTdataobj@norm1 <- as.data.frame(counts( cds, normalized=TRUE )[,1])
               NVTdataobj@norm2 <- as.data.frame(counts( cds, normalized=TRUE )[,2])

             } else {

               print ("DESeq package not available, data not normalized!")
               NVTdataobj@norm1 <- NVTdataobj@exp1
               NVTdataobj@norm2 <- NVTdataobj@exp2
             }

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
             NVTdataobj@norm1 <- as.data.frame(RPK1/(sum(RPK1)/10^6))

             ci2 <- NVTdataobj@exp2[,1]
             N2 <- sum(ci2)
             #RPKM2 <- 10^9*ci2/li*N2
             #NVTdataobj@norm2 <- ml*RPKM2/sum(RPKM2)*10^3
             RPK2 <- ci2/(li/1000)
             NVTdataobj@norm2 <- as.data.frame(RPK2/(sum(RPK2)/10^6))
           },
           G={
             print ("Normalization be given gene-set!")
              print(NVTdataobj@hklist)
             gn1 <- mean(NVTdataobj@exp1[NVTdataobj@hklist,])
             gn2 <- mean(NVTdataobj@exp2[NVTdataobj@hklist,])
             NVTdataobj@norm1 <- as.data.frame(NVTdataobj@exp1/gn1)
             NVTdataobj@norm2 <- as.data.frame(NVTdataobj@exp2/gn2)

           },
           {
             print ("No normalization!")

             NVTdataobj@norm1 <- NVTdataobj@exp1
             NVTdataobj@norm2 <- NVTdataobj@exp2
           }
    )

    #set names
    names(NVTdataobj@norm1)=names(NVTdataobj@exp1)
    names(NVTdataobj@norm2)=names(NVTdataobj@exp2)
    rownames(NVTdataobj@norm1)=rownames(NVTdataobj@exp1)
    rownames(NVTdataobj@norm2)=rownames(NVTdataobj@exp2)

    return(NVTdataobj)
  }else{
    stop("Not a valid NVTdata object!")
  }
}

#'Plot the MA-plot of a NVTobject
#'
#'@param NVTdataobj A previously initialized and normalized NVTobject
#'@return Plots the MA-plot with the housekeeping genes indicated
#'@examples
#'myexp1 = read.table("data/exp1.txt")
#'myexp2 = read.table("data/exp2.txt")
#'mylen = read.table("data/length.txt")
#'mylist1=c('gene1','gene2','gene3','gene10')
#'
#'mynvt <- NVTinit(mylist1,myexp1,myexp2,"N")
#'mynorm <- NVTnormalize(mynvt)
#'
#'NVTplot(mynorm)
NVTplot <- function(NVTdataobj) {
  if(check_expression_list(NVTdataobj@exp1) && check_expression_list(NVTdataobj@exp2)
     && check_hkgene_list(NVTdataobj@hklist) && check_method(NVTdataobj@method)
     && check_norm_list(NVTdataobj@norm1)  && check_norm_list(NVTdataobj@norm2)){

    l1 <- log(NVTdataobj@norm1[[1]])
    l2 <- log(NVTdataobj@norm2[[1]])
    names(l1) <- rownames(NVTdataobj@norm1)
    names(l2) <- rownames(NVTdataobj@norm2)

    #clean up
    l <- cbind(l1,l2)
    idx <- apply(l, 1, function(x) all(!is.infinite(x)))
    l <- l[idx,]
    idx <- apply(l, 1, function(x) all(!is.na(x)))
    l <- l[idx,]
    idx <- apply(l, 1, function(x) all(!is.nan(x)))
    l <- l[idx,]

    min <- min(l)
    max <- max(l)

    #only houskeeping genes
    m1 <- l1[NVTdataobj@hklist]
    m2 <- l2[NVTdataobj@hklist]

    #clean up for lm
    m <- cbind(m1,m2)
    idx <- apply(m, 1, function(x) all(!is.infinite(x)))
    m <- m[idx,]
    idx <- apply(m, 1, function(x) all(!is.na(x)))
    m <- m[idx,]
    idx <- apply(m, 1, function(x) all(!is.nan(x)))
    m <- m[idx,]

    plot(l1,l2,main=paste("MA-plot", names(NVTdataobj@norm1),"vs.",names(NVTdataobj@norm2)),
         xlab=paste("log( normalized expression",names(NVTdataobj@norm1),")"),
         ylab=paste("log( normalized expression",names(NVTdataobj@norm2),")")
         ,pch=20,col=rgb(193,205,205,90,maxColorValue=255),xlim=c(min, max),ylim=c(min, max))
    mtext(paste(NVTdataobj@method,"normalized"))

    points(m1,m2,col="blue",pch=19)


    fm <- lm(m[,2] ~ m[,1])

    #print lines
    abline(0, 1, col = "black", lwd = 1, lty = 2)
    abline(fm, col = "red")

  }else{
    stop("Not a valid NVTdata object with normalized values!")
  }
}

#'Calculate the pearson correclation of the housekeeping genes of an initialized and normalized NVTobject
#'
#'@param NVTdataobj A previously initialized and normalized NVTobject
#'@return Pearson correlation of the normalized housekeeping genes between the samples
#'@examples
#'myexp1 = read.table("data/exp1.txt")
#'myexp2 = read.table("data/exp2.txt")
#'mylen = read.table("data/length.txt")
#'mylist1=c('gene1','gene2','gene3','gene10')
#'
#'mynvt <- NVTinit(mylist1,myexp1,myexp2,"N")
#'mynorm <- NVTnormalize(mynvt)
#'
#'NVTpearson(mynorm)
NVTpearson <- function(NVTdataobj) {
  if(check_expression_list(NVTdataobj@exp1) && check_expression_list(NVTdataobj@exp2)
     && check_hkgene_list(NVTdataobj@hklist) && check_method(NVTdataobj@method)
     && check_norm_list(NVTdataobj@norm1)  && check_norm_list(NVTdataobj@norm2)){

    m1 <- NVTdataobj@norm1[1][NVTdataobj@hklist,]
    m2 <- NVTdataobj@norm2[1][NVTdataobj@hklist,]

    pearson <- cor(m1,m2,method="pearson")
    return(pearson)

  }else{
    stop("Not a valid NVTdata object with normalized values!")
  }
}


#'Calculate the pearson correclation of the housekeeping genes of an initialized NVTobject for all normalization methods
#'
#'@param NVTdataobj A previously initialzed and normalized NVTobject
#'@return Sorted pearson correlations of the normalized housekeeping genes between the samples
#'@examples
#'myexp1 = read.table("data/exp1.txt")
#'myexp2 = read.table("data/exp2.txt")
#'mylen = read.table("data/length.txt")
#'mylist1=c('gene1','gene2','gene3','gene10')
#'
#'mynvt <- NVTinit(mylist1,myexp1,myexp2,"N")
#'
#'NVTtestall(mynvt)
NVTtestall <- function(NVTdataobj) {
  tmpNVT <- NVTdataobj
  pv <- vector()

  #test all methods and extract pearson correlation
  for (n in method_v) {
    tmpNVT@norm_method <- n
    tmpnorm <- NVTnormalize(tmpNVT)
    p <- NVTpearson(tmpnorm)
    pv <-append(pv,p)
  }
  pearson <- as.data.frame(cbind(method_v,pv))

  return(pearson[order(pearson$pv,decreasing = T),])
}


#'Check housekeeping gene list
#'
#'@param hkgene_list Housekeeping gene list
#'@return true or false
check_hkgene_list <- function(hkgene_list) {
  if( is.character(hkgene_list) ){
    return(TRUE)
  }else{
    stop("Housekeeping-gene list is no list!")
  }
}

#'Check expression gene list
#'
#'@param exp_list Expression gene list
#'@return true or false
check_expression_list <- function(exp_list) {
  if( is.data.frame(exp_list) || is.list(exp_list) || is.vector(exp_list)){
    return(TRUE)
  }else{
    stop("Expression list is no data frame, list or vector!")
  }
}

#'Check normalized expression gene list
#'
#'@param norm_list Normalized expression gene list
#'@return true or false
check_norm_list <- function(norm_list) {
  if( is.data.frame(norm_list) || is.list(norm_list) || is.vector(norm_list)){
    if(length(norm_list)!=0){
      return(TRUE)
    }else{
      stop("Normalized-list is empty, use NVTnormalize() first!")
    }

  }else{
    stop("Expression list is no data frame, list or vector!")
  }
}

#'Check normalization method
#'
#'@param method Normalization method
#'@return true or false
check_method <- function(method) {
  if( method %in% method_v ){
    return(TRUE)
  }else{
    stop("Unknown method specified!")
  }
}
