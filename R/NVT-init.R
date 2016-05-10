#global methods availible
method_v <- c("N","TC","Med","TMM","UQ","UQ2","Q","RPKM","RPM","DEQ","TPM","G")
method_c <- c("p","rmsd","mae")

#'Initialize and load input data into a NVTobject
#'
#'@export
#'@param housekeeping_gene_list A list of housekeeping-genes
#'@param expression_list_1 The first data frame of expression values per gene
#'@param expression_list_2 Second data frame of expression values per gene
#'@param normalization_method The normalization method to use [N,TC,Med,TMM,UQ,UQ2,Q,RPKM,RPM,DEQ,TPM,G] N = No normalization, TC = Total count normalization, Med = Median normalization, TMM = Trimmed Mean of M-values normalization, UQ = Upper Quartile normalization , UQ2 = Upper Quartile normalization (from NOISeq), Q = Quantile normalization, RPKM = Reads Per Kilobase per Million mapped reads normalization, RPM = Reads per Million mapped reads normalization, DEQ = relative log expression method included in the DESeq package, TPM = transcripts per million normalization, G = use the provided genes to normalize
#'@param feat_length A data frame of length per gene/exon
#'@return An NVTobject ready for normalization
#'@examples
#'library("NVT")
#'data(myexp1)
#'data(myexp2)
#'data(mylen)
#'mylist1<-c("ENSG00000111640","ENSG00000163631","ENSG00000075624","ENSG00000172053",
#'"ENSG00000170950","ENSG00000165704","ENSG00000196839","ENSG00000168938","ENSG00000177700")
#'
#'mynvt <- NVTinit(mylist1,myexp1,myexp2,"N")
#'mynvt2 <- NVTinit(mylist1,myexp1,myexp2,"RPKM",mylen)
NVTinit <- function(housekeeping_gene_list, expression_list_1, expression_list_2, normalization_method, feat_length){

  if (missing('housekeeping_gene_list')){
    stop("No housekeeping-gene list specified!")
  }
  if (missing('expression_list_1')){
    stop("No expression list specified!")
  }
  if (missing('expression_list_2')){
    stop("No second expression list specified!")
  }
  if (missing('normalization_method')){
    stop("No method specified!")
  }

  if( all(rownames(expression_list_1) %in% rownames(expression_list_2)) && all(rownames(expression_list_2) %in% rownames(expression_list_1))  ){

    #check list of hk genes length
    if(length(housekeeping_gene_list) == 1){
      print("Warning: only one element in the housekeeping-gene-list this will disable plotting and the pearson correlation calculation")
    }

   #length missing
   if (missing('feat_length')){

      if(check_expression_list(expression_list_1) && check_expression_list(expression_list_2) && check_hkgene_list(housekeeping_gene_list) && check_method(normalization_method)) {

       #sort expression data via rownames
       expression_list_1 <- expression_list_1[order(rownames(expression_list_1)), , drop = FALSE]
       expression_list_2 <- expression_list_2[order(rownames(expression_list_2)), , drop = FALSE]

       return(new("NVTdata",exp1=expression_list_1,exp2=expression_list_2,hklist=housekeeping_gene_list,norm_method=normalization_method))
     }
   }

   #all here add length
   if(check_expression_list(expression_list_1) && check_expression_list(expression_list_2) && check_hkgene_list(housekeeping_gene_list) && check_method(normalization_method) ) {

     #check length data
     if( all(rownames(expression_list_1) %in% rownames(feat_length)) && all(rownames(expression_list_2) %in% rownames(feat_length)) ){

        #make sure its correctly sorted
        feat_length <- subset(feat_length, rownames(feat_length) %in% rownames(expression_list_1))
        feat_length <- feat_length[order(rownames(feat_length)), , drop = FALSE]
        expression_list_1 <- expression_list_1[order(rownames(expression_list_1)), , drop = FALSE]
        expression_list_2 <- expression_list_2[order(rownames(expression_list_2)), , drop = FALSE]

        return(new("NVTdata",exp1=expression_list_1,exp2=expression_list_2,hklist=housekeeping_gene_list,norm_method=normalization_method,length=feat_length))
     }else{
       stop("Not all expressed genes have a length in the provided data-set!")
     }
   }
  }else{
    stop("The gene/exon names in the expression data sets are unequal/different!")
  }
}

#'Normalize a NVTobject, the normalization method has been set already in the initialization step
#'
#'@export
#'@param NVTdataobj A previously initialized NVTobject
#'@param verbose mode on or off (T/F) [default=TRUE]
#'@return A normalized NVTobject
#'@examples
#'library("NVT")
#'
#'data(myexp1)
#'data(myexp2)
#'data(mylen)
#'mylist1<-c("ENSG00000111640","ENSG00000163631","ENSG00000075624","ENSG00000172053",
#'"ENSG00000170950","ENSG00000165704","ENSG00000196839","ENSG00000168938","ENSG00000177700")
#'
#'
#'mynvt <- NVTinit(mylist1,myexp1,myexp2,"N",mylen)
#'
#'mynorm <- NVTnormalize(mynvt)
NVTnormalize <- function(NVTdataobj, verbose=TRUE) {

  if(check_expression_list(NVTdataobj@exp1) && check_expression_list(NVTdataobj@exp2)
     && check_hkgene_list(NVTdataobj@hklist) && check_method(NVTdataobj@norm_method)
     && exists_hkgene_list(NVTdataobj,NVTdataobj@hklist)){
    switch(NVTdataobj@norm_method,
           N={
             if(verbose){print ("No normalization!")}
             NVTdataobj@norm_method_name="Not"

             NVTdataobj@norm1 <- NVTdataobj@exp1
             NVTdataobj@norm2 <- NVTdataobj@exp2
             NVTdataobj@is_norm = TRUE
           },
           TC={
             if(verbose){print ("Total count normalization!")}
             NVTdataobj@norm_method_name="Total count"

             ci1 <- NVTdataobj@exp1[,1]
             N1 <- sum(ci1)
             ci2 <- NVTdataobj@exp2[,1]
             N2 <- sum(ci2)
             m=mean(c(N1,N2))
             NVTdataobj@norm1 <-  as.data.frame(ci1/N1*m)
             NVTdataobj@norm2 <-  as.data.frame(ci2/N2*m)
             NVTdataobj@is_norm = TRUE
           },
           Med={
             if(verbose){print ("Median normalization!")}
             NVTdataobj@norm_method_name="Median"

             ci1 <- NVTdataobj@exp1[,1]
             N1 <- median(ci1[ci1>0])#non 0 values
             ci2 <- NVTdataobj@exp2[,1]
             N2 <- median(ci2[ci2>0])#non 0 values
             m=mean(c(N1,N2))
             NVTdataobj@norm1 <-  as.data.frame(ci1/N1*m)
             NVTdataobj@norm2 <-  as.data.frame(ci2/N2*m)
             NVTdataobj@is_norm = TRUE
           },
           TMM={
             if(verbose){print ("Trimmed Mean of M-values normalization!")}
             NVTdataobj@norm_method_name="Trimmed Mean of M-values (TMM)"

             if (requireNamespace("NOISeq", quietly = TRUE)) {
               #library("NOISeq")
               mymatrix <- as.matrix(cbind(NVTdataobj@exp1,NVTdataobj@exp2))

               mynormmatrix <- NOISeq::tmm(mymatrix)
               NVTdataobj@norm1 <- as.data.frame(mynormmatrix[,1])
               NVTdataobj@norm2 <- as.data.frame(mynormmatrix[,2])
               NVTdataobj@is_norm = TRUE
             } else {

               print ("NOISeq package not available, data not normalized!")
               NVTdataobj@norm1 <- NVTdataobj@exp1
               NVTdataobj@norm2 <- NVTdataobj@exp2
               NVTdataobj@is_norm = FALSE
             }

           },
           UQ={
             if(verbose){print ("Upper Quartile normalization!")}
             NVTdataobj@norm_method_name="Upper Quartile"

             ci1 <- NVTdataobj@exp1[,1]
             ci2 <- NVTdataobj@exp2[,1]

             ci <- cbind(ci1,ci2)
             #row_sub <- apply(ci, 1, function(row) all(row !=0 )) #remove all 0 rows
             #res <- ci[row_sub,]

             #N1 <- quantile(res[,1], 0.75)
             #N2 <- quantile(res[,2], 0.75)
             N1 <- quantile(ci1[ci1>0], 0.75) #non 0 values
             N2 <- quantile(ci2[ci2>0], 0.75) #non 0 values
             m=mean(N1,N2)

             NVTdataobj@norm1 <- as.data.frame( ci1/N1*m )
             NVTdataobj@norm2 <- as.data.frame( ci2/N2*m )
             NVTdataobj@is_norm = TRUE
           },
           UQ2={
             if(verbose){print ("Upper Quartile normalization (from NOISeq)!")}
             NVTdataobj@norm_method_name="Upper Quartile (from NOISeq)"

             if (requireNamespace("NOISeq", quietly = TRUE)) {
               #library("NOISeq")
               mymatrix <- as.matrix(cbind(NVTdataobj@exp1,NVTdataobj@exp2))

               mynormmatrix <- NOISeq::uqua(mymatrix, long = 1000, lc = 0, k = 0)
               NVTdataobj@norm1 <- as.data.frame(mynormmatrix[,1])
               NVTdataobj@norm2 <- as.data.frame(mynormmatrix[,2])
               NVTdataobj@is_norm = TRUE
             } else {

               print ("NOISeq package not available, data not normalized!")
               NVTdataobj@norm1 <- NVTdataobj@exp1
               NVTdataobj@norm2 <- NVTdataobj@exp2
               NVTdataobj@is_norm = FALSE
             }

           },
           Q={
             if(verbose){print ("Quantile normalization!")}
             NVTdataobj@norm_method_name="Quantile"

             if (requireNamespace("limma", quietly = TRUE)) {
               #library("limma")

               mymatrix <- as.matrix(cbind(NVTdataobj@exp1,NVTdataobj@exp2))
               mynormmatrix <- limma::normalizeBetweenArrays(mymatrix,"quantile")
               NVTdataobj@norm1 <- as.data.frame(mynormmatrix[,1])
               NVTdataobj@norm2 <- as.data.frame(mynormmatrix[,2])
               NVTdataobj@is_norm = TRUE

             } else {
               print ("limma package not available, data not normalized!")
               NVTdataobj@norm1 <- NVTdataobj@exp1
               NVTdataobj@norm2 <- NVTdataobj@exp2
               NVTdataobj@is_norm = FALSE
             }

           },
           RPKM={
             if(verbose){print ("RPKM normalization!")}

             if(nrow(NVTdataobj@length)==0){
               print ("No RPKM normalization possible without gene length, data not normalized!")
               NVTdataobj@norm1 <- NVTdataobj@exp1
               NVTdataobj@norm2 <- NVTdataobj@exp2
               NVTdataobj@is_norm = FALSE
             }else{

             NVTdataobj@norm_method_name="Reads Per Kilobase per Million mapped reads (RPKM)"
             li <- NVTdataobj@length

             ci1 <- NVTdataobj@exp1[,1]
             N1 <- sum(ci1)
             NVTdataobj@norm1 <- as.data.frame(ci1/(N1/10^9)/(li))
             #NVTdataobj@norm1 <- 10^9*ci1/li*N1

             ci2 <- NVTdataobj@exp2[,1]
             N2 <- sum(ci2)
             NVTdataobj@norm2 <- as.data.frame(ci2/(N1/10^9)/(li))
             #NVTdataobj@norm2 <- 10^9*ci2/li*N2
             NVTdataobj@is_norm = TRUE
            }
           },
           RPM={
             if(verbose){print ("RPM normalization!")}

             NVTdataobj@norm_method_name="Reads Per Million mapped reads (RPM)"

             ci1 <- NVTdataobj@exp1[,1]
             N1 <- sum(ci1)
             NVTdataobj@norm1 <- as.data.frame(ci1/(N1/10^9))

             ci2 <- NVTdataobj@exp2[,1]
             N2 <- sum(ci2)
             NVTdataobj@norm2 <- as.data.frame(ci2/(N1/10^9))
             NVTdataobj@is_norm = TRUE
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
           #  NVTdataobj@is_norm = TRUE
           #},
           DEQ={
             if(verbose){print ("DESeq normalization!")}

             if (requireNamespace("DESeq", quietly = TRUE)) {
               NVTdataobj@norm_method_name="Relative log expression method (DESeq)"
               if(verbose){print ("Using DESeq")}
               condition = factor( c( "untreated", "treated"))

               if(is.integer(NVTdataobj@exp1) && is.integer(NVTdataobj@exp1) ){
                 mymatrix <- as.matrix(cbind(NVTdataobj@exp1,NVTdataobj@exp2))
               }else{
                 if(verbose){print ("Input counts are not integer, converting them!")}
                 mymatrix <- as.matrix(cbind(as.integer(NVTdataobj@exp1[,1]),as.integer(NVTdataobj@exp2[,1])))
               }

               cds <- DESeq::newCountDataSet( mymatrix, condition )
               cds <- DESeq::estimateSizeFactors( cds )
               NVTdataobj@norm1 <- as.data.frame(DESeq::counts( cds, normalized=TRUE )[,1])
               NVTdataobj@norm2 <- as.data.frame(DESeq::counts( cds, normalized=TRUE )[,2])
               NVTdataobj@is_norm = TRUE

             } else if(requireNamespace("DESeq2", quietly = TRUE)){
                NVTdataobj@norm_method_name="Relative log expression method (DESeq)"
                if(verbose){print ("Using DESeq2")}
                condition = factor( c( "untreated", "treated"))

                if(is.integer(NVTdataobj@exp1) && is.integer(NVTdataobj@exp1) ){
                  mymatrix <- as.matrix(cbind(NVTdataobj@exp1,NVTdataobj@exp2))
                }else{
                  if(verbose){print ("Input counts are not integer, converting them!")}
                   mymatrix <- as.matrix(cbind(as.integer(NVTdataobj@exp1[,1]),as.integer(NVTdataobj@exp2[,1])))
                }

                dds <- DESeq2::DESeqDataSetFromMatrix( countData = mymatrix, colData = S4Vectors::DataFrame(condition),design = ~ condition)
                dds <- DESeq2::estimateSizeFactors(dds)

                NVTdataobj@norm1 <- as.data.frame(DESeq2::counts( dds, normalized=TRUE )[,1])
                NVTdataobj@norm2 <- as.data.frame(DESeq2::counts( dds, normalized=TRUE )[,2])
                NVTdataobj@is_norm = TRUE

             } else{
               print ("DESeq and DESeq2 packages not available, data not normalized!")
               NVTdataobj@norm1 <- NVTdataobj@exp1
               NVTdataobj@norm2 <- NVTdataobj@exp2
               NVTdataobj@is_norm = FALSE
             }

           },
           TPM={
             if(verbose){print ("TPM normalization!")}

             if(nrow(NVTdataobj@length)==0){
               print ("No TPM normalization possible without gene length, data not normalized!")
               NVTdataobj@norm1 <- NVTdataobj@exp1
               NVTdataobj@norm2 <- NVTdataobj@exp2
               NVTdataobj@is_norm = FALSE
             }else{

             NVTdataobj@norm_method_name="Transcripts Per Million mapped reads (TPM)"
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
             NVTdataobj@is_norm = TRUE
            }
           },
           G={
             if(verbose){
               print ("Normalization by given gene-set!")
               print (NVTdataobj@hklist)
             }
             NVTdataobj@norm_method_name="Gene-set"

             gn1 <- mean(NVTdataobj@exp1[NVTdataobj@hklist,])
             gn2 <- mean(NVTdataobj@exp2[NVTdataobj@hklist,])
             NVTdataobj@norm1 <- as.data.frame(NVTdataobj@exp1/gn1)
             NVTdataobj@norm2 <- as.data.frame(NVTdataobj@exp2/gn2)
             NVTdataobj@is_norm = TRUE
           },
           {
             if(verbose){print ("No normalization!")}
             NVTdataobj@norm_method_name="Not"

             NVTdataobj@norm1 <- NVTdataobj@exp1
             NVTdataobj@norm2 <- NVTdataobj@exp2
             NVTdataobj@is_norm = TRUE
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

#'Plot the XY-plot of a NVTobject
#'
#'@export
#'@param NVTdataobj A previously initialized and normalized NVTobject
#'@param cex Scaling of points and text relative to the default [default=1]
#'@param ... Arguments passed on to the plot function
#'@return Plots the Scatter-plot with the housekeeping genes indicated
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
#'NVTplot(mynorm,0.8)
NVTplot <- function(NVTdataobj, cex=1,  ...) {

  if(length(NVTdataobj@hklist) == 1){
    stop("Only one element in the housekeeping-gene-list, can not calculate linear model")
  }

  if(check_expression_list(NVTdataobj@exp1) && check_expression_list(NVTdataobj@exp2)
     && check_hkgene_list(NVTdataobj@hklist) && check_method(NVTdataobj@norm_method)
     && check_norm_list(NVTdataobj@norm1)  && check_norm_list(NVTdataobj@norm2)
     && exists_hkgene_list(NVTdataobj,NVTdataobj@hklist)){

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

    plot(l1,l2,main=paste("Scatter plot of:", names(NVTdataobj@norm1),"vs.",names(NVTdataobj@norm2)),
         xlab=paste("log( normalized expression",names(NVTdataobj@norm1),")"),
         ylab=paste("log( normalized expression",names(NVTdataobj@norm2),")")
         ,pch=20,col=rgb(193,205,205,50,maxColorValue=255),xlim=c(min, max),ylim=c(min, max),cex.lab = cex,cex.main = cex, cex.axis = cex,cex = cex, ...)
    mtext(paste(NVTdataobj@norm_method_name,"normalized"),cex=cex)

    fm <- lm(m[,2] ~ m[,1])
    #print lines
    abline(0, 1, lwd = 1, lty = 2, col=rgb(0,0,0,150,maxColorValue=255))
    abline(fm, col = "red")

    #add hk genes
    points(m1,m2,col="blue",pch=20,cex=cex)
    #hk gene names
    text(m1,m2,NVTdataobj@hklist, pos=4, col = "blue",cex=cex)
    #text(m1,m2,NVTdataobj@hklist,cex=0.6, pos=4, col = "blue")

  }else{
    stop("Not a valid NVTdata object with normalized values!")
  }
}


#'Plot the MA-plot of a NVTobject
#'
#'@export
#'@param NVTdataobj A previously initialized and normalized NVTobject
#'@param cex Scaling of points and text relative to the default [default=1]
#'@param ... Arguments passed on to the plot function
#'@return Plots the MA-plot with the housekeeping genes indicated
#'@examples
#'library("NVT")
#'data(myexp1)
#'data(myexp2)
#'data(mylen)
#'mylist1<-c("ENSG00000111640","ENSG00000163631","ENSG00000075624","ENSG00000172053",
#'"ENSG00000165704","ENSG00000196839","ENSG00000168938","ENSG00000177700")
#'
#'mynvt <- NVTinit(mylist1,myexp1,myexp2,"N")
#'mynorm <- NVTnormalize(mynvt)
#'
#'NVTmaplot(mynorm,0.8)
NVTmaplot <- function(NVTdataobj, cex=1,  ...) {

  if(length(NVTdataobj@hklist) == 1){
    stop("Only one element in the housekeeping-gene-list, can not calculate linear model")
  }

  if(check_expression_list(NVTdataobj@exp1) && check_expression_list(NVTdataobj@exp2)
     && check_hkgene_list(NVTdataobj@hklist) && check_method(NVTdataobj@norm_method)
     && check_norm_list(NVTdataobj@norm1)  && check_norm_list(NVTdataobj@norm2)
     && exists_hkgene_list(NVTdataobj,NVTdataobj@hklist)){

    l1 <- NVTdataobj@norm1[[1]]
    l2 <- NVTdataobj@norm2[[1]]
    #mean
    l <- cbind(l1,l2)
    mean <- log2(rowMeans(l))

    #log2 fold change
    fc <- log2(l1/l2)

    names(mean) <- rownames(NVTdataobj@norm1)
    names(fc) <- rownames(NVTdataobj@norm1)
    l <- cbind(fc,mean)

    #clean up
    idx <- apply(l, 1, function(x) all(!is.infinite(x)))
    l <- l[idx,]
    idx <- apply(l, 1, function(x) all(!is.na(x)))
    l <- l[idx,]
    idx <- apply(l, 1, function(x) all(!is.nan(x)))
    l <- l[idx,]

    #only houskeeping genes
    m <- cbind(fc,mean)[NVTdataobj@hklist,]

    #clean up for lm
    idx <- apply(m, 1, function(x) all(!is.infinite(x)))
    m <- m[idx,]
    idx <- apply(m, 1, function(x) all(!is.na(x)))
    m <- m[idx,]
    idx <- apply(m, 1, function(x) all(!is.nan(x)))
    m <- m[idx,]

    #,xlim=c(min, max),ylim=c(min, max)
    plot(l[,2],l[,1],main=paste("MA-plot of:", names(NVTdataobj@norm1),"vs.",names(NVTdataobj@norm2)),
         ylab=paste("log2( normalized expression",names(NVTdataobj@norm1),"/",names(NVTdataobj@norm2),")"),
         xlab=paste("log2(mean normalized expression)")
         ,pch=20,col=rgb(193,205,205,50,maxColorValue=255),cex.lab = cex,cex.main = cex, cex.axis = cex,cex = cex, ...)
    mtext(paste(NVTdataobj@norm_method_name,"normalized"),cex=cex)

    fm <- lm(m[,1] ~ m[,2])
    #print lines
    abline(fm, col = "red")

    #add hk genes
    points(m[,2],m[,1],col="blue",pch=20,cex=cex)
    #hk gene names
    text(m[,2],m[,1],NVTdataobj@hklist, pos=4, srt = 320, col = "blue",cex=cex)

  }else{
    stop("Not a valid NVTdata object with normalized values!")
  }
}

#'Plot the XY-plot of a NVTobject with ggplot2
#'
#'@export
#'@param NVTdataobj A previously initialized and normalized NVTobject
#'@param p_cex Point size factor [default = 1]
#'@param t_cex Title size factor [default = 1]
#'@param l_cex Label size factor [default = 1]
#'@return Plots the Scatter-plot with the housekeeping genes indicated
#'@examples
#'library("NVT")
#'data(myexp1)
#'data(myexp2)
#'data(mylen)
#'mylist1<-c("ENSG00000111640","ENSG00000163631","ENSG00000075624","ENSG00000172053",
#'"ENSG00000165704","ENSG00000196839","ENSG00000168938","ENSG00000177700")
#'
#'mynvt <- NVTinit(mylist1,myexp1,myexp2,"N")
#'mynorm <- NVTnormalize(mynvt)
#'
#'NVTadvancedplot(mynorm,2,2,2)
NVTadvancedplot <- function(NVTdataobj,p_cex=1,t_cex=1,l_cex=1) {

  if(length(NVTdataobj@hklist) == 1){
    stop("Only one element in the housekeeping-gene-list, can not calculate linear model")
  }

  if (requireNamespace("ggplot2", quietly = TRUE)  ) {
  #only needed for dsensity plots
  #&& requireNamespace("gridExtra", quietly = TRUE)

    if(check_expression_list(NVTdataobj@exp1) && check_expression_list(NVTdataobj@exp2)
     && check_hkgene_list(NVTdataobj@hklist) && check_method(NVTdataobj@norm_method)
     && check_norm_list(NVTdataobj@norm1)  && check_norm_list(NVTdataobj@norm2)
     && exists_hkgene_list(NVTdataobj,NVTdataobj@hklist)){

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

     myl <- as.data.frame(l)

     mysubl <-  myl[NVTdataobj@hklist,]
     myl[,3] <- rownames(myl)
     colnames(myl)[3] <- "names"

     #empty plot as spacing of the density plots
     #empty <- ggplot()+geom_point(aes(1,1), colour="white") +  theme(plot.background = element_blank(), panel.grid.major = element_blank(),  panel.grid.minor = element_blank(), panel.border = element_blank(),  panel.background = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank())

     d <-  ggplot2::ggplot(data=myl, ggplot2::aes(x = l1, y = l2)) + ggplot2::geom_point( alpha = 0.2,  color="gray", cex=p_cex)
     d <- d + ggplot2::geom_abline(intercept = 0, slope = 1, alpha = 0.5, linetype=2, color="black")
     d <- d + ggplot2::geom_smooth(data=subset(myl,names %in% NVTdataobj@hklist), method=lm,fullrange=TRUE, se=FALSE,alpha = 0.1, color="red",lwd = 0.5)
     d <- d + ggplot2::geom_point(data=subset(myl,names %in% NVTdataobj@hklist), ggplot2::aes(x = l1, y = l2), alpha = 0.8,  color="blue", cex=p_cex)

     #tested labeling with directlabels
     #if (requireNamespace("directlabels", quietly = TRUE)  ) {
     # d <- d + directlabels::geom_dl(mapping=ggplot2::aes(label=names), data=subset(myl,names %in% NVTdataobj@hklist), method="top.bumptwice",inherit.aes = TRUE, color="blue", size=2)
     #}else{
      # print("directlabels not found, printing normal labels")

     d <- d + ggplot2::geom_text(data=subset(myl,names %in% NVTdataobj@hklist), ggplot2::aes(label=names), hjust=-0.2, vjust=0.5, color="blue", cex=l_cex)
     d <- d + ggplot2::ggtitle(bquote(atop(.(paste("Scatter plot of:", names(NVTdataobj@norm1),"vs",names(NVTdataobj@norm2))), atop(italic(.(paste(NVTdataobj@norm_method_name,"normalized"))), ""))))
     d <- d + ggplot2::xlab(paste("log( normalized expression",names(NVTdataobj@norm1),")"))
     d <- d + ggplot2::ylab(paste("log( normalized expression",names(NVTdataobj@norm2),")"))
     d <- d + ggplot2::geom_rug(col="darkred",alpha=.1,position='jitter')
     d <- d + ggplot2::theme(axis.title.x = ggplot2::element_text(size = 10 * t_cex), axis.title.y = ggplot2::element_text(size = 10 * t_cex), plot.title = ggplot2::element_text(size = 10 * t_cex), axis.text.x = ggplot2::element_text(size = 10 * t_cex), axis.text.y = ggplot2::element_text(size = 10 * t_cex) )
     d

     #desnity plots
     #plot_top <- ggplot(data=myl, aes(l1)) +
     # geom_density(alpha=.5) +
     # scale_fill_manual(values = c("orange", "purple")) +
     # theme(legend.position = "none")
     #plot_right <- ggplot(data=myl, aes(l2)) +
     #  geom_density(alpha=.5) +
     #  coord_flip() +
     # scale_fill_manual(values = c("orange", "purple")) +
     #  theme(legend.position = "none")
     #grid.arrange(plot_top, empty, d, plot_right, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))

   }else{
     stop("Not a valid NVTdata object with normalized values!")
   }
  }else{
    print("ggplot2 not found, please use the NVTplot function!")
  }
}

#'Plot the MA-plot of a NVTobject with ggplot2
#'
#'@export
#'@param NVTdataobj A previously initialized and normalized NVTobject
#'@param p_cex Point size factor [default = 1]
#'@param t_cex Title size factor [default = 1]
#'@param l_cex Label size factor [default = 1]
#'@return Plots the MA-plot with the housekeeping genes indicated
#'@examples
#'library("NVT")
#'data(myexp1)
#'data(myexp2)
#'data(mylen)
#'mylist1<-c("ENSG00000111640","ENSG00000163631","ENSG00000075624","ENSG00000172053",
#'"ENSG00000165704","ENSG00000196839","ENSG00000168938","ENSG00000177700")
#'
#'mynvt <- NVTinit(mylist1,myexp1,myexp2,"N")
#'mynorm <- NVTnormalize(mynvt)
#'
#'NVTadvancedmaplot(mynorm,2,2,2)
NVTadvancedmaplot <- function(NVTdataobj,p_cex=1,t_cex=1,l_cex=1) {

  if(length(NVTdataobj@hklist) == 1){
    stop("Only one element in the housekeeping-gene-list, can not calculate linear model")
  }

  if (requireNamespace("ggplot2", quietly = TRUE)  ) {

    if(check_expression_list(NVTdataobj@exp1) && check_expression_list(NVTdataobj@exp2)
       && check_hkgene_list(NVTdataobj@hklist) && check_method(NVTdataobj@norm_method)
       && check_norm_list(NVTdataobj@norm1)  && check_norm_list(NVTdataobj@norm2)
       && exists_hkgene_list(NVTdataobj,NVTdataobj@hklist)){
      l1 <- NVTdataobj@norm1[[1]]
      l2 <- NVTdataobj@norm2[[1]]

      #mean
      l <- cbind(l1,l2)
      mean <- log2(rowMeans(l))

      #log2 fold change
      fc <- log2(l1/l2)

      names(mean) <- rownames(NVTdataobj@norm1)
      names(fc) <- rownames(NVTdataobj@norm1)
      l <- cbind(fc,mean)

      #clean up
      idx <- apply(l, 1, function(x) all(!is.infinite(x)))
      l <- l[idx,]
      idx <- apply(l, 1, function(x) all(!is.na(x)))
      l <- l[idx,]
      idx <- apply(l, 1, function(x) all(!is.nan(x)))
      l <- l[idx,]

      #only houskeeping genes
      m <- cbind(fc,mean)[NVTdataobj@hklist,]
      idx <- apply(m, 1, function(x) all(!is.infinite(x)))
      m <- m[idx,]
      idx <- apply(m, 1, function(x) all(!is.na(x)))
      m <- m[idx,]
      idx <- apply(m, 1, function(x) all(!is.nan(x)))
      m <- m[idx,]

      myl <- as.data.frame(l)

      mysubl <-  myl[NVTdataobj@hklist,]
      myl[,3] <- rownames(myl)
      colnames(myl)[3] <- "names"

      mysubdata=subset(myl,names %in% NVTdataobj@hklist)
      rownames(mysubdata) <- NULL

      d <-  ggplot2::ggplot(data=myl, ggplot2::aes(x = mean, y = fc)) + ggplot2::geom_point( alpha = 0.2,  color="gray", cex=p_cex)
      d <- d + ggplot2::geom_smooth(data=mysubdata, ggplot2::aes(y =fc, x = mean), method=lm, fullrange=TRUE, se=FALSE, alpha = 0.1, color="red", lwd = 0.5)
      d <- d + ggplot2::geom_point(data=mysubdata, ggplot2::aes(y = fc, x = mean), alpha = 0.8, color="blue", cex=p_cex)

      d <- d + ggplot2::geom_text(data=mysubdata, ggplot2::aes(label =names), hjust=-0.2, vjust=0.5, color="blue", cex=l_cex, angle = 320)

      d <- d + ggplot2::ggtitle(bquote(atop(.(paste("MA-plot of:", names(NVTdataobj@norm1),"vs",names(NVTdataobj@norm2))), atop(italic(.(paste(NVTdataobj@norm_method_name,"normalized"))), ""))))
      d <- d + ggplot2::xlab(paste("log2( mean expression )"))
      d <- d + ggplot2::ylab(paste("log2( normalized expression",names(NVTdataobj@norm1),"/",names(NVTdataobj@norm2),")"))
      d <- d + ggplot2::geom_rug(col="darkred",alpha=.1,position='jitter')
      d <- d + ggplot2::theme(axis.title.x = ggplot2::element_text(size = 10 * t_cex), axis.title.y = ggplot2::element_text(size = 10 * t_cex), plot.title = ggplot2::element_text(size = 10 * t_cex), axis.text.x = ggplot2::element_text(size = 10 * t_cex), axis.text.y = ggplot2::element_text(size = 10 * t_cex) )
      d
    }else{
      stop("Not a valid NVTdata object with normalized values!")
    }
  }else{
    print("ggplot2 not found, please use the NVTmaplot function!")
  }
}

#'Returns the linear model a normalized NVTobject
#'
#'@export
#'@param NVTdataobj A previously initialized and normalized NVTobject
#'@return Returns the linear model a normalized NVTobject
#'@examples
#'library("NVT")
#'data(myexp1)
#'data(myexp2)
#'data(mylen)
#'mylist1<-c("ENSG00000111640","ENSG00000163631","ENSG00000075624","ENSG00000172053",
#'"ENSG00000165704","ENSG00000196839","ENSG00000168938","ENSG00000177700")
#'
#'mynvt <- NVTinit(mylist1,myexp1,myexp2,"N")
#'mynorm <- NVTnormalize(mynvt)
#'
#'mylm <- NVTlm(mynorm)
#'
#'summary(mylm)
NVTlm <- function(NVTdataobj) {

  if(length(NVTdataobj@hklist) == 1){
    stop("Only one element in the housekeeping-gene-list, can not calculate linear model")
  }

  if(check_expression_list(NVTdataobj@exp1) && check_expression_list(NVTdataobj@exp2)
     && check_hkgene_list(NVTdataobj@hklist) && check_method(NVTdataobj@norm_method)
     && check_norm_list(NVTdataobj@norm1)  && check_norm_list(NVTdataobj@norm2)
     && exists_hkgene_list(NVTdataobj,NVTdataobj@hklist)){

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
    sample1 <- m[,1]
    sample2 <- m[,2]
    fm <- lm(sample2 ~ sample1)
    return(fm)

  }else{
    stop("Not a valid NVTdata object with normalized values!")
  }
}


#'Calculate the pearson correclation of the housekeeping genes of an initialized and normalized NVTobject
#'
#'@export
#'@param NVTdataobj A previously initialized and normalized NVTobject
#'@return Pearson correlation of the normalized housekeeping genes between the samples
#'@examples
#'library("NVT")
#'data(myexp1)
#'data(myexp2)
#'data(mylen)
#'mylist1<-c("ENSG00000111640","ENSG00000163631","ENSG00000075624","ENSG00000172053",
#'"ENSG00000165704","ENSG00000196839","ENSG00000168938","ENSG00000177700")
#'
#'mynvt <- NVTinit(mylist1,myexp1,myexp2,"TMM")
#'mynorm <- NVTnormalize(mynvt)
#'
#'NVTpearson(mynorm)
NVTpearson <- function(NVTdataobj) {

  if(length(NVTdataobj@hklist) == 1){
    stop("Only one element in the housekeeping-gene-list, can not calculate Pearson correlation")
  }

  if(check_expression_list(NVTdataobj@exp1) && check_expression_list(NVTdataobj@exp2)
     && check_hkgene_list(NVTdataobj@hklist) && check_method(NVTdataobj@norm_method)
     && check_norm_list(NVTdataobj@norm1)  && check_norm_list(NVTdataobj@norm2)
     && exists_hkgene_list(NVTdataobj,NVTdataobj@hklist)){
    if(NVTdataobj@is_norm ){

      m1 <- NVTdataobj@norm1[1][NVTdataobj@hklist,]
      m2 <- NVTdataobj@norm2[1][NVTdataobj@hklist,]

      #pearson <- cor(m1,m2,method="pearson")
      pearson <- cor.test(m1,m2,method="pearson")
      names(pearson$p.value) <- "p-value"
      names(pearson$estimate) <- "pearson"

      #return(pearson)
      return(c(pearson$estimate,pearson$p.value))

    }else{
      return(c("NA","NA"))
    }

  }else{
    stop("Not a valid NVTdata object with normalized values!")
  }
}

#'Calculate the root mean squared deviation (RMSD) of the housekeeping genes of an initialized and normalized NVTobject
#'
#'@export
#'@param NVTdataobj A previously initialized and normalized NVTobject
#'@return Root mean squared deviation (RMSD) of the normalized housekeeping genes between the samples
#'@examples
#'library("NVT")
#'data(myexp1)
#'data(myexp2)
#'data(mylen)
#'mylist1<-c("ENSG00000111640","ENSG00000163631","ENSG00000075624","ENSG00000172053",
#'"ENSG00000165704","ENSG00000196839","ENSG00000168938","ENSG00000177700")
#'
#'mynvt <- NVTinit(mylist1,myexp1,myexp2,"TMM")
#'mynorm <- NVTnormalize(mynvt)
#'
#'NVTrmsd(mynorm)
NVTrmsd <- function(NVTdataobj) {
  if(check_expression_list(NVTdataobj@exp1) && check_expression_list(NVTdataobj@exp2)
     && check_hkgene_list(NVTdataobj@hklist) && check_method(NVTdataobj@norm_method)
     && check_norm_list(NVTdataobj@norm1)  && check_norm_list(NVTdataobj@norm2)
     && exists_hkgene_list(NVTdataobj,NVTdataobj@hklist)){
    if(NVTdataobj@is_norm ){

      m1 <- NVTdataobj@norm1[1][NVTdataobj@hklist,]
      m2 <- NVTdataobj@norm2[1][NVTdataobj@hklist,]

      mydif <- m1 - m2
      rmsd <- sqrt(mean(mydif^2))

      names(rmsd) <- "RMSD"

      return(rmsd)

    }else{
      return("NA")
    }

  }else{
    stop("Not a valid NVTdata object with normalized values!")
  }
}


#'Calculate the mean absolute error (MAE) of the housekeeping genes of an initialized and normalized NVTobject
#'
#'@export
#'@param NVTdataobj A previously initialized and normalized NVTobject
#'@return Mean absolute error (MAE) of the normalized housekeeping genes between the samples
#'@examples
#'library("NVT")
#'data(myexp1)
#'data(myexp2)
#'data(mylen)
#'mylist1<-c("ENSG00000111640","ENSG00000163631","ENSG00000075624","ENSG00000172053",
#'"ENSG00000165704","ENSG00000196839","ENSG00000168938","ENSG00000177700")
#'
#'mynvt <- NVTinit(mylist1,myexp1,myexp2,"TMM")
#'mynorm <- NVTnormalize(mynvt)
#'
#'NVTmae(mynorm)
NVTmae <- function(NVTdataobj) {
  if(check_expression_list(NVTdataobj@exp1) && check_expression_list(NVTdataobj@exp2)
     && check_hkgene_list(NVTdataobj@hklist) && check_method(NVTdataobj@norm_method)
     && check_norm_list(NVTdataobj@norm1)  && check_norm_list(NVTdataobj@norm2)
     && exists_hkgene_list(NVTdataobj,NVTdataobj@hklist)){
    if(NVTdataobj@is_norm ){

      m1 <- NVTdataobj@norm1[1][NVTdataobj@hklist,]
      m2 <- NVTdataobj@norm2[1][NVTdataobj@hklist,]

      mydif <- m1 - m2
      mae <- mean(abs(mydif))

      names(mae) <- "MAE"

      return(mae)

    }else{
      return("NA")
    }

  }else{
    stop("Not a valid NVTdata object with normalized values!")
  }
}

#'Calculate the chosen correclation of the housekeeping genes of an initialized NVTobject for all normalization methods
#'
#'@export
#'@param NVTdataobj A previously initialzed and normalized NVTobject
#'@param cmethod Pearson correlation (p), root mean square deviation (rmsd) or mean absolute error (mae) [default="p"]
#'@param verbose mode on or off (T/F) [default=TRUE]
#'@return Sorted pearson correlation coefficients, RMSD or MAE of the normalized housekeeping genes between the samples
#'@examples
#'library("NVT")
#'data(myexp1)
#'data(myexp2)
#'data(mylen)
#'mylist1<-c("ENSG00000111640","ENSG00000163631","ENSG00000075624","ENSG00000172053",
#'"ENSG00000165704","ENSG00000196839","ENSG00000168938","ENSG00000177700")
#'
#'mynvt <- NVTinit(mylist1,myexp1,myexp2,"N",mylen)
#'
#'NVTtestall(mynvt,"p")
NVTtestall <- function(NVTdataobj, cmethod="p", verbose=TRUE) {
  tmpNVT <- NVTdataobj
  first=TRUE
  check_cmethod(cmethod)

  #test all methods and extract correlation value
  for (n in method_v) {
    tmpNVT@norm_method <- n
    tmpnorm <- NVTnormalize(tmpNVT,verbose = verbose)

    switch(cmethod,
           p={
             p <- NVTpearson(tmpnorm)
           },
           rmsd={
             p <- NVTrmsd(tmpnorm)
           },
           mae={
             p <- NVTmae(tmpnorm)
           },
           {
             stop("Invalid correlation calculation method!")
           }
    )

    if(first){
      pf <- p
      first=FALSE
    }else{
      pf <- rbind(pf,p)
    }
  }

  rownames(pf) <- method_v

  if(cmethod == "p"){#order pearson
    results <- pf[order(pf[,1,drop=FALSE],pf[,2,drop=FALSE],decreasing = T),]
  }else{#order other methods with just one value
    results <- pf[order(pf[,1,drop=FALSE],decreasing = F),]
  }

  return(results)
}


#'Load a gff2 (gtf) or gff3 file, this may take a while depending on the file-size / number of features
#'
#'@export
#'@param gff_file The annotation in gff format
#'@param gff_version The version of the provided gff file [gff1,gff2,gff3,gtf]
#'@param gff_feature The feature to use [default: exon]
#'@param gff_name The name to use [default: gene_id]
#'@return List of gene/exon names and their length
#'@examples
#'library("NVT")
#'
#'#get test GFF-file provided by this library
#'mygffpath<-system.file("extdata", "Ctr-D-UW3CX.gff", package = "NVT")
#'
#'mylen <- NVTloadgff(mygffpath,"gff3","gene","locus_tag")
NVTloadgff <- function(gff_file, gff_version, gff_feature, gff_name) {

  if (requireNamespace("GenomicRanges", quietly = TRUE) && requireNamespace("rtracklayer", quietly = TRUE) && requireNamespace("S4Vectors", quietly = TRUE)) {

    if(is.null(gff_feature)) gff_feature <- "exon"
    if(is.null(gff_name)) gff_name <- "gene_id"

    mygff <- rtracklayer::import(gff_file, format=gff_version, feature.type=gff_feature )
    #mygff <- rtracklayer::import.gff(gff_file, format=gff_version, feature.type=gff_feature )

    mygred <- GenomicRanges::reduce(S4Vectors::split(mygff, GenomicRanges::elementMetadata(mygff)[[gff_name]]))

    mygffred <- GenomicRanges::unlist(mygred, use.names=T)
    GenomicRanges::elementMetadata(mygffred)$gene_id <- rep(names(mygred), S4Vectors::elementLengths(mygred))
    GenomicRanges::elementMetadata(mygffred)$widths <- GenomicRanges::width(mygffred)

    get_length <- function(x) { sum(GenomicRanges::elementMetadata(x)$widths) }

    mylen <- t(sapply(S4Vectors::split(mygffred, GenomicRanges::elementMetadata(mygffred)$gene_id), get_length))
    mylen <- t(mylen)
    colnames(mylen) <- c("Length")

    return(as.data.frame(mylen))

  }else{
    print("GenomicRanges and/or rtracklayer not found, please install these packages or extract gene length manually!")
  }
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

#'Check if housekeeping gene list is valid
#'
#'@param NVTobject A initialized NVTdata object
#'@param hkgene_list Housekeeping gene list
#'@return true or false
exists_hkgene_list <- function(NVTobject,hkgene_list) {
  for ( gene in hkgene_list ) {
    if((gene %in% rownames(NVTobject@exp1)) && (gene %in% rownames(NVTobject@exp2))){

    }else{
      stop(paste("Housekeeping-gene",gene,"is not present in expression list!"))
    }
  }
  return(TRUE)
}

#'Check expression gene list
#'
#'@param exp_list Expression gene list
#'@return true or false
check_expression_list <- function(exp_list) {
  if( is.data.frame(exp_list) || is.list(exp_list) || is.vector(exp_list) || is.matrix(exp_list)){
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
#'@param norm_method Normalization method
#'@return true or false
check_method <- function(norm_method) {
  if( norm_method %in% method_v ){
    return(TRUE)
  }else{
    stop("Unknown method specified!")
  }
}

#'Check correlation method
#'
#'@param c_method Normalization method
#'@return true or false
check_cmethod <- function(c_method) {
  if( c_method %in% method_c ){
    return(TRUE)
  }else{
    stop("Unknown method specified!")
  }
}

#' List of length of each gene
#'
#' A data set with each gene name and its length in nucleotides
#'
#' @format A data frame with 64102 rows and 1 variable:
#' \itemize{
#'   \item{Length}
#' }
"mylen"

#' List of expression data of each gene
#'
#' A data set with each gene name and its expression
#'
#' @format A data frame with 64102 rows and 1 variable:
#' \itemize{
#'   \item{GSM1275862}
#' }
"myexp1"

#' List of expression data of each gene
#'
#' A data set with each gene name and its expression
#'
#' @format A data frame with 64102 rows and 1 variable:
#' \itemize{
#'   \item{GSM1275863}
#' }
"myexp2"

#' List of human housekeeping genes
#'
#' A data set of 10 human housekeeping genes
#'
#' @format A list of 10 elements:
#' \itemize{
#'   \item{HK_genes}
#' }
"mylist_hs"

#' List of ensembl human gene-length
#'
#' A data set with each human ensemble gene name and its length in nucleotides
#'
#' @format A list of 60619 elements:
#' \itemize{
#'   \item{Length}
#' }
"myusecaselen"

#' Dataframe of expression of six human RNA-Seq samples
#'
#' A data set with each gene name and its expression in the six samples
#'
#' @format A data frame with 60619 rows and 6 variables:
#' \itemize{
#'   \item{GSM1464282}
#'   \item{GSM1464283}
#'   \item{GSM1464284}
#'   \item{GSM1464289}
#'   \item{GSM1464290}
#'   \item{GSM1464291}
#' }
"myusecaseexp"
