## ----knitr, echo=FALSE, results="hide"-----------------------------------
library("knitr")
opts_chunk$set(
  tidy=FALSE,
  dev="png",
  fig.show="hide",
  #fig.width=4, fig.height=4.5,
  #cache=TRUE,
  message=FALSE)

## ----style-knitr, eval=TRUE, echo=FALSE, results="asis"---------------------------------
BiocStyle::latex()

## ----installNVT, eval=FALSE-------------------------------------------------------------
#  install.packages(file.path("/home/user/Downloads/","NVT_1.0.tar.gz"),
#  repos=NULL, type="source")

## ----loadNVT----------------------------------------------------------------------------
library("NVT")

## ----loadExData-------------------------------------------------------------------------
data(mylen)
data(myexp1)
data(myexp2)

#show just the first four elements of the loaded data
head(mylen,4)
head(myexp1,4)

## ----createList-------------------------------------------------------------------------
mylist1<-c("ENSG00000111640","ENSG00000163631","ENSG00000075624",
           "ENSG00000172053","ENSG00000165704","ENSG00000196839",
           "ENSG00000177700")

## ----loadLengthData---------------------------------------------------------------------
#this line gets the path to the gff file provided in the NVT package
#(annotation from Chlamydia trachomatis, ACCESSION: NC_000117 )
mygffpath <- system.file("extdata", "Ctr-D-UW3CX.gff", package = "NVT")

#this function loads the gff file from the gffpath
mylen1 <- NVTloadgff(mygffpath,"gff3","gene","locus_tag")

head(mylen1)

## ----generateNVTobj1--------------------------------------------------------------------
mynvt <- NVTinit(mylist1,myexp1,myexp2,"TMM")

## ----generateNVTobj2--------------------------------------------------------------------
mynvt <- NVTinit(mylist1,myexp1,myexp2,"RPKM",mylen)

## ----normlaizeNVTobj--------------------------------------------------------------------
mynorm <- NVTnormalize(mynvt)

## ----getnorm----------------------------------------------------------------------------
mynvalues <- show(mynorm)

head(mynvalues)

## ----simpleplot, dev="pdf", fig.width=3.5, fig.height=3.5-------------------------------
NVTplot(mynorm,0.4)

## ----simplemaplot, dev="pdf", fig.width=3.5, fig.height=3.5-----------------------------
NVTmaplot(mynorm,0.4)

## ----advancedplot, dev="pdf", fig.width=3.5, fig.height=3.5-----------------------------
mynvt <- NVTinit(mylist1,myexp1,myexp2,"TMM",mylen)
mynorm <- NVTnormalize(mynvt)
NVTadvancedplot(mynorm,1,1,1)

## ----advancedmaplot, dev="pdf", fig.width=3.5, fig.height=3.5---------------------------
NVTadvancedmaplot(mynorm,1,1,1)

## ----lm---------------------------------------------------------------------------------
mylm <- NVTlm(mynorm)

summary(mylm)

## ----pearson----------------------------------------------------------------------------
NVTpearson(mynorm)

## ----rmsd-------------------------------------------------------------------------------
NVTrmsd(mynorm)

## ----mae--------------------------------------------------------------------------------
NVTmae(mynorm)

## ----testall----------------------------------------------------------------------------
NVTtestall(mynorm,"p")

## ----sessInfo, results="asis", echo=FALSE-----------------------------------------------
toLatex(sessionInfo())

## ----resetOptions, results="hide", echo=FALSE-------------------------------------------
options(prompt="> ", continue="+ ")

