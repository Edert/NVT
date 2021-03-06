% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NVT-classes.R
\docType{methods}
\name{show}
\alias{show}
\alias{show,NVTdata-method}
\title{Show method, to retrieve the normalized expression values}
\usage{
\S4method{show}{NVTdata}(object)
}
\arguments{
\item{object}{A previously initialized and normalized NVTobject}
}
\value{
Returns the normalized expression values of both samples
}
\description{
Show method, to retrieve the normalized expression values
}
\examples{
library("NVT")
data(myexp1)
data(myexp2)
data(mylen)
mylist1<-c("ENSG00000111640","ENSG00000163631","ENSG00000075624","ENSG00000172053",
"ENSG00000170950","ENSG00000165704","ENSG00000196839","ENSG00000168938","ENSG00000177700")

mynvt <- NVTinit(mylist1,myexp1,myexp2,"N")
mynorm <- NVTnormalize(mynvt)

show(mynorm)
}

