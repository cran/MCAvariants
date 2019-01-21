\name{plot.MCAvariants}
\alias{plot.MCAvariants}
\title{Main plot function for classical and ordered multiple correspondence analysis}
\description{
This function allows the analyst to produce the suitable graphical displays with respect to the  classical and ordered multiple
correspondence analysis. 
The main plot function called from the main function \code{MCAvariants}. 
It produces classical graphical displays for \code{catype = "mca"} and \code{catype = "omca"}.\cr
}
\usage{
\S3method{plot}{MCAvariants}(x, catype = "mca", firstaxis = 1, lastaxis = 2, cex = 0.8, 
cex.lab = 0.8, prop = 1, M=2,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Represents the set of the output parameters of the main function \code{MCAvariants} of the \code{R} 
object class \code{mcacorporateris}.
}
  \item{catype}{
The input parameter specifying what variant of correspondence analysis is requested.
}
  \item{firstaxis}{
The dimension reflected along  the horizontal axis.
}
  \item{lastaxis}{
The dimension reflected along the vertical axis.
}
  \item{cex}{
The parameter that specifies the size of character labels of points in graphical displays. By default, it is equal to 1.
}
\item{cex.lab}{
The parameter \code{cex.lab} that specifies the size of character labels of axes in graphical displays. By default, \code{cex.lab = 0.8}.
}
  \item{prop}{
The scaling parameter for specifying the limits of the plotting area. By default, it is equal to 1. 
}
\item{M}{
The number of axes \code{M} considered when portraying  the elliptical confidence regions.  \cr
By default, it is equal to \code{M = 2}. 
}
\item{...}{Further arguments passed to or from other methods.}
}
\details{
It produces classical  graphical displays. Further when \code{catype} is equal to \code{"omca"},  
the individual clusters are portrayed.
}
\references{
Lombardo R and Meulman JJ (2010) Journal of Classification, 27, 191-210. \cr
Beh EJ Lombardo R (2014) Correspondence Analysis, Theory, Practice and New Strategies. Wiley}
\author{
Rosaria Lombardo and Eric J Beh}
\examples{
data(satisfaction)
ris1=MCAvariants(satisfaction, catype = "mca")
plot.MCAvariants(ris1)
print.MCAvariants(ris1)
ris2=MCAvariants(satisfaction, catype = "omca")
plot.MCAvariants(ris2)
print.MCAvariants(ris2)
}
\keyword{ graphical }% __ONLY ONE__ keyword per line
