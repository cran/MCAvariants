\name{mcaellipse}
\alias{mcaellipse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Algebraic elliptical confidence regions 
}
\description{
It produces elliptical confidence regions when multiple correspondence analysis is performed. 
This function allows the analyst to superimpose confidence ellipses into a graphical display when the input parameter \code{catype}
of the main code \code{MCAvariants} is equal to \code{"mca"} or \code{"omca"}. 
It is called internally  from the main plot function \code{plot.MCAvariants}. 
It uses the function \code{ellipse}.
}
\usage{
mcaellipse(Xtable, a1 = 1, a2 = 2, alpha = 0.05, cols = c(2, 4), np=5,
nmod=c(5,5,5,5,5), M = 2, cex = 0.8, cex.lab = 0.8,  
mar = c(5, 4, 4, 2) + 0.1, prop = 0.8, Imass, 
Jmass, a, b, g, fr, dmu, inertiapc, pos = 2,  graphy = TRUE, ell = FALSE)
}
\arguments{
\item{Xtable}{
The two-way contingency table.
}
\item{a1}{
The dimension reflected along the horizontal axis.
}
  \item{a2}{
The dimension reflected along the vertical axis.
}
  \item{alpha}{
The confidence level of the elliptical regions. By default, \code{alpha = 0.05}.
}
  \item{cols}{
The graphical parameter for setting the colours of the points in the graphical displays.
}
  \item{np}{The input parameter for specifying the number of categorical variables.
}
  \item{nmod}{The input parameter (it is a vector) for specifying the number of categories of each categorical variables.
}
  \item{M}{
The number of  axes used when drawing ellypse. By default,  \code{M = 2}. Its maximum value is equal to the rank of the data matrix.}
  \item{cex}{
The parameter for setting the size of character labels of points in graphical displays.  By default,  \code{cex = 0.8}. 
}
  \item{cex.lab}{
The parameter for setting the size of character labels of axes in graphical displays.  By default,  \code{cex.lab = 0.8}.
}
  \item{mar}{
The parameter for setting the size of the plotting area. 
}
  \item{prop}{
The scaling parameter for specifying the limits of the plotting area. By default, \code{prop = 1}. 
}
   \item{Imass}{
The weight matrix of the row variable.
}
  \item{Jmass}{
The weight matrix of the column variable.
}
  \item{a}{
The row principal or polynomial axes.
}
  \item{b}{
The column principal or polynomial axes.
}
  \item{fr}{
The row coordinates.
}
  \item{g}{
The column coordinates.
}
  \item{dmu}{
The squared singular values or inertia of each axis.
}
  \item{inertiapc}{
The percentage of explained inertia.
}
  \item{pos}{
The parameter that specifies the position of  point symbols in the graphical displays. By default, \code{pos = 2}.}
\item{graphy}{The flag parameter used for  diplaying ellipses when using \code{plot.CAvariants} and not displaying when using \code{print.CAvariants} function.
It is set as \code{graphy = TRUE} in \code{plot.CAvariants}  and as \code{graphy = FALSE}  in \code{print.CAvariants}.}
\item{ell}{The flag parameter used for displaying ellipse. By default, \code{ell = TRUE}.}
}
\details{
The output values of this function.
}
\value{
\item{eccentricity}{Value of ellipse eccentricity. This is the distance between the ellipse center and either of its two foci, 
which can be thought of as a measure of how much the conic section deviates from being circular (when  it is equal to zero then the region becomes circular). }
\item{HL Axis 1}{Value of ellipse semi-axis 1 for each row and column points.}
\item{HL Axis 2}{Value of ellipse semi-axis 2 for each row and column points.}
\item{Area}{Ellipse area  for each row and column points.}
\item{pvalcol}{P-value for each row and column points.}
}
\references{
Beh EJ and Lombardo R 2014 Correspondence Analysis: Theory, Practice and New Strategies. John Wiley & Sons.
}
\author{
Rosaria Lombardo and Eric J Beh
}
\note{
This function is called from the main plot function \code{plot.MCAvariants} and is executed when 
 \code{ell = TRUE}.
}
\keyword{ multivariate }% __ONLY ONE__ keyword per line
