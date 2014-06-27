\name{mcacorporaterisclass-class}
\Rdversion{1.1}
\docType{class}
\alias{mcacorporaterisclass-class}
\alias{mcacorporaterisclass}
\title{Class \code{"mcacorporaterisclass"}}
\description{
The Class of objects to show the graphical and numeric results of one of the  variants of multiple 
correspondence analysis
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{mcacorporaterisclass(...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
\describe{
\item{\code{br}:}{Objects of class \code{"mcabasicresultsclass"}: The objects of the basic type of multiple correspondence analysis performed.}
\item{\code{DataMatrix}:}{Object of type \code{"matrix"}: Data Matrix.}
\item{\code{rows}:}{Object of type \code{"numeric"}: Number of rows of data matrix.}
\item{\code{cols}:}{Object of type \code{"numeric"}: Number of columns of data matrix.}
\item{\code{np}:}{Object of type \code{"numeric"}: The number of variables.}
\item{\code{nmod}:}{Object of type \code{"numeric"}: The number of variable categories.}
\item{\code{rowlabels}:}{Object of type \code{"character"}: The labels of rows.}
\item{\code{collabels}:}{Object of type \code{"character"}: The labels of columns.}
\item{\code{Rprinccoord}:}{Object of type \code{"matrix"}: The row coordinates.}
\item{\code{Cprinccoord}:}{Object of type \code{"matrix"}: The columns coordinates.}
\item{\code{inertiaXsum}:}{Object of type \code{"numeric"}: The total inertia of the super-indicator table.}
\item{\code{inertiaBurtsum}:}{Object of type \code{"numeric"}: The total inertia of the Burt table.}
\item{\code{inertiasAdjusted}:}{Object of type \code{"numeric"}: The adjusted inertia values.}
\item{\code{inertias}:}{Object of type \code{"matrix"}: The inertia adjusted inertia values, percentages and cumulative values.}
\item{\code{printdims}:}{Object of type \code{"numeric"}: The  dimension of a matrix in print. By default it is equal to 3.}
\item{\code{maxaxes}:}{Object of type \code{"numeric"}: The maximum number of axes.}
\item{\code{catype}:}{Object of type \code{"character"}: The kind of multiple correspondence analysis chosen, classical or ordered, that is "catype" is \dQuote{mca} or \dQuote{omca}.}
\item{\code{comp}:}{Object of type \code{"matrix"}: The polynomial components of inertia when "mcatype" is \dQuote{omca}. The total inertia is partitioned in terms of polynomial components.}
\item{\code{componentpvalue1}:}{Object of type \code{"matrix"}: The p-value of the polynomial components of total inertia, when "catype" is \dQuote{omca}.}
\item{\code{degreef}:}{Object of type \code{"vector"}: Degree of freedom of polynomial components of total inertia when, "catype" is \dQuote{omca}.}
}
}
\section{Methods}{
No methods defined with class "mcacorporaterisclass" in the signature.
}
\references{
Lombardo R and Meulman JJ (2010) Journal of Classification, 27, 191-210. \cr
Beh EJ and Lombardo R (2014) Correspondence analysis, Theory, Practice and New Strategies. Wiley.
}
\author{
Rosaria Lombardo
}
%\note{
%%  ~~further notes~~
%}
%% ~Make other sections like Warning with \section{Warning }{....} ~
%\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
%}
\examples{
showClass("mcacorporaterisclass")
}
\keyword{classes}
