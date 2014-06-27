\name{mcabasicresultsclass}
\Rdversion{1.1}
\docType{class}
\alias{mcabasicresultsclass-class}
\title{Class \code{"mcabasicresultsclass"}}
\description{
The Class of objects common to the two variants of multiple correspondence analysis
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{mcabasicresultsclass(...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
\describe{
\item{\code{RX}:}{Object of type \dQuote{matrix}. The weighted super-indicator table.}
\item{\code{CX}:}{Object of type \dQuote{matrix}. The weighted transposed super-indicator table.}
\item{\code{Rweights}:}{Object of type \dQuote{matrix}. The weights of row space.}
\item{\code{Cweights}:}{Object of type \dQuote{matrix}. The weights of column space.}
\item{\code{nmod}:}{Object of type \dQuote{numeric}. The number of categories of each variables.}
\item{\code{tmod}:}{Object of type \dQuote{numeric}. The total number of categories in the super-indicator table.}
\item{\code{np}:}{Object of type \dQuote{numeric}. The total number of variables.}
\item{\code{Raxes}:}{Object of type \dQuote{matrix}. The axes of row space.}
\item{\code{Caxes}:}{Object of type \dQuote{matrix}. The axes of column space.}
\item{\code{mu}:}{Object of type \dQuote{numeric}. The singular values when \code{catype}=\code{"mca"} or eigenvalues when \code{catype}=\code{"omca"}.}
\item{\code{dj}:}{Object of type \dQuote{matrix}. The supe-diagonal column weight matrix.}
\item{\code{xo}:}{Object of type \dQuote{matrix}. The original data table.}
\item{\code{listBpoly}:}{Object of type \dQuote{list}. The orthogonal polynomials related to each ordered variable
when \code{catype}=\dQuote{omca}.}  
\item{\code{LinearPercentage}:}{Object of type \dQuote{matrix}. The percentage number of the individuals in clusters when
\code{catype}=\dQuote{omca}.}
\item{\code{BURT}:}{Object of type \dQuote{list}. The contingency tables composing the Burt table.}
}
}
\section{Extends}{
Class \code{"\linkS4class{mcacorporaterisclass}"}, directly.
}
\section{Methods}{
No methods defined with class "mcabasicresultsclass" in the signature.
}
\references{
Lombardo R and Meulman JJ (2010) Journal of Classification, 27, 191-210. \cr
Beh EJ Lombardo R (2014) Correspondence Analysis, Theory, Practice and New Strategies. Wiley}
\author{
Rosaria Lombardo}
%\note{
%%  ~~further notes~~
%}
%% ~Make other sections like Warning with \section{Warning }{....} ~
%\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
%}
\examples{
showClass("mcabasicresultsclass")
}
\keyword{classes}

