\name{alligator.dat}
\alias{alligator.dat}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
Alligator in American lakes }
\description{
The data set is a three-way contingency table. 
It consists of 2 rows (alligators'size), 5 columns (alligators'food) 
by 4 tubes (alligators'lake). The table should be converted in reduced code table,
using the function \code{tableconvert} for getting \code{alligatormca}. 
}
\usage{data(alligator.dat)}
\format{
A data frame with 300 alligators on the following 3 variables.
  \describe{
    \item{Size}{A numeric vector of  categories ranging from 1 to 2 (small and large).}
    \item{Food}{A numeric vector of  categories ranging from 1 to 5 
(type of food: fish, invertebrate, reptile, bird, other.}
    \item{Lake}{a numeric vector of  categories ranging from 1 to 4 for the 
four American lakes: Hancock, Oklawaha, Trafford, George.}
 }
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
Agresti (2007), p. 270\cr
Agresti A and Gottard A 2007 Independence in multiway contingency tables: S.n. roys breakthroughs and later
developments. Journal of Statistical Planning and Inference, 137:3126--3226.
}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(alligator.dat)
#dim(alligator.dat)
#dimnames(alligator.dat)
}
\keyword{datasets}
