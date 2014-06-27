\name{orthopoly}
\alias{orthopoly}
\title{
Orthogonal polynomials}
\description{
This function is called  from the function
\code{omca}. 
It allows the analyst to compute the orthogonal polynomials of each ordered categorical variable. 
The number of the polynomials is equal to the variable category less one.
The function computes the polynomial transformation of the ordered categorical variable.}
\usage{
orthopoly(marginals, scores)
}
\arguments{
  \item{scores}{
 The ordered scores of an ordered variable. By default \code{mj=NULL}, the natural scores (1,2,...) are computed.
}
  \item{marginals}{
 The marginals, relative frequencies of the ordered variable. 
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Describe the value returned
\item{B}{the matrix of the orthogonal polynomials without the trivial polynomial.}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Beh EJ and Lombardo R 2014 Correspondence analysis, Theory, Practice and New Strategies. Wiley.
}
\author{
Rosaria Lombardo and Eric J Beh
}
\note{
Note that the sum of the marginals of the ordered variables should be one.
At the end, the various polynomial matrices will be stored in a super-diagonal matrix. 
}


%% ~Make other sections like Warning with \section{Warning }{....} ~
%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
orthopoly(marginals=c(.1,.2,.3,.2,.2), scores=c(1,2,3,4,5))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
\keyword{nonparametric}% __ONLY ONE__ keyword per line
