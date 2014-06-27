\name{mcafun}
\alias{mcafun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Classical multiple correspondence analysis 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function is used in the  secondary function \code{mcabasic} when 
the input parameter of \code{MCAvariants} is \code{catype="mca"}. 
It performs the singular value decomposition of the weighted super-indicator matrix and
compute principal axes, coordinates, weights of rows and columns and
total inertia.

}
\usage{
mcafun(XO, Burt, np, idr, idc, nmod)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{XO}{
The super-indicator data table.}
  \item{Burt}{
The Burt data table.}
  \item{np}{
The number of categorical variables.}
  \item{idr}{
The row labels of data table.}
  \item{idc}{
The column labels of data table.}
  \item{nmod}{
The category number of each variable.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
%\examples{
%    list(xo = xo, xc = xc, autovet = autovet, autovetn = autovetn, 
%        values = sing, valuesn = valuesn, pc1 = pc1, pc0 = pc0, 
%        pc2 = pc2, dj = dj, totin = totin, tot = tot, sing = sing, 
%        nr = nr, Burt = Burt, Raxes = autovet, Caxes = autovetn[, 
%            -1], mu = sing[-1], valuesn = valuesn, R = pc1, C = pc2, 
%        Rweights = uni1, Cweights = uni2)
%  }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
\keyword{multivariate }% __ONLY ONE__ keyword per line
