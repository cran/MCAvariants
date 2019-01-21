\name{tableconvert}
\alias{tableconvert}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Convert contingency table in table of reduced code
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This simple piece of R code converts a two-way or three-way 
contingency table into what is required to analyse MCA (table of reduced code: n by number of variables).
}
\usage{
tableconvert(N)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
A two-way or three-way contingency table to convert in a table n by np, 
where np is the number of the categorical variables.}
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
Beh EJ Lombardo R (2014) Correspondence Analysis, Theory, Practice and New Strategies. Wiley}
\author{
Rosaria Lombardo and Eric J Beh}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
alligatormca<-tableconvert(alligator.dat)
dimnames(alligatormca)<-list(paste("a", 1:300,sep = ""),c("Size","Food","Lake"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
\keyword{ multivariate}% __ONLY ONE__ keyword per line
