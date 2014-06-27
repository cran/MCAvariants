\name{miocount}
\alias{miocount}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
The counting function
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The function that counts the number of individuals in each clusters automatically generated in ordered multiple correpsondence analysis.
}
\usage{
miocount(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The coordinates of axes}
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
\note{
This function is used in the function \code{omcabasic} when in the main function \code{MCAvariants} 
the input parameter is \code{catype="omca"}. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
%\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
\keyword{ multivariate}% __ONLY ONE__ keyword per line
