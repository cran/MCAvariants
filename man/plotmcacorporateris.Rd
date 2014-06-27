\name{plotmcacorporateris}
\alias{plotmcacorporateris}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Main plot function}
\description{
This function allows the analyst to produce the suitable graphical displays with respect to the  classical and ordered multiple
correspondence analysis. 
The main plot function called from the main function \code{MCAvariants}. 
It produces classical graphical displays for \code{catype = "mca"} and \code{catype = "omca"}.\cr
}
\usage{
plotmcacorporateris(x, catype = "mca", firstaxis = 1, lastaxis = 2, cex = 0.8, prop = 1)
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
  \item{prop}{
The scaling parameter for specifying the limits of the plotting area. By default, it is equal to 1. 
}
%  \item{plottype}{
%The input parameter specifying what kind of graphical display is requested. It is equal 
%to \code{plottype = "classic"} as, at the moment, only the classical plot is allowed.
%}
}
\details{
It produces classical  graphical displays. Further when \code{catype} is equal to \cr
\code{"omca"},  
the individual clusters are portrayed.\cr}
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
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graphical }% __ONLY ONE__ keyword per line
