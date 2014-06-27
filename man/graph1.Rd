\name{graph1}
\alias{graph1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Secondary graphical function}
\description{
This function produces a classical plot when the input parameter of the main code \code{MCAvariants} 
has been set as \code{plottype="classic"}.
This function allows the analyst to plot category and individual points in graphical displays, using symbols and labels for each point. 
It is called internally from the main graphical function \code{plotmcacorporateris}.}
\usage{
graph1(f, prop = 1, cex = 1,  cex.lab = 0.8, cols = c(2, 4), 
 main = "Plot", mar = c(5, 4, 4, 2) + 0.1, a1, a2, inertiapc, Inames)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{f}{
The coordinates of the variable category or the coordinates of the individuals (objects)}
  \item{prop}{
The scaling parameter for specifying the limits of the plotting area. By default, it is equal to 1. 
}
  \item{cex}{
The parameter that specifies the size of character labels of points in graphical displays. By default, it is equal to 1.
}
  
  \item{cex.lab}{ The parameter for setting the size of character labels of axes in graphical displays.
}
  \item{cols}{The parameter for setting the column size of data table}
 
  \item{main}{
The tite of the plot.
}
  \item{mar}{
The parameter for setting the size of plot area in graphical displays.
}
  \item{a1}{
The dimension reflected along the horizontal (polynomial or principal) axis in the graphical display.
}
  \item{a2}{
The dimension reflected along the vertical (polynomial or principal) axis in the graphical display.}
  \item{inertiapc}{
The percentage of explained inertia.
}
  \item{Inames}{
The labels of the row variable.}
}
\details{
This graphical function is used internally by the main graphical function \code{plotmcacorporateris}.}
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
%This graphical function is used by the main graphical function printmcacorporateris
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
% R documentation directory.
\keyword{ graphic }
