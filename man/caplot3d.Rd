\name{caplot3d}
\alias{caplot3d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Three dimensional correspondence plot
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function is used in the plot function \code{plot.CAvariants} when the logical parameter is 
\code{plot3d = TRUE}. 
It produces a 3-dimensional visualization of the association. 
}
\usage{
caplot3d(coordR, coordC, inertiaper, firstaxis = 1, lastaxis = 2, thirdaxis = 3)
}
\arguments{
\item{coordR}{
The row principal or standard coordinates.
}
\item{coordC}{
The column principal or standard coordinates.
}
\item{inertiaper}{
The percentage of the total inertia explained inertia by each dimension.
}
\item{firstaxis}{
The first axis number. By default, \code{firstaxis = 1}.
}
\item{lastaxis}{
The second axis number. By default, \code{lastaxis = 2}.
}
\item{thirdaxis}{
The third axis number. By default, \code{thirdaxis = 3}.
}
}
\references{
Beh EJ and Lombardo R 2014 Correspondence Analysis: Theory, Practice and New Strategies. Wiley.\cr
Lombardo R Beh EJ 2016 Variants of Simple Correspondence Analysis. The R Journal, 8 (2), 167--184.
}
\author{
Rosaria Lombardo and Eric J. Beh
}
\note{
This function depends on the \code{R} library \code{plotly}.
}
%\examples{
%data(asbestos)
%cabasic(asbestos)
%}
\keyword{plot}
