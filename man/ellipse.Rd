\name{ellipse}
\alias{ellipse}
\title{
Ellipse}
\description{
This function is called internally from the graphical functions \code{caellipse} and \code{nscaellipse}.
It constructs the algebraic ellipses of confidence.
}
%\usage{
%ellipse(hlaxa = 1, hlaxb = 1, theta = 0, xc = 0, yc = 0, npoints = 100, col)
%}
\arguments{
  \item{hlaxa}{
The length of the horizontal ellipse semi-axis.}
  \item{hlaxb}{
The length of the vertical ellipse semi-axis.}
%  \item{theta}{
%%     ~~Describe \code{theta} here~~
%}
  \item{xc}{
The coordinate of the generic point on the horizontal axis.}
  \item{yc}{
The coordinate of the generic point on the vertical axis.}
%  \item{npoints}{
%%     ~~Describe \code{npoints} here~~
%}
 \item{col}{
The color of ellipses.
}
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
Beh EJ and Lombardo R 2014 Correspondence Analysis: Theory, Practice and New Strategies. John Wiley & Sons. \cr
Macdonald PDM 2002 Drawing an ellipse in Splus or R. \cr
Available at www.math.mcmaster.ca/peter/s4m03/s4m03_0102/ellipse.html (accessed August 18, 2016).
}
\author{
Rosaria Lombardo and Eric J Beh
}
\note{
This function is called from the secondary graphical function \code{caellipse} or \code{nscaellipse}, which is 
called from the main plot function \code{plot.CAvariants} and it can be executed 
for all variants of correspondence analysis.}
%% ~Make other sections like Warning with \section{Warning }{....} ~
%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
%\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
%## The function is currently defined as
%function (hlaxa = 1, hlaxb = 1, theta = 0, xc = 0, yc = 0, npoints = 100, 
%    col) 
%{
%    xp <- NULL
%    yp <- NULL
%    for (i in 0:npoints) {
%        a <- (2 * pi * i)/npoints
%        x <- hlaxa * cos(a)
%        y <- hlaxb * sin(a)
%        if (theta != 0) {
%            alpha <- angle(x, y)
%            rad <- sqrt(x^2 + y^2)
%            x <- rad * cos(alpha + theta)
%            y <- rad * sin(alpha + theta)
%        }
%        xp <- c(xp, x + xc)
%        yp <- c(yp, y + yc)
%    }
%    lines(xp, yp, type = "l", col = col)
%    invisible()
%  }
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
\keyword{ univariate }% __ONLY ONE__ keyword per line
