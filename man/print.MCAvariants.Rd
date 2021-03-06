\name{print.MCAvariants}
\alias{print.MCAvariants}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Main printing function}
\description{
This function prints results of classical or ordered multiple correspondence analysis. 
The input parameter is the name of the output of the main function  \code{MCAvariants}. 
}
\usage{
\S3method{print}{MCAvariants}(x, printdims = 2,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{
The output of the main function \code{CAvariants}.}
\item{printdims}{The number of dimensions, \code{printdims}, that are used to generate the correspondence plot 
and for summarising the numerical output of the analysis. By default, \code{printdims  =  2}.}
\item{...}{Further arguments passed to or from other methods.}
}
\details{
This function uses another function (called \code{printwithaxes}) for specifying the number of matrix dimensions to print. 
}
\value{
The value of output returned depends on the kind of multiple correspondence analysis performed.
\item{DataTable }{The Burt data table.}
\item{Row coordinates}{Rows in principal coordinates: the first 10.}
\item{Column coordinates }{ Column in principal coordinates.}
\item{Polynomials}{Polynomial functions of each variable. 
When \code{catype} is \code{omca}.}
\item{Linear Percentage of Clusters }{The percentage of individuals belonging to each cluster. When \code{catype} is \code{omca}.}
\item{Polynomial Components  of Total Inertia}{The decomposition of total inertia via orthogonal polynomials. When \code{catype} is \code{omca}. }
\item{Degree of Freedom}{Degree of Freedom of Polynomial Component.  When \code{catype} is \code{omca}.}
\item{Inertia values}{Inertia values of super-indicator and Burt table.}
\item{Benzecri's Inertia values}{Adjusted Inertia values, percentages and cumulative.}
\item{Total Degree of Freedom}{The degree of freedom of total inertia.}
\item{Total inertia of X}{Total inertia of Super-Indicator table}
\item{Total inertia of B}{Total inertia of BURT table.}
\item{Chi-square values}{Chi-square values of  BURT Inertia.}
\item{Total Chi-square values}{Chi-square values of total Inertia of Burt table.}
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
%##---- Should be DIRECTLY executable !! ----
res=MCAvariants(satisfaction, catype = "omca", np = 5, vordered=c(TRUE,TRUE,TRUE,TRUE,TRUE))
print(res)
   }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
\keyword{ multivariate }% __ONLY ONE__ keyword per line
