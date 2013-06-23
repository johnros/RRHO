\name{pvalRRHO}
\alias{pvalRRHO}
\title{
	Compute the significance of the overlap
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pvalRRHO(RRHO.obj, replications, FUN = max)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{RRHO.obj}{
  The output object of the \code{\link{RRHO}} function.
}
  \item{replications}{
  The number of samples from the distribution of the aggregating test statistic.
}
  \item{FUN}{
  The statistic aggregating the infomation from the whole overlap matrix. Typically the \eqn{min} pvalue, or \eqn{max} on \eqn{-log(pval)} scale.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\author{
Jonathan Rosenblatt
}
\note{
Might take a long time to run. Depending on the number of replications, the item (gene) count and the \code{stepsize}.
}

\seealso{
\code{\link{RRHO}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

\keyword{htest}

