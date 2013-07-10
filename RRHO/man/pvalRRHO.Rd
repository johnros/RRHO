\name{pvalRRHO}
\alias{pvalRRHO}
\title{
	Compute the significance of the overlap between two lists
}
\description{
 Computes the significance of the agreements between lists as returned by \code{\link{RRHO}} using resampling.
}
\usage{
pvalRRHO(RRHO.obj, replications, stepsize=RRHO.obj$stepsize, FUN= max)
}

\arguments{
  \item{RRHO.obj}{
  The output object of the \code{\link{RRHO}} function.
}
  \item{replications}{
  The number of samples to be taken from the distribution of the aggregated test statistic.
}
\item{stepsize}{
  Controls the resolution of the test: how many items between any two overlap tests (i.e., netween any two $i$s and two $j$s.)
}
  \item{FUN}{
  The function aggregating infomation from the whole overlap matrix into one summary statistic. Typically the \eqn{min} pvalue, or \eqn{max} on \eqn{-log(pval)} scale.
}
}
\details{
  The distribution of \eqn{FUN(-log(pval))} is computed using resampling. 
  
  The aggregating function will typically be the \code{max} function, corresponding to the maximal \code{-log(pvalue)}, i.e., the most significant agreement over all sublists. 
  
  The distribution is computed by resampling pairs of null sequences, computing the significances of all the overlaps as done in the reference, applying the aggregating function supplied by the user, and returning the permutation based significance.
}
\value{
\item{pval}{The FWER corrected significance of observed aggregated pvalue.}
\item{FUN.ecdf}{The simulated sampling distribution of the aggregated pvalues. }
\item{FUN}{The matrix aggregation function used. typicall \code{max} for minimal p-value.}
\item{n.items}{Length of lists.}
\item{stepsize}{See \code{\link{RRHO}}}
\item{replications}{The number of simulation replications.}
\item{call}{The function call.}
}

\author{
Jonathan Rosenblatt
}
\note{
Might take a long time to run. Depending on the number of \code{replications}, the item (gene) count and the \code{stepsize}.

Also note that the significance returned is a conservative value (by a constant of 1/\code{replications}). 
}

\seealso{
\code{\link{RRHO}}
}
\examples{
list.length <- 100
list.names <- paste('Gene',1:list.length, sep='')
gene.list1<- data.frame(list.names, sample(list.length))
gene.list2<- data.frame(list.names, sample(list.length))
RRHO.example <-  RRHO(gene.list1, gene.list2)
pval.testing <- pvalRRHO(RRHO.example,50) 
}

\keyword{htest}

