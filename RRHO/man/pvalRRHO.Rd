\name{pvalRRHO}
\alias{pvalRRHO}
\title{
	Compute the significance of the overlap
}
\description{
 Computes the significance of the agreements between lists as treturned by \code{\link{RRHO}} using resampling.
}
\usage{
pvalRRHO(RRHO.obj, replications, stepsize=RRHO.obj$stepsize, FUN= max)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{RRHO.obj}{
  The output object of the \code{\link{RRHO}} function.
}
  \item{replications}{
  The number of samples from the distribution of the aggregating test statistic.
}
\item{stepsize}{
  Controls the resolution of the test: how many items between any two overlap tests. 
}
  \item{FUN}{
  The statistic aggregating the infomation from the whole overlap matrix. Typically the \eqn{min} pvalue, or \eqn{max} on \eqn{-log(pval)} scale.
}
}
\details{
  The distribution of some user supplied function of the RRHO pvalues is computed using simulation. 
  
  The aggregating function will typically be the \code{max} function, corresponding to the maximal \code{-log(pvalue)}, i.e., the most significant agreement over all sublists. 
  
  The distribution is computed by resampling pairs of null sequences, computing the significances of all the overlaps as done in the reference, applying the aggregating function supplied by the user, and returning the permutation based significance.
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

Also note that the significance returned is a conservative value (by a constant of 1/\code{replications}). The conservativeness vanishes are more replication are performed. 
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

