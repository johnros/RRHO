\name{RRHO}
\alias{RRHO}
\alias{pvalRRHO}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Rank-Rank Hypergeometric Overlap Test}

\description{The function tests for significant overlap between two sorted lists using the method in the reference.}

\usage{
RRHO(list1, list2, stepsize = defaultftepSize(list1, list2), labels, plots = FALSE, outputdir = NULL, BY = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{list1}{
	 data.frame. First colume is the element (possibly gene) identifier, and the second is its value.
}
  \item{list2}{
	 data.frame. First colume is the element (possibly gene) identifier, and the second is its value.
}
  \item{stepsize}{
	  Controls the resolution of the test: how many items between any two overlap tests. 
}
  \item{labels}{
	  CCharacter vector with two elements: the labels of the two lists.
}
  \item{plots}{
	  Logical. Should output plots be returned?
}
  \item{outputdir}{
	  Path name where plots ae returned.
}
\item{n.items}{
	Number of items (genes) in lists.
}
\item{stepsize}{
	See in arguments.
}
  \item{BY}{
	  Logical. Should Benjamini-Yekutieli FDR corrected pvalues be computed?
}
}
\details{
	Following the method in the reference, the function computes the number of overlapping elements in the first \eqn{i} and \eqn{j} elements of each list, and return the observed significance of this overlap using a hypergeometric test (see \code{\link{fisher.exact}}).
	The output is returned as a list of matrices including: the overlap in the first \eqn{i,j} elements and the significance of this overlap. Optional outputs include plots of these matrices in .jpg format.
	The reference suggests using the Benjamini-Yekutieli, which is implemented as an option. This option however defaults to FALSE, as the second author of this package feels (weak) FWER control is more appropriate than FDR. This is bacause all the p-values originate from the \emph{same hypothesis} being tested, even if using different test statistics. 
}

\value{
\item{hypermat}{Matrix of \eqn{-log(pvals} of the test for the first \eqn{i,j} elements of the lists.}
\item{hypermat.counts}{Counts of the number of agreements in the first \eqn{i,j} elements of the lists.}
\item{hypermat.by}{An optional output of the B-Y corrected p-values of \code{hypermat} }
}

\references{
Plaisier, Seema B., Richard Taschereau, Justin A. Wong, and Thomas G. Graeber. “Rank–rank Hypergeometric Overlap: Identification of Statistically Significant Overlap Between Gene-expression Signatures.” Nucleic Acids Research 38, no. 17 (September 1, 2010): e169–e169. doi:10.1093/nar/gkq636.
TODO: add Benjamini-Yekutieli
}
\author{
	Jason Stein and Jonathan Rosenblatt
}
% \note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

% \seealso{
% ~~objects to See Also as \code{\link{help}}, ~~~
% }
\examples{
	list.length <- 100
	list.names <- paste('Gene',1:list.length, sep='')
	gene.list1<- data.frame(list.names, sample(100))
	gene.list2<- data.frame(list.names, sample(100))
	RRHO.example <-  RRHO(gene.list1, gene.list2)
	image(RRHO.example$hypermat)
}
\keyword{htest}
