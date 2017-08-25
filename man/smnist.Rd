\name{smnist}
\alias{smnist}
\docType{data}
\title{
A subsampled MNIST data
}
\description{
Given MNIST training data, 1000 images are randomly sampled with their labels.
}
\usage{data("smnist")}
\format{
  A list with following two object.
  \describe{
    \item{\code{strain}}{a data frame with 1000 observations on 784 features (pixels).}
    \item{\code{slab}}{a numeric label vector of length of 1000.}
  }
}
\source{
http://yann.lecun.com/exdb/mnist/
}
\references{
Y. LeCun, L. Bottou, Y. Bengio, and P. Haffner. "Gradient-based learning applied to document recognition." Proceedings of the IEEE, 86(11):2278-2324, November 1998.
}
\examples{
data(smnist)
## maybe str(smnist) ; ...
}
\keyword{datasets}
