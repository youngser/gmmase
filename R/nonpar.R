#'
#' Nonparametric two-sample testing using kernel-based test statistic
#'
#' This is a simple implementation of the kernel-based test statistic for the nonparametric
#' two-sample testing problem of given \eqn{X_1, X_2, \dots, X_n} i.i.d. \eqn{F} and
#' \eqn{Y_1, Y_2, \dots, Y_m} i.i.d. \eqn{G}, test the null hypothesis of \eqn{F = G} against
#' the alternative hypothesis of \eqn{F \not = G}. The test statistic is based on embedding
#' \eqn{F} and \eqn{G} into a reproducing kernel Hilbert space and then compute a distance between
#' the resulting embeddings. For this primitive, the Hilbert space is associated with the
#' Gaussian kernel.
#'
#' @param Xhat1 a \eqn{n} x \eqn{d} matrix
#' @param Xhat2 a \eqn{n} x \eqn{d} matrix
#' @param sigma a bandwidth for the Gaussian kernel
#'
#' @return \code{T} A scalar value \eqn{T} such that \eqn{T} is near 0 if the rows of
#' \eqn{X} and \eqn{Y} are from the same distribution and \eqn{T} far from 0 if the rows of
#' \eqn{X} and \eqn{Y} are from different distribution.
#'
#' @author Youngser Park <youngser@jhu.edu>
#' @export

nonpar <- function(Xhat1,Xhat2,sigma=0.5)
{
    #    Tij <- find.transform(Xhat1, Xhat2)
    #    Sij <- kernel.stat(Xhat1 %*% Tij, Xhat2, sigma)
    Sij1 <- kernel.stat(Xhat1, Xhat2, sigma)
#    Sij2 <- kernel.stat(Xhat1 %*% diag(c(-1,1)), Xhat2, sigma)
#    Sij3 <- kernel.stat(Xhat1 %*% diag(c(1,-1)), Xhat2, sigma)
#    Sij4 <- kernel.stat(Xhat1 %*% diag(c(-1,-1)), Xhat2, sigma)
#    Sij <- min(Sij1, Sij2, Sij3, Sij4)

    return(Sij1)
}

kernel.stat <- function(X,Y,sigma=0.2){
    n <- nrow(X)
    m <- nrow(Y)

    tmpXX <- sum(exp(-(as.matrix(stats::dist(X))^2)/(2*sigma^2))) - n
    tmpYY <- sum(exp(-(as.matrix(stats::dist(Y))^2)/(2*sigma^2))) - m
    tmpXY <- sum(exp(-(rect.dist(X,Y))/(2*sigma^2)))

    tmp <- tmpXX/(n*(n-1)) + tmpYY/(m*(m-1)) - 2*tmpXY/(m*n)

    return((m+n)*tmp)
}

rect.dist <- function(X,Y){
    X <- as.matrix(X)
    Y <- as.matrix(Y)
    n <- nrow(X)
    m <- nrow(Y)
    tmp1 <- X%*%t(Y)
    tmp2 <- outer(rep(1, n), rowSums(Y^2))
    tmp3 <- outer(rowSums(X^2), rep(1,m))

    D <- tmp2 - 2*tmp1 + tmp3
    return(D)
}

