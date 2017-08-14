#'
#' omni embedding
#'
#' Get the omnibus embedding and compute the test statistic,
#' which is given by the squared Frobenius norm between the two embeddings.
#'
#' @param A a \eqn{n} x \eqn{n} matrix
#' @param B a \eqn{n} x \eqn{n} matrix
#' @param d an embedding dimension
#' @param nbs a number of bootstrapping
#' @param do.null a boolean to decide to calculate null distribution of the test statistics
#' @param runQ a boolean to decide to calculate null distribution of the test statistics for the given graphs
#' @param pvec an edge probability vector for ER random graphs for the null statistics
#'
#' @author Youngser Park <youngser@jhu.edu>
#' @export
#' @import igraph

omni <- function( A, B, d=2, nbs=500, do.null=TRUE, runQ=TRUE, pvec=0.5) {
    # Tobs
    if( (nrow(A) != ncol(A)) | (nrow(B) != ncol(B)) ) {
        stop('Matrices must be square.');
    }
    if( (nrow(A) != nrow(B)) | (ncol(A) != ncol(B)) ) {
        stop('Matrix dimensions disagree.');
    }
    nn <- nrow(A);

    # construct the omni matrix and compute test statistic.
    Tobs <- omni.computeT2(A,B,d);

    if (do.null) {

        if (runQ) {
            # Q1
            out <- ase( A, d );
            XhatA <- out$X
            Q1 <- XhatA %*% t(XhatA);

            # Q2
            out <- ase( B, d );
            XhatB <- out$X
            Q2 <- XhatB %*% t(XhatB);

            # Q12
            out <- ase( (A+B)/2, d );
            XhatAB <- out$X
            Q12 <- XhatAB %*% t(XhatAB);
        }

        # We're going to draw from this Q's as though it were the truth.
        BSsamps <- matrix(0,3+length(pvec),nbs);
        rownames(BSsamps) <- c("Q1","Q2","Q12",paste0("ER-p",pvec))
        for( i in 1:nbs ) {
            if (runQ) {
                Abs <- rg.sample(Q1);
                Bbs <- rg.sample(Q1);
                BSsamps[1,i] <- omni.computeT2(Abs,Bbs,d);
                Abs <- rg.sample(Q2);
                Bbs <- rg.sample(Q2);
                BSsamps[2,i] <- omni.computeT2(Abs,Bbs,d);
                Abs <- rg.sample(Q12);
                Bbs <- rg.sample(Q12);
                BSsamps[3,i] <- omni.computeT2(Abs,Bbs,d);
            }
            for (j in 1:length(pvec)) {
                Abs <- as.matrix(random.graph.game(nn,pvec[j])[])
                Bbs <- as.matrix(random.graph.game(nn,pvec[j])[])
                BSsamps[3+j,i] <- omni.computeT2(Abs,Bbs,d);
            }
        }
    } else {
        BSsamps <- NULL
    }

    return(list(Tobs=Tobs,Tbs=BSsamps))
}

#' @import irlba

ase <- function(A, dim){
    diag(A) <- rowSums(A) / (nrow(A)-1)

    if(nrow(A) >= 400){
        A.svd <- irlba::irlba(A, nu = dim, nv = dim)
        A.svd.values <- A.svd$d[1:dim]
        A.svd.vectors <- A.svd$v[,1:dim,drop=F]
        if(dim == 1)
            A.coords <- sqrt(A.svd.values) * A.svd.vectors
        else
            A.coords <- A.svd.vectors %*% diag(sqrt(A.svd.values))
    } else{
        A.svd <- svd(A)
        A.svd.values <- A.svd$d[1:dim]
        if(dim == 1)
            A.coords <- matrix(A.svd$v[,1] * sqrt(A.svd$d[1]),ncol=dim)
        else
            A.coords <- A.svd$v[,1:dim] %*% diag(sqrt(A.svd$d[1:dim]))
    }

    return(list(X=A.coords,D=A.svd.values))
}

ase2 <- function(A, dim){
    diag(A) <- rowSums(A) / (nrow(A)-1)

    if(nrow(A) >= 400){
        A.svd <- irlba::irlba(A, nu = dim, nv = dim)
        D <- A.svd$d[1:dim]
        X <- A.svd$u[,1:dim,drop=F]
        Y <- A.svd$v[,1:dim,drop=F]
        if(dim == 1) {
            X <- sqrt(D) * X
            Y <- sqrt(D) * Y
        }
        else {
            X <- X %*% diag(sqrt(D))
            Y <- Y %*% diag(sqrt(D))
        }
    } else{
        A.svd <- svd(A)
        D <- A.svd$d[1:dim]
        if(dim == 1) {
            X <- matrix(A.svd$u[,1,drop=F] * sqrt(D[1]),ncol=dim)
            Y <- matrix(A.svd$v[,1,drop=F] * sqrt(D[1]),ncol=dim)
        }
        else {
            X <- A.svd$u[,1:dim] %*% diag(sqrt(D[1:dim]))
            Y <- A.svd$v[,1:dim] %*% diag(sqrt(D[1:dim]))
        }
    }

    return(list(X=X,Y=Y,D=D))
}

ase3 <- function(A, dim, scaling=TRUE){
    diag(A) <- rowSums(A) / (nrow(A)-1)

    A.svd <- irlba::irlba(A, nu = dim, nv = dim)
    A.svd.values <- A.svd$d[1:dim]
    A.svd.vectors <- A.svd$v[,1:dim]

    if (scaling) {
        if(dim == 1)
            A.coords <- sqrt(A.svd.values) * A.svd.vectors
        else
            A.coords <- A.svd.vectors %*% diag(sqrt(A.svd.values))
    } else {
        A.coords <- A.svd.vectors
    }

    return(list(X=A.coords,D=A.svd.values))
}

rg.sample <- function(P){
    n <-  nrow(P)
    U <- matrix(0, nrow = n, ncol = n)
    U[col(U) > row(U)] <- stats::runif(n*(n-1)/2)
    U <- (U + t(U))
    A <- (U < P) + 0 ;
    diag(A) <- 0
    return(A)
}

omni.construct <- function( A11, A12, A21, A22 ) {
    # Take these four matrices and put them in an omnibus matrix,
    # with layout [ A11 A12 ]
    #		[ A21 A22 ]
    if( (nrow(A11) != nrow(A12))
        | (nrow(A21) != nrow(A22))
        | (ncol(A11) != ncol(A21))
        | (ncol(A12) != ncol(A22)) ) {
        stop('Matrix dimensions incompatible for forming omnibus matrix.');
    }
    nn <- nrow(A11);

    M <- matrix( nrow=2*nn, ncol=2*nn );
    M[1:nn, 1:nn] <- A11;
    M[(nn+1):(2*nn), (nn+1):(2*nn)] <- A22;
    M[1:nn, (nn+1):(2*nn)] <- A12;
    M[(nn+1):(2*nn), 1:nn] <- A21;

    return(M);
}

omni.computeT2 <- function( AA, BB, d ) {
    if( (nrow(AA) != ncol(AA)) | (nrow(BB) != ncol(BB)) ) {
        stop('Matrices must be square.');
    }
    if( (nrow(AA) != nrow(BB)) | (ncol(AA) != ncol(BB)) ) {
        stop('Dimension mismatch in matrices.');
    }
    nn <- nrow(AA);
    M <- omni.construct( AA, (AA+BB)/2, (AA+BB)/2, BB );
    # Get the omnibus embedding and compute the test statistic,
    # which is given by the squared Frobenius norm between the two embeddings.
    #  d <- ifelse(is.null(d), nrow(M), min(d,nrow(M)))
    out <- ase( M, d );
    #  dhat <- ifelse(d==2, d, max(2,getElbows(out$D,plot=FALSE)[2]))
    Zhat <- out$X#[,1:dhat]
    T <- norm(Zhat[1:nn,,drop=FALSE] - Zhat[(nn+1):(2*nn),,drop=FALSE], type="F");
    return(T^2)
}
