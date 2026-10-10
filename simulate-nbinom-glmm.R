simulate_nbinom_glmm <- function(
    n_clusters,
    cluster_size,
    beta,
    Gamma,
    theta,
    x = NULL
) {
    stopifnot(
        length(beta) == 2L,
        is.matrix(Gamma),
        identical(dim(Gamma), c(2L, 2L)),
        theta > 0,
        n_clusters >= 1,
        cluster_size >= 1
    )

    n <- n_clusters * cluster_size

    # Cluster identifiers
    id <- rep(
        seq_len(n_clusters),
        each = cluster_size
    )

    # Covariate
    if (is.null(x)) {
        x <- rep(
            seq(0, 1, length.out = cluster_size),
            times = n_clusters
        )
    }

    stopifnot(length(x) == n)

    # Generate random effects
    # Works for positive semidefinite Gamma,
    # including Gamma = 0.
    eig <- eigen(
        (Gamma + t(Gamma)) / 2,
        symmetric = TRUE
    )

    if (min(eig$values) < -1e-10) {
        stop("Gamma must be positive semidefinite.")
    }

    A <- eig$vectors %*%
        diag(
            sqrt(pmax(eig$values, 0)),
            nrow = 2
        )

    B <- matrix(
        rnorm(n_clusters * 2L),
        nrow = n_clusters,
        ncol = 2L
    ) %*% t(A)

    # Linear predictor
    eta <- beta[1] +
        beta[2] * x +
        B[id, 1] +
        B[id, 2] * x

    # Conditional mean
    mu <- exp(eta)

    # Negative binomial responses
    y <- rnbinom(
        n = n,
        size = theta,
        mu = mu
    )

    data.frame(
        id = factor(id),
        x = x,
        y = y
    )
}
