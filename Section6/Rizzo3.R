# Generating Random Variates
## from:  Statistical Computing with R, Maria Rizzo Chap.3

# inverse transform method 

rlogarithmic <- function (n, theta) {
    u <- runif(n)
    N <- ceiling(-16 / log10(theta))
    k <- 1:N
    a <- -1/log(1-theta)
    fk <- exp(log(a) + k * log(theta) - log(k))
    Fk <- cumsum(fk)
    x <- integer(n)
    for (i in 1:n) {
        x[i] <- as.integer(sum(u[i] > Fk))
        while (x[i]==N) {
            logf <- log(a) + (N+1) * log(theta) - log(N+1)
            fk <- c(fk, exp(logf))
            Fk <- c(Fk, Fk[N] + fk[N+1])
            N < N + 1
            x[i] <- as.integer(sum(u[i] > Fk))
        }
    }
    x + 1
}


n <- 1000
theta <- 0.5
x <- rlogarithmic(n, theta)

# get density
k <- sort(unique(x))
p <- -1 / log(1-theta) * theta^k / k
se <- sqrt(p*(1-p)/n)

tbl <- round(rbind(table(x)/n, p, se), 3)
tbl

## Acceptance-Rejection Method

n <- 1000
k <- 0
j <- 0
y <- numeric(n)

while (k < n) {
    u <- runif(1)
    j <- j + 1
    x <- runif(1)
    if (x * (1-x) > u) {
        k <- k + 1
        y[k] <- x
    }
}

j # number of iterations required to generate

p <- seq(.1, .9, .1)
Qhat <- quantile(y, p)
Q <- qbeta(p, 2, 2)
se <- sqrt(p * (1-p)/(n * dbeta(Q, 2, 2)))

tbl <- round(rbind(Qhat, Q, se), 3)


## MultiVariate - Choleski

rmnv.Choleski <- function(n, mu, Sigma) {
    # generate n random vectors from MVN(mu, Sigma)
    # dimension is inferred from mu and Sigma
    d <- length(mu)
    Q <- chol(Sigma) #Choleski factorization of Sigma
    Z <- matrix(rnorm(n * d), nrow=n, ncol=d)
    X <- Z %*% Q + matrix(mu, n, d, byrow=T)
    X
}

y <- subset(x=iris, Species=="virginica")[, 1:4]
mu <- colMeans(y)
Sigma <- cov(y)

#now generate MVN data with this mean and covariance
X <- rmnv.Choleski(200, mu, Sigma)
pairs(X)

## Symmetric Random Walk
n <- 400
incr <- sample(c(-1, 1), size = n, replace=T)
S <- as.integer(c(0, cumsum(incr)))
plot(0:n, S, type="l", main = "", xlab="i")