# Monte Carlo Methods in Inference
## from:  Statistical Computing with R, Maria Rizzo Chap.5

## estimate mean and standard error
m <- 1000
g <- numeric(m)
for (i in 1:m) {
    x <- rnorm(2)
    g[i] <- abs(x[1] - x[2])
}
est <- mean(g)

## pop mean and (n-1) doesn't matter with large sample like above
se <- sqrt(sum((g-mean(g))^2)) / m
se

# trimmed mean 
n <- 20 
m <- 1000
tmean <- numeric(m)
for (i in 1:m) {
    x <- sort(rnorm(n))
    tmean[i] <- sum(x[2: (n-1)]) / (n-2)
}
mse <- mean(tmean^2)
se <- sqrt(sum((tmean - mean(tmean)^2))) / m

n <- 20
K <- n/2 - 1
m <- 1000
mse <- matrix(0, n/2, 6)
trimmed.mse <- function(n, m, k, p) {
    #MC est of mse for k-level trimmed mean of 
    # contaminated normal pN(0,1) + (1-p)N(0,100)
    tmean <- numeric(m)
    for (i in 1:m){
        sigma <- sample(c(1,10), size=n,
                        replace=T, prob=c(p, 1-p))
        x <- sort(rnorm(n, 0 , sigma))
        tmean[i] <- sum(x[(k+1):(n-k)]) / (n-2*k)
    }
    mse.est <- mean(tmean^2)
    se.mse <- sqrt(mean((tmean-mean(tmean))^2)) / sqrt(m)
    return(c(mse.est, se.mse))
}

for (k in 0:K) {
    mse[k+1, 1:2] <- trimmed.mse(n=n, m=m, k=k, p=1)
    mse[k+1, 3:4] <- trimmed.mse(n=n, m=m, k=k, p=.95)
    mse[k+1, 5:6] <- trimmed.mse(n=n, m=m, k=k, p=.9)
}

# MC experiment to estimate a confidence level
n <- 20
alpha <- .05
UCL <- replicate(1000, expr = {
    x <- rnorm(n, mean=0, sd=2)
    (n-1) * var(x) / qchisq(alpha, df=n-1)
})
sum(UCL > 4)
mean(UCL > 4) # 95% of the UCLs are > 4