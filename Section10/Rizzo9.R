# Markov Chain Monte Carlo Methods
## from:  Statistical Computing with R, Maria Rizzo Chap.9

# Metropolis-Hastings sampler (example with Rayleigh Distribution)
f <- function(x, sigma) {
    if (any(x<0)) return (0)
    stopifnot(sigma > 0)
    return((x / sigma^2) * exp(-x^2/ (2*sigma^2)))
}

m <- 10000
sigma <- 4
x <- numeric(m)
x[1] <- rchisq(1, df=1)
k <- 0
u <- runif(m)


for (i in 2:m) {
    xt <- x[i-1]
    y <- rchisq(1, df = xt)
    num <- f(y, sigma) * dchisq(xt, df=y)
    den <- f(xt, sigma) * dchisq(y, df=xt)
    if (u[i] <= num/den) x[i] <- y else {
        x[i] <- xt
        k <- k + 1
        }
}

index <- 5000:5100
y1 <- x[index]
plot(index, y1, type="l", main="", ylab="x")

# compare MCMC generated with Raleigh(sigma=4) / intended distr
b <- 2001
y <- x[b:m]
a <- ppoints(100)
QR <- sigma * sqrt(-2 * log(1-a))
Q <- quantile(x, a)

qqplot(QR, Q)
hist(y, breaks="scott", freq=F)
lines(QR, f(QR,4))

# Random Walk Metropolis
## t distribution

rw.Metropolis <- function(n, sigma, x0, N) {
    x <- numeric(N)
    x[1] <- x0
    u <- runif(N)
    k <- 0
    for (i in 2:N) {
        y <- rnorm(1, x[i-1], sigma)
            if (u[i] <= (dt(y,n) / dt(x[i-1], n)))
                x[i] <- y else {
                    x[i] <- x[i-1]
                    k <- k+1
                }
    }
    return(list(x=x, k=k))
}

n <- 4
N <- 2000
sigma <- c(.05, .5, 2, 16)

x0 <- 25
rw1 <- rw.Metropolis(n, sigma[1], x0, N)
rw2 <- rw.Metropolis(n, sigma[2], x0, N)
rw3 <- rw.Metropolis(n, sigma[3], x0, N)
rw4 <- rw.Metropolis(n, sigma[4], x0, N)

#num rejected points
print(c(rw1$k, rw2$k, rw3$k, rw4$k))

# plot
par(mfrow=c(2,2))
refline <- qt(c(.025, .975), df=n)
rw <- cbind(rw1$x, rw2$x, rw3$x, rw4$x)
for (j in 1:4) {
    plot(rw[, j], type='l', 
        xlab=bquote(sigma == .(round(sigma[j], 3))),
        ylab = "X", ylim=range(rw[,j]))
    abline(h=refline)
}
par(mfrow=c(1,1))

a <- c(.05, seq(.1, .9, .1), .95)
Q <- qt(a, n)
rw <- cbind(rw1$x, rw2$x, rw3$x, rw4$x)
mc <- rw[501:N, ]
Qrw <- apply(mc, 2, function(x) quantile(x, a))
print(round(cbind(Q, Qrw), 3))

## simple investment model

b <- .2
w <- .25
m <- 5000
burn <- 1000
days <- 250
x <- numeric(m)
i <- sample(1:5, size=days, replace=TRUE, 
            prob=c(1, 1-b, 1-2*b, 2*b, b))
win <- tabulate(i)

prob <- function(y, win) {
    if (y < 0 || y >= 0.5)
        return(0)
    return((1/3)^win[1] * 
               ((1-y)/3)^win[2] * ((1-2*y)/3)^win[3] * 
               ((2*y)/3)^win[4] * (y/3)^win[5])
}

u <- runif(m)
v <- runif(m, -w, w)
x[1] <-.25
for (i in 2:m) {
    y <- x[i-1] + v[i]
    if (u[i] <= prob(y, win) / prob(x[i-1], win))
        x[i] <- y else
            x[i] <- x[i-1]
}

print(win)
print(round(win/days, 3))
print(round(c(1, 1-b, 1-2*b, 2*b, b)/3, 3))
xb <- x[(burn+1):m]
print(mean(xb))

# Independance Sampler

M <- 5000
xt <- numeric(m)
a <- 1
b <- 1
p <- .2
n <- 30
mu <- c(0, 5)
sigma <- c(1, 1)

i <- sample(1:2, size=n, replace=T, prob=c(p, 1-p))
x <- rnorm(n, mu[i], sigma[i])

u <- runif(m)
y <- rbeta(m, a, b)
xt[1] <- .5

for (i in 2:m) {
    fy <- y[i] * dnorm(x, mu[1], sigma[1]) +
        (1-y[i]) * dnorm(x , mu[2], sigma[2])
    fx <- xt[i-1] * dnorm(x, mu[1], sigma[1]) +
        (1-xt[i-1]) * dnorm(x , mu[2], sigma[2])
    
    r <- prod(fy / fx) * 
        (xt[i-1] ^ (a-1) * (1-xt[i-1]) ^(b-1)) /
            (y[i]^(a-1) * (1-y[i]) ^(b-1))
    
    if (u[i] <= r) xt[i] <- y[i] else
        xt[i] <- xt[i-1]
        
}

plot(xt, type="l", ylab="p")
hist(xt[101:m], main="", xlab="p", prob=T)
print(mean(xt[101:m]))

## Gibbs Sampler

N <- 5000
burn <- 1000
X <- matrix(0, N, 2)

rho <- -.75
mu1 <- 0
mu2 <- 2
sigma1 <- 1
sigma2 <- .5
s1 <- sqrt(1-rho^2)*sigma1
s2 <- sqrt(1-rho^2)*sigma2

X[1, ] <- c(mu1, mu2)

for (i in 2:N) {
    x2 <- X[i-1, 2]
    m1 <- mu1 + rho * (x2-mu2) * sigma1/sigma2
    X[i, 1] <- rnorm(1, m1, s1)
    x1 <- X[i, 1]
    m2 <- mu2 + rho * (x1 - mu1) * sigma2/sigma1
    X[i, 2] <- rnorm(1, m2, s2)
}

b <- burn + 1
x <- X[b:N,]
colMeans(x)
cov(x)
cor(x)

plot(x, cex=.5, xlab=bquote(X[1]), ylab=bquote(X[2]), ylim=range(x[,2]))

## Gelman Rubin method for monitoring MCMC Convergence

# in example target distribution = normal(0,1)
# proposal distribution is Normal(X_t, sigma^2)

Gelman.Rubin <- function(psi){
    ps <- as.matrix(psi)
    n <- ncol(psi)
    k <- nrow(psi)
    
    psi.means <- rowMeans(psi)
    B <- n * var(psi.means)
    psi.w <- apply(psi, 1, "var")
    W <- mean(psi.w)
    v.hat <- W*(n-1)/n + (B/n)
    r.hat <- v.hat / W
    return(r.hat)
}

normal.chain <- function(sigma, N, X1) {
    x <- rep(0, N)
    x[i] <- X1
    u <- runif(N)
    
    for (i in 2:N) {
        xt <- x[i-1]
        y <- rnorm(1, xt, sigma)
        r1 <- dnorm(y, 0, 1)
        r2 <- dnorm(xt, 0, 1)
        r <- r1 / r2
        if (u[i] <= r) x[i] <- y else
            x[i] <- xt
    }
    return(x)
}

sigma <- .2
k <- 4
n <- 15000
b <- 1000

# choose overdispersed initial values
x0 <- c(-10, -5, 5, 10)

X <- matrix(0, nrow=k, ncol=n)
for (i in 1:k)
    X[i, ] <- normal.chain(sigma, n, x0[i])

psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
    psi[i, ] <- psi[i,] / (1:ncol(psi))
print(Gelman.Rubin(psi))

par(mfrow=c(2,2))
for (i in 1:k) {
    plot(psi[i, (b+1):n], type="l",
         xlab=i, ylab=bquote(psi))}
par(mfrow=c(1,1))

rhat <- rep(0, n)
for(j in (b+1):n){
    rhat[j] <- Gelman.Rubin(psi[, 1:j])
}
plot(rhat[(b+1):n], type="l", xlab="", ylab="R")
abline(h=1.1, lty=2)

