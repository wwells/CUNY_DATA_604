# Monte Carlo Integration and Variance Reduction
## from:  Statistical Computing with R, Maria Rizzo Chap.5

#simple integration of theta = integral_0^1 e^(-1) dx
m <- 10000
x <- runif(m)
theta.hat <- mean(exp(-x))

theta.hat
1 - exp(-1)

#simple integration of theta = integral_2^4 e^(-1) dx
x <- runif(m, min=2, max=4)
theta.hat <- mean(exp(-x)) * 2

theta.hat
exp(-2) - exp(-4)

#monte carlo integration of standard normal cdf
x <- seq(.1, 2.5, length=10)
u <- runif(m)
cdf <- numeric(length(x))
for (i in 1:length(x)) {
    g <- x[i] * exp(-(u * x[i])^2 / 2)
    cdf[i] <- mean(g) / sqrt(2 * pi) + 0.5
}

Phi <- pnorm(x)
round(rbind(x, cdf, Phi), 3)

# control variance in MC sim for integral_0^1 e^(-x)/(1+x^2) dx
f <- function(u)
    exp(-.5)/(1+u^2)
g <- function(u)
    exp(-u)/(1+u^2)

set.seed(510)
u <- runif(10000)
B <- f(u)
A <- g(u)
cor(A, B) # initial variance
a <- -cov(A, B) / var(B)
a

T1 <- g(u)
T2 <- T1 + a * (f(u) - exp(-.5)* pi/4)
c(mean(T1), mean(T2))
(var(T1) - var(T2)) / var(T1)

# importance sampling, determine optimal function
m <- 10000
theta.hat <- se <- numeric(5)
g <- function(x) {
    exp(-x - log(1+x^2)) * (x > 0) * (x < 1)
}

x <- runif(m)
fg <- g(x)
theta.hat[1] <- mean(fg)
se[1] <- sd(fg)

x <- rexp(m, 1)
fg <- g(x) / exp(-x)
theta.hat[2] <- mean(fg)
se[2] <- sd(fg)

x <- rcauchy(m)
i <- c(which(x > 1), which(x < 0 ))
x[i] <- 2
fg <- g(x) / dcauchy(x)
theta.hat[3] <- mean(fg)
se[3] <- sd(fg)

u <- runif(m)
x <- -log(1 - u * (1 - exp(-1)))
fg <- g(x) / (exp(-x) / (1 - exp(-1)))
theta.hat[4] <- mean(fg)
se[4] <- sd(fg)

u <- runif(m)
x <- tan(pi * u / 4)
fg <- g(x) / (4 / ((1 + x^2) * pi))
theta.hat[5] <- mean(fg)
se[5] <- sd(fg)

rbind(theta.hat, se)