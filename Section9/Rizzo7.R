# Bootstrap and Jackknife
## from:  Statistical Computing with R, Maria Rizzo Chap.7

# bootstrap method to estimate SE
library(bootstrap)
cor(law$LSAT, law$GPA)
cor(law82$LSAT, law82$GPA)

B <- 200
n <- nrow(law)
R <- numeric(B)

for (b in 1:B) {
    #randomly select indices
    i <- sample(1:n, size=n, replace=T)
    LSAT <- law$LSAT[i]
    GPA <- law$GPA[i]
    R[b] <- cor(LSAT, GPA)
}

se.R <- sd(R)
hist(R, prob=T)

#using boot function
library(boot)
r <- function(x, i) {
    #get cor of columsn 1 and 2
    cor(x[i, 1], x[i, 2])
}

obj <- boot(data=law, statistic=r, R=2000)
obj

## Jackknife
# like leave one out cross validation

# jackknife estimate of bias
data(patch, package='bootstrap')
n <- nrow(patch)
y <- patch$y
z <- patch$z
theta.hat <- mean(y) / mean(z)

theta.jack <- numeric(n)
for (i in 1:n) {
    theta.jack[i] <- mean(y[-i]) / mean(z[-i])
}
bias <- (n-1) * (mean(theta.jack) - theta.hat)

# jacknife failure - when population mean not smooth
# smoothness: small changes in the data = small changes in static 
# median, potentially not smooth

n <- 10
x <- sample(1:100, size=n)

#jackknife estimate of se
M <- numeric(n)
for (i in 1:n) {
    y <- x[-i]
    M[i] <- median(y)
}
Mbar <- mean(M)
print(sqrt((n-1)/n * sum((M-Mbar)^2)))

#bootstrap estimate of se
Mb <- replicate(1000, expr = {
    y <- sample(x, size=n, replace=T)
    median(y)
})
print(sd(Mb))

## jacknife after bootstrap

data(patch, package="bootstrap")
n <- nrow(patch)
y <- patch$y
z <- patch$z
B <- 2000
theta.b <- numeric(B)
indices <- matrix(0, nrow=B, ncol=n)

#step 1: run bootstrap
for (b in 1:B) {
    i <- sample(1:n, size=n, replace=T)
    y <- patch$y[i]
    z <- patch$z[i]
    theta.b[b] <- mean(y) / mean(z)
    indices[b, ] <- i
}

#step 2 jackknife after boostrap to est se
se.jack <- numeric(n)
for (i in 1:n) {
    keep <- (1:B)[apply(indices, MARGIN = 1, 
                        FUN = function(k) {!any(k==i)})]
    se.jack[i] <- sd(theta.b[keep])
}

print(sd(theta.b))
print(sqrt((n-1) * mean((se.jack - mean(se.jack))^2)))