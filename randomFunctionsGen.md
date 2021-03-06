# Gaussians Animals


```r
require(knitr)
require(MASS)

require(data.table)
require(reshape2)

require(ggplot2)
```


## Definitions

_Gaussian process_ is a stochastic process where any finite number of
random variables have a joint Gaussian distribution.

A Gaussian process is completely specified by its _mean function_
\(m(x)\) and _covariance function_ \(k(x,x')\). They are defined by

\[m(x) = \mathcal{E}[f(x)],\]

\[ k(x,x') = \mathcal{E}[(f(x) − m(x))(f(x') − m(x'))],\]

and will write the Gaussian process as \( f(x) \sim GP(m(x), k(x,
x'))\).



## Squared exponential


\[ cov(f(x_p), f (x_q)) = k(x_p , x_q) =
    \exp\Big(−\frac{1}{2} |x_p − x_q |\Big)^2.\]


```r
squared.exp <- function(x,y,l=1) {
  exp(-0.5* (norm(x-y, type="2")/l)^2)
}

Sigma <- function(x1, x2, cov.fun=squared.exp) {
  n1 <- length(x1)
  n2 <- length(x2)
  S <- matrix(rep(0, n1*n2), nrow=n1)
  for (i in 1:n1) {
    for (j in 1:n2) {
      S[i,j] <- cov.fun(x1[i],x2[j]) }}
  S
}

melted.simulations <- function(x.star, sim.number, ms, S) {
  sims <- cbind(data.table(x=x.star), t(mvrnorm(sim.number, ms, S)))
  melt(sims, id=1)
}

gg.simulations <- function(melted.sims) {
  ggplot() +
      geom_line(data = melted.sims, aes(x=x,y=value,
                    group=factor(variable)), alpha=.6) +
          theme_bw()
}


x.star <- seq(-5,5,0.1)
sim.number <- 5
S <- Sigma(x.star,x.star)
msims <- melted.simulations(x.star, sim.number,
                            rep(0, length(x.star)), S)
gg.simulations(msims)
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2-1.png) 

## Constant

\[cov(f(x_p), f (x_q)) = k(x_p , x_q) = \sigma_0^2 \]


```r
sigma0 <- 0.1
cov.fun <- function(x1, x2) {sigma0^2}
S <- Sigma(x.star,x.star, cov.fun)
msims <- melted.simulations(x.star, sim.number,
                            rep(0, length(x.star)), S)
gg.simulations(msims)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png) 

## Exponential

\[cov(f(x_p), f (x_q)) = k(x_p , x_q) = \exp(-|x_p-x_q|/l) \]



```r
x.star <- seq(-5,5,0.1)
sim.number <- 5

exp.cov <- function(x,y,l=1) {
  exp(-norm(x-y, type="2")/l)
}


S <- abs(Sigma(x.star, x.star, exp.cov))
msims <- melted.simulations(x.star, sim.number,
                            rep(0, length(x.star)), S)
gg.simulations(msims)
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png) 

## Linear

\[cov(f(x_p), f (x_q)) = k(x_p , x_q) = \sigma_1\cdot x_p\cdot x_q \]


```r
x.star <- seq(-5,5,0.1)
sim.number <- 5

lin.cov <- function(x,y,sigma1=1) {
  (sigma1^2)*x*y
}


S <- abs(Sigma(x.star, x.star, lin.cov))
msims <- melted.simulations(x.star, sim.number,
                            rep(0, length(x.star)), S)
gg.simulations(msims)
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png) 

## Rational Quadratic

\[ cov(f(x_p), f (x_q)) = k(x_p , x_q) =
    \Big(1+\frac{r^2}{2\alpha l^2}\Big)^{-\alpha} \]



```r
alpha <- 2
l <- 1
rq.cov <- function(x,y,sigma1=1) {
  (1+(norm(x-y, type="2")^2/(2*alpha*l^2)))^(-alpha)
}

sim.number <- 9
S <- abs(Sigma(x.star, x.star, rq.cov))
msims <- melted.simulations(x.star, sim.number,
                            rep(0, length(x.star)), S)
ggplot() +
    geom_line(data = msims, aes(x=x,y=value,
                  colour=factor(variable)), alpha=.6) +
        theme_bw() + scale_color_brewer(palette = "Blues") +
            theme(legend.position="none")
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-1.png) 

## Brownian Motion


```r
brow <- function(x1,x2) {
  min(x1,x2)
}

x.star <- seq(0, 5, 0.005)
S <- Sigma(x.star, x.star, brow)
ms <- rep(0, length(x.star))
msims <- melted.simulations(x.star, sim.number, ms, S)
gg.simulations(msims)
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png) 

```r
bm <- data.table(t(mvrnorm(2, rep(0, length(x.star)), S)))
ggplot(bm, aes(V1,V2))+
    geom_path() + theme_bw()
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-2.png) 

## Predictions

\[\overline{y}^* = K^*(\lambda^2\cdot I+K)^{-1}\cdot y \]

\[C = K^{**}-K^*\cdot(\lambda^2\cdot I+K)^{-1}\cdot(K^*)^t \]


```r
y.predict <- function(x, y, x.star, cov.fun, sigma=0) {
  K <- Sigma(x,x)
  I <- diag(length(x))
  K.star <- Sigma(x.star, x)
  K.star %*% solve(K+sigma^2*I) %*% y
}

cov.predict <- function(x, y, x.star, cov.fun, sigma=0) {
  K <- Sigma(x,x)
  I <- diag(length(x))
  K.star <- Sigma(x.star, x)
  K.star.star <- Sigma(x.star, x.star)
  K.star.star - K.star %*% solve(K+sigma^2*I) %*% t(K.star)
}

gg.simulations2 <- function(x, y, x.star, cov.fun, sim.number, sigma=0) {
  observations <- data.table(x, y)
  y.star.bar <- y.predict(x, y, x.star, cov.fun, sigma)[,1]
  predictors <- data.table(x=x.star, y=y.star.bar)
  S <- cov.predict(x, y, x.star, cov.fun, sigma)
  ms <- y.star.bar
  msims <- melted.simulations(x.star, sim.number, ms, S)
  d <- abs(diag(S))
  predictors$sigma <- sqrt(d)
  
  ggplot() +
      geom_line(data = msims, aes(x=x, y=value,
                    group=factor(variable)), alpha=.3, colour="green") +
          geom_point(data=observations,
                     aes(x,y), size=3, colour="red")+
              geom_line(data=predictors, aes(x, y), colour='blue') +
                  geom_ribbon(data=predictors,
                              aes(x, ymax=y+2*sigma,
                                  ymin=y-2*sigma), alpha=0.1)+
                      theme_bw()
}                 

x <- c(-4,-3,-1,0,2)
y <- c(-2,0, 1, 2, -1)
x.star <- seq(-5,5,0.1)
sim.number <- 20
cov.fun <- squared.exp;

gg.simulations2(x, y, x.star, cov.fun, sim.number)
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png) 

## Brownian bridge


```r
x <- c(0,1)
y <- c(0,0)
x.star <- seq(0,1,0.01)
sim.number <- 20
cov.fun <- brow;
gg.simulations2(x, y, x.star, cov.fun, sim.number)
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-1.png) 

## Sigma non zero


```r
x <- c(-4,-3,-1,0,2)
y <- c(-2,0, 1, 2, -1)
x.star <- seq(-5,5,0.1)
sim.number <- 20
cov.fun <- squared.exp;
sigma <- 0.1

gg.simulations2(x, y, x.star, cov.fun, sim.number, sigma) +
    geom_errorbar(data=data.table(x,y),
                  aes(x,y, ymin=y-2*sigma, ymax=y+2*sigma), width=0.2)
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10-1.png) 

## Large n

If $n$ is large, for example $n\geq 10000$, then to solve the equation
is not practical. Therefore is is usefull to develope approximate
solution.
