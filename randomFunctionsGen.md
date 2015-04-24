# Title


```r
require(knitr)
```


## Definitions

_Gaussian process_ is a collection of random variables, any
finite number of which have a joint Gaussian distribution.

A Gaussian process is completely specified by its mean function and
co- variance function. We define mean function \(m(x)\) and the
covariance function \(k(x,x')\) of a real process \(f(x)\) as
_covariance_ and _mean_ function
\[m(x) = \mathcal{E}[f(x)],\]
\[ k(x,x') = \mathcal{E}[(f(x) − m(x))(f(x') − m(x'))],\]
and will write the Gaussian process as
\( f(x) \sim GP(m(x), k(x, x'))\).

## Squared exponential

\[ cov(f(x_p), f (x_q)) = k(x_p , x_q) =
    \exp\Big(−\frac{1}{2} |x_p − x_q |\Big)^2.\]


```r
2
```

```
## [1] 2
```
