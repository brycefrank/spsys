<img src="docs/logo.png" width="60"></img><br> **sys** is an R package that allows users to explore systematic variance estimators.

Many environmental surveys use systematic sampling to produce estimates of population parameters. Estimating the precision of these quantities has proven a difficult task, with several systematic variance estimators proposed over the past several decades. While not exhaustive, `sys` implements several different variance estimators and provides diagnostic and simulation tools to allow analysts to select an appropriate variance estimator for their population.

## Standardization Format

This is a temporary idea for standard input for variance estimation using `sys` given a set of arbitrary (plot-level) attributes.

```
x    - x Coordinate
y    - y Coordinate
z_1  - First attribute
...
z_R  - Rth attribute
pi_i - The first order inclusion probability
```


## Variance Estimators

### Local Variance Estimators

- Aune-Lundberg & Strand (2014)
- Matern (1980)
- Stevens and Olsen (2003)

### Other

- Dorazio (2003)