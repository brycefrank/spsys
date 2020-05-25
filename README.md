<img src="docs/logo.png" width="60"></img><br> **sys** is an R package that allows users to explore systematic variance estimators.

Many environmental surveys use systematic sampling to produce estimates of population parameters. Estimating the precision of these quantities has proven a difficult task, with several systematic variance estimators proposed over the past several decades. While not exhaustive, `sys` implements several different variance estimators and provides diagnostic and simulation tools to allow analysts to select an appropriate variance estimator for their population.

## Getting Started

`sys` operates on the `sp` package `SpatialPointsDataFrame` class. 

The first step is to load in your `SpatialPointsDataFrame` that represents your systematic sample using the proper input function. The input function will standardize your dataset to a `SysFrame` class, upon which further analysis is conducted.

For example, we can load in a set of points from a hexagonal grid:

```{r}

hex_points <- readOGR('my_hex_points.shp')
hex_frame <- HexFrame(hex_points, c('vol', 'ba'))
```

`HexFrame` takes two arguments: a `SpatialPointsDataFrame` and a vector of column names that indicate attributes we are interested in conducting the analyses on. Here we indicate volume (vol) and basal area (ba) as our attributes of interest. `HexFrame`, and its sister class `RectFrame` (for rectangular systematic samples), achieve several things. First, the points are cast onto a set of polygons representing the tesselation of the study area. Second, they implement a standardized indexing system that is used in several variance estimation functions. Finally, they provide a standardized interface for different types of analysis functions.

The `hex_frame` variable now represents an object of the `HexFrame` class.

### Estimating Variances

Given a `HexFrame` or `RectFrame` that represents a sample of the study area, it is elementary to estimate variances using a host of available estimators. For starters we can estimate the variance assuming simple random sampling without replacement. Assume the population size is equal to 10,000. We obtain

```
N <- 10000
var_srs(hex_frame, N=N)
```

### Subsampling

One way to assess the performance of variance estimators is to treat an existing `SysFrame` as a population, and subsample repeatedly from it. `sys` enables subsampling using the following

```{r}
subsample(hex_frame, c(1,1), 3)
```

where `c(1,1)` is the starting position in index space and `3` is the sampling interval. Fans of the `dplyr` package may want to use pipe operators. That is also possible

```{r}
hex_frame %>% subsample(c(1,1), 3)
```

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