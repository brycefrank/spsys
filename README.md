<img src="docs/logo.png" width="60"></img><br> **sys** is an R package that allows users to explore systematic variance estimators.

Many environmental surveys use systematic sampling to produce estimates of population parameters. Estimating the precision of estimates from systematic sampling designs has proven a difficult task, with several systematic variance estimators proposed over the past several decades. While not exhaustive, `sys` implements several different variance estimators and provides diagnostic and simulation tools to allow analysts to select an appropriate variance estimator for their population. More specifically, `sys` provides variance estimation for surveys that rely on point estimates of attributes of interest that use the Horvitz-Thompson estimator in two-dimensional settings.

`sys` supports both hexagonal and rectangular grids and has variance estimators implemented specifically for both configurations.

## Installation

Install this package directly from GitHub using `devtools`

```{r}
library(devtools)
devtools::install_github('https://github.com/brycefrank/sys)
```

## Getting Started

`sys` operates on a modified version of the now ubiquitious `sp` package `SpatialPointsDataFrame` class. The entry point into `sys` are the `HexFrame` and `RectFrame` classes, which represent hexagonal systematic and rectangular systematic sampling configurations respectively.

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
var_srs(hex_frame, N=N, fpc=TRUE)
```

Please see the vignette for formal descriptions of the variance estimators, their required arguments and other relavent details.

### Variance Assessment

`sys` provides two ways to assess the behavior of variance estimators: via synthetic populations and via subsampling.

One way to assess the performance of variance estimators is to treat an existing `SysFrame` as a population, and subsample repeatedly from it. `sys` enables subsampling using the following

```{r}
a <- 3
subsample(hex_frame, c(1,1), a)
```

where `c(1,1)` is the starting position in index space and `a=3` is the sampling interval. Fans of the `dplyr` package may want to use pipe operators. That is also possible

```{r}
hex_frame %>% subsample(c(1,1), a)
```

In most assessments we will be interested in all possible subsamples. Here we iterate over all possible subsamples and compute the simple random sampling with replacement estimator

```{r}
all_starts <- subsample_starts(a)

for(i in 1:nrow(all_starts)) {
  subsample(hex_frame, all_starts[i,]) %>%
    var_srs(fpc=FALSE)
}
```

## Development

This package is currently in development. Interested collaborators can email the author at bryce.frank@oregonstate.edu