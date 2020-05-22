<img src="docs/logo.png" width="60"></img><br> **sys** is an R package that allows users to explore systematic variance estimators.

Many environmental surveys use systematic sampling to produce estimates of population parameters. Estimating the precision of these quantities has proven a difficult task, with several systematic variance estimators proposed over the past several decades. While not exhaustive, `sys` implements several different variance estimators and provides diagnostic and simulation tools to allow analysts to select an appropriate variance estimator for their population.

## Getting Started

`sys` operates on the `sp` package `SpatialPointsDataFrame` class. 

The first step is to load in your `SpatialPointsDataFrame` that represents your systematic sample using the proper input function. The input function will standardize your dataset to a `SysFrame` class, upon which further analysis is conducted.

For example, we can load in a set of points from a hexagonal grid:

```{r}

hex_points <- readOGR('my_hex_points.shp')
hex_frame <- load_hex(hex_points, c('vol', 'ba'))
```

The input function `load_hex` takes two arguments. A `SpatialPointsDataFrame` and a vector of column names that indicate attributes we are interested in conducting the analyeses on. Here we indicate volume (vol) and basal area (ba) as our attributes of interest. `load_hex`, and its sister function `load_rect` (for rectangular systematic samples), achieve several things. First, the points are cast onto a set of polygons representing the tesselation of the study area. Second, they implement a standardized indexing system that is used in several variance estimation functions. Finally, they provide a standardized interface for different types of analysis functions.

The `hex_frame` variable now represents an object of the `SysFrame` class.

### Subsampling

One way to assess the performance of variance estimators is to treat an existing `SysFrame` as a population, and subsample repeatedly from it.

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