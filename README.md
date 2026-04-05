
# MFRC

<!-- badges: start -->
<!-- badges: end -->

**MFRC** is the R package for mineral formula recalculation. This includes:

- Basic functions to convert wt% into atom% (`Mass2Cation()`) and atom% into wt% (`Cation2MassPc()`)
- Site fraction calculators for each mineral

## Installation

You can install this package from [GitHub](https://github.com/) with:

``` r
install.packages("remotes")
remotes::install_github("Ryo-fkushima/MFRC")
```

## How to use

1. Prepare a dataframe of wt%- or atom%-based compositional table.
2. If the input is wt%, use `Mass2Cation()` to convert it into atomic proportions.
3. Apply a site-fraction calculator.


## Example

``` r
install.packages("clipr")

library("MFRC")
library("clipr")

# Copy pyroxene wt% table from somewhere (e.g., Excel sheet) and run the command below.
# The data in the clipboard will be converted into the site fraction table.
read_clip_tbl() |> Mass2Cation(MFRC::ElementList_MnNCKFMASTCr) |> Pyroxene_MinPlot() |> write_clip()

# Paste the result somewhere to check it.

```

## Implemented site-fraction calculators (as of v0.1.0)

- `Omp_AlT_Omit()`: Omphacite. Al(T) was ignored and Fe3+ = Na - Al.
- `Pyroxene_MinPlot()`: Pyroxene. Equivalent to the method used in MinPlot (Walters, 2022, Mineralogia)


## Input table format
 
Examples of the input table are presented below. 

wt% table:

|SpotNo|Phase|Na2O|SiO2|FeO|Al2O3|MgO|CaO|
:---|:---|:---|:---|:---|:---|:---|:---|
|p1|Omp|7.20|56.09|5.49|10.78|7.31|13.04|
  
atom% table:

|SpotNo|Phase|Si|Al|Fe2|Mn2|Mg|Ca|Na|
:---|:---|:---|:---|:---|:---|:---|:---|:---|
|p1|Omp|8.03|1.82|0.66|0.00|1.56|2.00|2.00|

## Note

- The input table data does not have to be normalized to 100%.
- One does not have to remove non-concentration columns beforehand.
- Note that `FeO` should be used in the wt% input table, and `Fe2`/`Mn2` should be used for Fe/Mn in the atom% input table!
- `Fe3` is implemented in the package, but it is mainly for the `Cation2MassPc()` function.



## Author

Ryo Fukushima (rpifukushima@gmail.com)

## Updates

Apr 05, 2026: v0.1.0 was released.