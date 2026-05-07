
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
# The data in the clipboard will be converted into the site fraction table, and saved again in the clipboard (with overwriting)
read_clip_tbl() |> Mass2Cation(MFRC::ElementList_MnNCKFMASTCr) |> Pyroxene_MinPlot() |> write_clip()

# Paste the clipboard data somewhere to check results.

```

## Implemented site-fraction calculators (as of v0.2.0)

- `Omp_AlT_Omit()`: Omphacite. Al(T) is ignored and Fe3+ = Na - Al (e.g., Fukushima et al., 2024, Phys. Earth Planet. Inter.)
- `Pyroxene_MinPlot()`: Pyroxene. Equivalent to the method used in MinPlot (Walters, 2022, Mineralogia)


## Input table format
 
Examples of the input table are presented below. 

wt% table:

|SpotNo|Phase|Na2O|SiO2|FeO|Al2O3|MgO|CaO|
:---|:---|:---|:---|:---|:---|:---|:---|
|p1|Omp|7.20|56.09|5.49|10.78|7.31|13.04|

(other available elements: TiO2, Cr2O3, Y2O3, NiO, ZnO, MnO, K2O, BaO, H2O)
  
atom% table:

|SpotNo|Phase|Si|Al|Fe|Mn|Mg|Ca|Na|
:---|:---|:---|:---|:---|:---|:---|:---|:---|
|p1|Omp|8.03|1.82|0.66|0.00|1.56|2.00|2.00|

(other available elements: Ti, Cr, Y, Ni, Zn, K, Ba, H)

## Note

- The input table data does not have to be normalized to 100%.
- One does not have to remove non-concentration columns beforehand.
- Note that `FeO` should be used in the wt% input table.



## Author

Ryo Fukushima (rpifukushima@gmail.com)

## Updates
May 07, 2026 (v0.2.0): `Fe2O3`/`Fe3` option was removed for simplicity. `H2O` was implemented.

Apr 05, 2026: v0.1.0 was released.