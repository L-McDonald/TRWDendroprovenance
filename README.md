# Tree-Ring Width Dendroprovenancing
Project using oak ring-width patterns from across Europe to determine the
provenance of unknown samples.

## Instructions to reproduce

## Scripts
### 00_Functions.R
Various functions used throughout the project

### 01_SpeciesDistributionRanges.R
Plots the species distribution ranges for Quercus robur, Q. petraea, Q. cerris, 
Q. frainetto, and Q. pubescens from
[EuForGen](https://data.mendeley.com/datasets/hr5h2hcgg4/19). Determines the
maximum extent of species occurrence for downloading from the ITRDB.

### 02_DownloadChronologiesITRDB.R
Downloads all ITRDB chronologies of the above species and saves RWL files along
with metadata.

### 03_PrepareCTRL.R
Reads the rwl files for the CTRL chronologies and saves as a nested tibble with
metadata.

### 04_PrepareWazny.R
Reads the rwl files for the Wazny chronologies and saves as a nested tibble with
metadata.

### 05_RemoveHistorical.R
Removes historical sequences from the ITRDB chronologies.

### 06_MergeDatasets.R
Merges the ITRDB, CTRL, and Wazny chronologies into a single nested tibble. Also
adds columns for coordinates in EPSG 3035.

### 07_CleanChronologies.R
Cleans the merged dataset in four steps: (1) checks for and removes chronologies 
with coordinates located at sea, (2) combines chronologies located within 1 km of 
each other, (3) checks for and removes problem series/chronologies according to 
the method of van Sluijs et al., (4) removes sites with chronologies shorter than 
40 years and also removes some specific series which could not be detrended using
the spline approach due to negative values. Tests showed that removing these series
had negligible effect on the indices (calculated in script 08).

### 08_CalculateChronologies.R
Detrends the rwl measurements using a smoothing spline. Then calculates three
mean chronologies: std, res, and ars. ars is Cook's ARSTAN chronology. The ring
width indices for six sites (DEU304, DEU306, DEU312, DEU314, SWIT232, and SWIT317)
have an AR model order of zero which prevents ars being calculated. For these an 
AR model of order 1 was fitted which tests showed resulted in a time series very
closely correlated with the std chronology.

### 09_CalculateInterseriesCorrelations.R
Calculates pairwise Pearson correlations (with the number of overlapping years)
between all chronologies for each of the std, res, and ars series, and pairwise
Euclidean distances (m) between chronology locations using EPSG 3035 coordinates.
Combines these into a single long-format table with one row per unique site
pair.
