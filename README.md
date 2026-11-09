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
maximum extent of species occurrence for downloading from the ITRDB. Merges the
five SDRs and saves in EPSG 4326 and 3035.

### 02_PrepareEnvironmentalVariables.R
Downloads coastline from natural earth, crops to approximate region, converts to
EPSG 3035 and saves. Takes environmental rasters from CRU, SoilGrids, and 
elevation from worldclim. Crops them to the study region with a buffer, converts
to EPSG 3035 and saves.

### 03_DownloadChronologiesITRDB.R
Downloads all ITRDB chronologies of the above species and saves RWL files along
with metadata.

### 04_PrepareCTRL.R
Reads the rwl files for the CTRL chronologies and saves as a nested tibble with
metadata.

### 05_PrepareWazny.R
Reads the rwl files for the Wazny chronologies and saves as a nested tibble with
metadata.

### 06_RemoveHistorical.R
Removes historical sequences from the ITRDB chronologies.

### 07_MergeDatasets.R
Merges the ITRDB, CTRL, and Wazny chronologies into a single nested tibble. Also
adds columns for coordinates in EPSG 3035.

### 08_CleanChronologies_Step1.R
Cleans the merged dataset in four steps: (1) checks for and removes chronologies 
with coordinates located at sea, (2) combines chronologies located within 1 km of 
each other, (3) checks for and removes problem series/chronologies according to 
the method of van Sluijs et al., (4) removes sites with chronologies shorter than 
40 years and also removes some specific series which could not be detrended using
the spline approach due to negative values. Tests showed that removing these series
had negligible effect on the indices (calculated in script 08).

### 09_CalculateChronologies_Step1.R
Detrends the rwl measurements using a smoothing spline. Then calculates three
mean chronologies: std, res, and ars. ars is Cook's ARSTAN chronology. The ring
width indices for six sites (DEU304, DEU306, DEU312, DEU314, SWIT232, and SWIT317)
have an AR model order of zero which prevents ars being calculated. For these an 
AR model of order 1 was fitted which tests showed resulted in a time series very
closely correlated with the std chronology.

### 10_CalculateInterseriesCorrelations_Step1.R
Calculates pairwise Pearson correlations (with the number of overlapping years)
between all chronologies for each of the std, res, and ars series, and pairwise
Euclidean distances (m) between chronology locations using EPSG 3035 coordinates.
Combines these into a single long-format table with one row per unique site
pair.

### 11_CleanChronologies_Step2.R
Examines the pariwise correlations calculated in script 10. Identifies sites with
correlation greater than 0.7 and determines whether they should be combined or
removed. After removing these reruns the check for problem series using the
method of van Sluijs et al. then saves.

### 12_CalculateChronologies_Step2.R
Same as script 9 but after the second cleaning step.

### 13_CacluateInterseriesCorrelations_Step2.R
Same as script 10 after the second cleaning step. Also calculates a correlation
for each chronology with itself which is the maximum of rbar and the highest
correlation that chronology has with another chronology.

### 14_DefinePredictionGrid.R
Generates a 37.5 Km grid across the study region. Saves it as a raster and a csv
of coordinates.

### 15_GeneratePredictionTable.R
Extracts values from the environmental variables at the grid coordinates defined
in script 14. Saves as a csv.

### 16_GenerateSiteEnvironmentalVariables.R
Extracts values from the environmental variables rasters at the site coordinates
and saves as a csv. Upscales the soil data by a factor of six to better account
for positional uncertainty in the ITRDB coordinates.

### 17_GenerateTrainingTables.R
For each site generate a table of correlations with all other sites along with 
the environmental predictors. Save as nested Rdata tibble and as individual csvs.

### 18_CompareModels.R
For each sites full training table, 6 different possible interpolation models 
are calculated: 1 random forest and 6 universal Kriging R2 and RMSE are 
calculated and the models compared. Based on these results three models are 
tested with leave-one-out cross-validation in the next step: random forest, 
universal Kriging with a second order model on the coordinates, and universal
Kriging with a linear model on the coordinates and also including distance from 
the central site. Saves the resulting maps in Figures/ModelComparisons.

### 19_LOOCV_Models.R
Treats each site in turn as the centre site and builds correlation maps with one
additional site left out each time. Three models tested: random forest and two 
universal Kriging models. Maps generated here will be used for assignment in 
later step.