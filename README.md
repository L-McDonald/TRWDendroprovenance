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
Merges the ITRDB, CTRL, and Wazny chronologies into a single nested tibble.