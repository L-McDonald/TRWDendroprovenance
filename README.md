# Tree-Ring Width Dendroprovenancing
Project using oak ring-width patterns from across Europe to determine the
provenance of unknown samples.

## Instructions to reproduce

## Scripts
### 00_Functions.R
Various functions used throughout the project

### 01_SpeciesDistributionRanges.R
Plots the species distribution ranges for Quercus robur, Q. petraea, Q. cerris,
and Q. frainetto from
[EuForGen](https://data.mendeley.com/datasets/hr5h2hcgg4/19). Determines the
maximum extent of species occurrence for downloading from the ITRDB.

### 02_DownloadChronologiesITRDB.R
Downloads all ITRDB chronologies of the above species and saves RWL files along
with metadata.

### 03_RemoveHistorical.R
Removes historical sequences from the ITRDB chronologies.