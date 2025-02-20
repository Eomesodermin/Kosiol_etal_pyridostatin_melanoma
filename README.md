# Stabilizing G-quadruplexes with Pyridostatin in murine melanoma augments response to checkpoint blockade by modulation of the tumor microenvironment

## Project Summary

Understand the .....

## To-Do

-   [ ] Placeholder

## Dataset Info

-   FILL OUT

    
## Development Notes

-   FILL OUT

## File Structure

The repository contains the following top-level directories:

- **`data/`** _(ignored by Git)_  
  - Contains output from `cellranger`  
  - Will be uploaded to Zenodo for accessibility. (already on zenodo)

- **`results/`** _(ignored by Git)_  
  - Stores processed analysis results.  
  - Not tracked in version control.

- **`saves/`** _(ignored by Git)_  
  - Contains intermediate objects such as Seurat objects.  
  - These files will also be uploaded to Zenodo. (already on zenodo)

- **`scripts/`**  
  - Includes all analysis and processing scripts.
  - Tracked in version control.
    - **`01_Preprocessing.Rmd`**
      - Reads in `cellranger` output and performs the following 
      - QC, Ambient RNA removal, Doublet detection, Normalisation, Integration, Dim Reduction, cluster calling

## Data availability

- Raw Data (Upload to GEO and provide link)  
- `cellranger` output (upload zenodo link)
- `seurat.objects` (upload zenodo link)

## Author Information

-   [Dillon Corvino](https://github.com/Eomesodermin)
