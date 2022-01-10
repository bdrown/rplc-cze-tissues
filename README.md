# Mapping the Proteoform Landscape of Five Human Tissues

This repository contains supporting data analysis for manuscript focused on 
identifying tissue-specific proteoforms and selectivity of reverse phase
liquid chromatography.

## Directory Structure

- `data` - supporting data such as PFRs found in the Human Proteoform Atlas, 
           mapping for specific PTMs to general categories, and list of histones
           curated in TDportal 4.
- `rawtools` - results from RawTools
- `tdreport_dumps` - resutls from TDportal 4
- `figures` - destination for figures
- `output_data` - destination for tables summarizing results

## Data Availability

The input data is too large to be stored in a git repo, so it was deposited in
MassIVE under accession [MSV000088565](ftp://massive.ucsd.edu/MSV000088565/). 
Download the `.tdReport` and `.csv` files from this dataset and put them in 
`tdreport_dumps` directory.

## Data Preprocessing

### Physiochemical Properties

Each row in the "Hits List" `csv` files correspond to each proteoform-spectral
match (PfSM). The physiochemical properties of proteoforms are calculated  by
`calc_protein_props.py`. Python dependencies can be handled with a conda
environment and 

```sh
conda env create -f environment.yml
conda activate hubmap_analysis
make calc_props
```

### RData

The SQL queries of the .tdreport files to count supporting fragment ions take
quite a lot of time, so all the work of merging data and executing queries is
split off into `load_data.R`. It will generate `preprocessed.RData`, which is
used by all the R Notebooks.

```
Rscript load_data.R
```

## Dependencies

### Python

Python dependencies are listed in `environment.yml`. Physiochemical property
calculation requires Python 3.9, pandas, Biopython, and progressbar.

### R 

Most of the figures are generated with R 4.1.0 in `.Rmd` notebooks. 
Dependencies include `Tidyverse`, `VennDiagram`, `RSQLite`, and `Peptides`.

```r
>sessionInfo()
R version 4.1.0 (2021-05-18)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19043)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252 LC_NUMERIC=C                          
[5] LC_TIME=English_United States.1252    

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] Peptides_2.4.4      car_3.0-11          carData_3.0-4       venn_1.10           VennDiagram_1.6.20  futile.logger_1.4.3 forcats_0.5.1       stringr_1.4.0       dplyr_1.0.7        
[10] purrr_0.3.4         readr_1.4.0         tidyr_1.1.3         tibble_3.1.2        ggplot2_3.3.5       tidyverse_1.3.1     RSQLite_2.2.7      

loaded via a namespace (and not attached):
  [1] colorspace_2.0-2       seqinr_4.2-8           ellipsis_0.3.2         rio_0.5.27             XVector_0.32.0         fs_1.5.0               rstudioapi_0.13        farver_2.1.0          
  [9] bit64_4.0.5            fansi_0.5.0            lubridate_1.7.10       xml2_1.3.2             splines_4.1.0          cachem_1.0.5           knitr_1.36             ade4_1.7-18           
 [17] jsonlite_1.7.2         broom_0.7.9            gridBase_0.4-7         Rmpfr_0.8-5            ashr_2.2-47            dbplyr_2.1.1           compiler_4.1.0         httr_1.4.2            
 [25] backports_1.2.1        assertthat_0.2.1       Matrix_1.3-3           fastmap_1.1.0          cli_3.0.1              formatR_1.11           htmltools_0.5.1.1      admisc_0.16           
 [33] tools_4.1.0            gmp_0.6-2              gtable_0.3.0           glue_1.4.2             GenomeInfoDbData_1.2.6 Rcpp_1.0.7             jquerylib_0.1.4        cellranger_1.1.0      
 [41] Logolas_1.3.1          vctrs_0.3.8            Biostrings_2.60.2      ape_5.5                nlme_3.1-152           xfun_0.24              CVXR_1.0-9             openxlsx_4.2.4        
 [49] rvest_1.0.2            lifecycle_1.0.1        irlba_2.3.3            zlibbioc_1.38.0        MASS_7.3-54            scales_1.1.1           hms_1.1.1              parallel_4.1.0        
 [57] RColorBrewer_1.1-2     lambda.r_1.2.4         yaml_2.2.1             curl_4.3.2             memoise_2.0.0          stringi_1.7.3          SQUAREM_2021.1         S4Vectors_0.30.2      
 [65] BiocGenerics_0.38.0    zip_2.2.0              truncnorm_1.0-8        GenomeInfoDb_1.28.4    rlang_0.4.11           pkgconfig_2.0.3        bitops_1.0-7           evaluate_0.14         
 [73] lattice_0.20-44        invgamma_1.1           labeling_0.4.2         bit_4.0.4              tidyselect_1.1.1       magrittr_2.0.1         R6_2.5.1               IRanges_2.26.0        
 [81] generics_0.1.0         DBI_1.1.1              pillar_1.6.4           haven_2.4.1            foreign_0.8-81         withr_2.4.2            mgcv_1.8-35            abind_1.4-5           
 [89] RCurl_1.98-1.5         mixsqp_0.3-43          modelr_0.1.8           crayon_1.4.1           futile.options_1.0.1   utf8_1.2.1             rmarkdown_2.11         readxl_1.3.1          
 [97] data.table_1.14.2      blob_1.2.2             reprex_2.0.1           digest_0.6.27          stats4_4.1.0           munsell_0.5.0 
 ```