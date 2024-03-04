HTTP Pipeline Code (http.pl)
===========================

System Requirements
-------------------

**This package provides scripts and APIs in R and works on both Windows and Linux systems**


### R

All code in this package has been tested on R 3.6.0, but any version of R 3.x or 4.x should be compatible. Other package versions may have an impact due to changes in functionality or API.

The following R packages are required to pull existing HTTP data from MongoDB:

+ mongolite - Required for DB access in R, code has been tested with v2.1.0.
+ jsonlite - Required for DB access and handling JSON files, code has been tested with v1.6.0
+ readr - Required to read in tables and delimited files, code has been tested with version v2.1.4

Additional R packages are required to support pipelining of new HTTP data:

+ doParallel - Required for parallization of certain functions for users on Linux systems, code has been tested with version v1.0.17
+ foreach - Required for parallization of certain functions for users on Linux systems, code has been tested with version v1.5.2
+ parallel - Required for parallization of certain functions for users on Linux systems, code has been tested with version v3.6.0
+ tcplfit2 - Required for concentration-response modeling, code has been tested with version v0.1.5


### Additional R packages are also required to support formatting, plotting and timing for output legibility:

+ data.table - code has been tested with version v1.14.8
+ dplyr - code has been tested with version v1.1.3
+ ggplot2 - code has been tested with version v3.4.3
+ plyr - code has been tested with version v1.8.8
+ rlist - code has been tested with version v0.4.6.2
+ stringr - code has been tested with version v1.5.0
+ tibble - code has been tested with version v3.2.1
+ tictoc - code has been tested with version v1.2
+ tidyr - code has been tested with version v1.3.0


Installing htpp.pl
------------------

The htpp.pl package is currently implemented as a proto-R package and can be installed using the *devtools* R package.

1. Ensure that R and Rstudio are properly installed (including Rtools)
2. Clone or download the package from github (TBD):
```bash
git clone (TBD)
```
3. Move to the htpp.pl installation location and load the package using *devtools* functions in R/Rstudio:
```r

#load devtools library
library(devtools)

#change working directory to htpp.pl repo
setwd("path/to/htpp.pl/")

#use devtools functions to load package and documentation
devtools::load_all() #loads the htpp.pl package similar to library(htpp.pl)
devtools::document() #builds documentation, creates
devtools::build_vignettes() #builds package vignette in '/htpp.pl/doc'
devtools::build_manual() #build manual and saves PDF in parent directory of the package
```


### Version History

**v0.3 (TBD)**

+ Various code improvements
+ Improved roxygen2 documentation
+ devtools package load functions work as intended
+ Added vignette
+ Added README


**v0.2 (11/2/2023)**

+ Various code improvements throughout
+ Add parallel capabilities to several functions
+ Added external data
+ Refined documentation


**v0.1.1 (11/1/2023)**

+ roxygenized all functions and formatted for R package
+ General code improvements throughout

**v0.1 (8/3/2023)**

+ Initial version of code adapted from J. Nyffeler

