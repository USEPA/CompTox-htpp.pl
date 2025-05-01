HTTP Pipeline Code (http.pl)
===========================

System Requirements
-------------------

**This package provides scripts and APIs in R and works on both Windows and Linux systems**


### R

All code in this package has been tested on R 3.6.0 and R 4.4.1, but any version of R 3.x or 4.x should be compatible. Other package versions may have an impact due to changes in functionality or API.

The following R packages are required to pull existing HTTP data from MongoDB:

+ mongolite - Required for DB access in R
+ jsonlite - Required for DB access and handling JSON files
+ readr - Required to read in tables and delimited files

Additional R packages are required to support pipelining of new HTTP data:

+ doParallel - Required for parallelization of certain functions for users on Linux systems
+ foreach - Required for parallelization of certain functions for users on Linux systems
+ parallel - Required for parallelization of certain functions for users on Linux systems
+ tcplfit2 - Required for concentration-response modeling


### Additional R packages are also required to support formatting, plotting and timing for output legibility:

+ data.table
+ dplyr
+ ggplot2
+ plyr
+ rlist
+ stringr
+ tibble
+ tictoc
+ tidyr


Installing htpp.pl
------------------

The htpp.pl package is an R package and can be installed using the *devtools* R package. R package dependencies, described above, should be installed prior to installing *htpp.pl*


This package can be installed from the public EPA Github or installed locally after cloning the code repository:

**Install from Github**

1. Run devtools function to install directly from Github:
```r
#load devtools library
library(devtools)

#install from Github
devtools::install_github("USEPA/CompTox-htpp.pl")

library(htpp.pl)
```

**Local installation**

1. Clone or download the package from github (TBD):
```bash
git clone  https://github.com/USEPA/CompTox-htpp.pl.git
```
2. Move to the htpp.pl installation location and load the package using *devtools* functions:
```r
#load devtools library
library(devtools)

#use devtools functions to load package and documentation
devtools::install_local(path="path/to/CompTox-htpp.pl/", build_vignettes=TRUE) #installs package and builds vignette

library(htpp.pl)
```


### Version History

**v0.4 (4/30/2025)**

+ Package installation can use `devtools::install_github` and `devtools::install_local`
+ Added DB-free functionality for full htpp.pl -- see vignette for details
+ Updated Vignette
+ Various code improvements and bug fixes

**v0.3 (3/4/2024)**

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

