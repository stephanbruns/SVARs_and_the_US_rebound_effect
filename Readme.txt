This is the replication package of "Estimating the Economy-Wide Rebound Effect Using 
Empirically Identified Structural Vector Autoregressions"

This package contains data and code to replicate the findings reported in the article.

#------------
# Data
#------------
The data folder contains one csv with monthly data and one csv with quarterly data. 
These time series are constructed as outlined in the data appendix and they are 
deseasonalized using the x11 procedure in RATS (except for population). 

Data before these data processing steps is available from the authors. 

#------------
# Code
#------------
# Data overview (Figure 3)
This code produces Figure 3.

# Three_variables_dcov_ngml.r
This code produces the estimates for the SVARs with three variables and the 
identification methods: dcov and ngml (monthly and quarterly data).

# Five_variables_dcov_ngml.r
This code produces the estimates for the SVARs with five variables and the
identification methods: dcov and ngml (monthly and quarterly data).

# Three_variables_fastica.R
This code produces the estimates for the SVARs with three variables and the 
identification method: fastICA (monthly and quarterly data).

# Three_variables_hddcdd_fastica.R
This code produces the estimates for the SVARs with three variables, with the 
addition of two exogenous variables: HDD and CDD, and the 
identification method: fastICA (monthly data).

# Five_variables_fastica.r
This code produces the estimates for the SVARs with five variables and the
identification method: fastICA (monthly and quarterly data).

# Five_variables_hddcdd_fastica.R
This code produces the estimates for the SVARs with five variables, with the 
addition of two exogenous variables: HDD and CDD, and the 
identification method: fastICA (monthly data).

# Three_variables_lingam.R
This code produces the estimates for the SVARs with three variables and the 
identification method: LiNGAM (monthly and quarterly data).

# Five_variables_lingam.R
This code produces the estimates for the SVARs with five variables and the 
identification method: LiNGAM (monthly and quarterly data).

#------------
# Results
#------------
This folder contains *.RDATA files that are produced by the code and then
used as an input for other parts of the code.
This folder also contains some of the final results, but most of the final
results are directly shown in R.