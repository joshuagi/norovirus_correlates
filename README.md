# Norovirus correlates of protection analysis pipeline

This repository is to replicate analyses from the manuscript titled "An effective oral norovirus vaccine generates mucosal immunity and reduces viral shedding in a phase 2 placebo-control challenge study". Data and code are provided for the totality of evidence analysis (Figure 2) and the correlates of protection machine learning analysis (Figure 5 and Supplmemental Figure 6).


## Running the analyses

- Clone the repository: 

```
git clone https://github.com/joshuagi/norovirus_correlates.git
```

- Then, update the paths to the data files in `totality_of_evidence.R` and `ml_pipeline.R`. The data files are available in the supplement of the manuscript (Supplementary Data File).

```
immunogenicity <-  read.xlsx("your_file_path_here", sheet = 3)
endpoints_raw_data <- read.xlsx("your_file_path_here", sheet = 2)
```

- Start R in the root directory... 

```
cd norovirus_correlates
R
```

- ... and restore the project environment and run the R scripts:


```
renv::restore()
source("code/totality_of_evidence.R")
source("code/ml_pipeline.R")
```


## File Guide

### code

This directory contains scripts to replicate analyses with the provided data:

- `totality_of_evidence.R` is for the totality of evidence analysis

- `ml_pipeline.R` & `pipeline.utils.R` are for the correlates of protection analysis

### results

This directory contains the output from `ml_pipeline.R`. Model results presented in the manuscript are available; running the analysis pipeline again will save output here.

### figures

This directory contains figures presented in the manuscript; running the analyses above will save the figures here.


