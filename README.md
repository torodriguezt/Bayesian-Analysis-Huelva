# Bayesian Analysis Huelva

Bayesian spatial analysis of the Huelva topsoil survey. For each contamination
source (element signature), the models separate the part of the signal that is
mappable at the survey scale (a smooth Matern field) from the patchy,
site-level part, and assess what this implies for compliance with the
reference levels (NGR) of Decreto 18/2015. Models are fitted with INLA.

## Requirements

R 4.4 (tested with R 4.4.2, INLA 24.12.11, fmesher 0.5.0).

```r
install.packages(c("dplyr", "tidyr", "readr", "sf", "Matrix", "gstat", "units",
                   "ggplot2", "patchwork", "ggspatial"))
install.packages("INLA", repos = c(getOption("repos"),
                 INLA = "https://inla.r-inla-download.org/R/stable"),
                 dependencies = TRUE)
```

## Run

From the repository root, in a fresh R session:

```r
source("main.R")
```

`main.R` runs five stages in order:

| Stage | Script | Output |
| --- | --- | --- |
| Data preparation | `R/Data_preparation/Data_preparation.R` | `results/data/` |
| Exploratory analysis | `R/Exploratory/Exploratory.R` | `results/tables/exploratory/` |
| Models F0-F3 | `R/All_Models/All_Models.R` | `results/fits/` |
| Simulation check | `R/Robustness/Robustness.R` | `results/tables/simulation/` |
| Tables and figures | `R/Results/Results.R` | `results/tables/`, `figures/` |

Fitting F3 can take a long time and a lot of memory.

## Data

`data/raw/` holds the input shapefiles: sample results for Huelva, Aroche and
Campillo, the prediction grid, the mesh boundary, and land-use distance layers.
Paths and study settings (elements, NGR thresholds, CRS) are in
`R/Data_preparation/config.R`.

## License

Copyright (C) 2026 Anonymous Authors. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
