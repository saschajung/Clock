# Transcriptional Aging Predictor
This repository stores the scripts used for building the transcriptional aging predictor used in Plesa AM et al. Transcriptomic reprogramming screen identifies SRSF1 as a regulator of aging. Nature Aging. 2026. The paper can be read here: https://rdcu.be/AzbXBE4JzJ6a

# Reproducibility Information

## R and package versions

The results reported in the manuscript were generated and tested using R 4.2.1 on a MacBook Pro (2019) with an Intel Core i5 processor, using the following package versions:

- R 4.2.1
- sva 3.46.0
- h2o 3.44.0.3
- DESeq2 1.38.3
- plotly 4.12.1
- plyr 1.8.9
- stringr 1.6.0
- dplyr 1.2.1
- readr 2.2.0

**Important considerations for reproducibility**

**R version and numerical differences:** R 4.2.1 is the environment used to generate and test the results reported in the manuscript. Beginning with R 4.4.0, changes to `stats::density()` affect the calculation of surrogate variables in the SVA step of our pipeline. In comparisons between the tested environments, differences in the surrogate variables were on the order of 1e-6 to 1e-5 and were amplified in the corrected gene expression matrix, with maximum absolute differences reaching approximately 1e-1. These differences propagate to the estimated coefficients of the transcriptional aging predictor, although the overall model performance appeared stable in our comparisons.

For the closest reproduction of the published results, we recommend using R 4.2.1 together with the package versions listed above. Other combinations of R and package versions have not necessarily been validated but may be permissible if they are compatible with R versions < 4.4.0.

**Processor architecture:** Small numerical differences may also arise from the processor architecture and numerical computations. For example, in our comparisons between Intel and Apple Silicon systems, differences in TPM values reached approximately 1e-14 in absolute terms. These differences had negligible effects on the estimated surrogate variables and predictor coefficients in our tests.

### biomaRt compatibility (September 2026)

Changes to the Ensembl BioMart service in September 2026 affected queries made using older versions of the biomaRt package. In our environment, versions below 2.68.1 could no longer perform the required queries, whereas version 2.68.1 requires R >= 4.5.0 and is therefore incompatible with the R version used for the published analyses.

To avoid this dependency, the Ensembl lookup table was downloaded and saved as EnsemblIdToSymbol_Sep2025.rds. With the exception of the obsolete `GetGenesForRemoval.R` script, the remaining scripts now use this lookup table rather than querying BioMart directly. The generated genes_to_filter.txt file is also included in the repository, so GetGenesForRemoval.R does not need to be run.

## Training and application data availability

The data required to reproduce the analyses reported in the manuscript are available from [Zenodo](https://zenodo.org/records/23570472).

Download the archive and extract its contents into the main repository directory, alongside the scripts. The archive contains the required 'Data' and 'ProcessSelectionData' subdirectories.

## Note on validation data performance

To avoid evaluating the predictor on samples that had already contributed to process selection, the reported validation performance was calculated after excluding samples from the 'Fleischer' and 'LIB049469' libraries, which were used during process selection.

The validation metric stored in the fitted model object still includes these samples and therefore does not correspond to the validation performance reported in the manuscript. The reported performance is based on the recalculation excluding these samples.
