# Mutational signature dynamics across the disease course in oesophageal adenocarcinoma

## Description

This code has been developed to investigate the prevalence and dynamics of mutational signatures across Barrett Oesophagus, primary tumours and metastases is oesophageal adenocarcinoma.

## System Requirements
Operating system(s): Unix (Linux, Mac OS)
Programming Language: R 4.1.2

Package versions utilised:

desconstructSigs v1.9.0

MutationalPatterns v3.10.0

MutationTimer v1.00.0 

dNdScv v0.1.0 

xgboost v1.7.3.1

randomForest v4.7-1.1

glmnet v4.1-7

## Installation guide
Make sure you are running R version 4.1.2 or above, and install the required package dependencies as indicated in the System Requirements section above.

Typical install time on a "normal" desktop computer: negligible, mostly reliant on dependent packages.

## Data

The data employed for this analysis can be found in the **/data** folder or is available upon request.

## Code
The scripts employed in this analysis can be found in the **/code** folder as follows:

-**MutationalSignatureAnalysis:**  
Scripts for the inference of mutational signatures and their clonality.

-**PositiveSelectionAnalysis:**  
Positive selection analysis using the dndscv method for defined subgroups of the cohort.

-**GenomicAssociations:**  
Script employed to link signatures with genomic markers and pathway activities.

-**ModellingDiseaseStages:**  
Scripts employed for the gradient boost and multinomial regression classifiers of disease stages based on mutational signatures.

-**DDRsignatures:**  
NMF methodology to infer pathway-level signatures of mutations across DDR processes.

-**SignaturePrognosis:**  
Scripts employed to investigate the relation between mutational signatures and clinical outcomes.

## Demo / Instructions for use

To run a particular script, navigate to the corresponding folder and run the script using the command **Rscript name_of_script.R**. Outputs include results of statistical tests, plots and newly generated result files.

Expected run times: 
Most scripts should run instantaneously, except for the MutationTimer analysis that will depend on the size of your input dataset. The mutation signature-based classifiers can take 2-5 minutes to run on average.

# Copyright
This code is free and is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY. See the GNU General Public License for more details.
