# Mutational signature dynamics shaping the evolution of oesophageal adenocarcinoma

## Description

This code has been developed to investigate the prevalence and dynamics of mutational signatures across Barrett Oesophagus, primary tumours and metastases in oesophageal adenocarcinoma.

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

The data employed for this analysis can be found in the **/data** folder, the Source Data file of the corresponding paper or are available upon request.

## Code
The scripts employed in this analysis can be found in the **/code** folder.

### **MutationalSignatureAnalysis:**  
Scripts for the inference of mutational signatures and their clonality.

### **PositiveSelectionAnalysis:**  
Positive selection analysis using the dndscv method for defined subgroups of the cohort.

### **AllOtherAnalyses:**  
Script employed to generate the results and figures presented in the paper, numbered according to the figures in the corresponding paper. They include the following analyses:

**- Linking signatures with genomic markers and pathway activities:**

2_compareSignaturesStages.R

4_oncoprintDrivers_sigcutoff05.R

6_CINandS17.R

6_analyseClonality.tracksig.R

6_signaturePresenceCorrWithProgrammes.R

s1_TMB.R

s21_upsetPlots.R

s22-23_PIK3CA_KRAS_APC.R

s24_clonalityPurity.R

s27_signature17.R

s3_compareSignaturesStages_indels.R


**- Gradient boost and multinomial regression classifiers of disease stages based on mutational signatures:**

3_multipleLinearReg-oct2021.R

3_xgboost_withcv_april2022.R

s25_xgboost_withcv_jan2023_plusindels.R

s26_final_xgboostsimpleWithClonality_primvsmets.R


**- DDR process analyses:**

5_ddr.snvs.2023.R

5_genericDDRpatterns.R

7_SBS30_BER.R

s4-5_MMRandImmunity.R


**- Investigating the relation between mutational signatures and clinical outcomes:**

7_linkTreatmentResponse.R

s14_produceClinicalPlot.R

s15-19_clinicalAssociations_smoking.R

s6_barretts_heterogeneity.R


**- Strelka vs Mutect comparison:**

s29_mutectVSStrelka.R


## Demo / Instructions for use

To run a particular script, navigate to the corresponding folder and run the script using the command **Rscript name_of_script.R**. Outputs include results of statistical tests, plots and newly generated result files.

Expected run times: 
Most scripts should run instantaneously, except for the MutationTimer analysis that will depend on the size of your input dataset. The mutation signature-based classifiers can take 2-5 minutes to run on average.

# Publication

This code has been generated for the following publication: 

**Mutational signature dynamics shaping the evolution of oesophageal adenocarcinoma.** (2023) Sujath Abbas, Oriol Pich, Ginny Devonshire, Shawn A. Zamani, Annalise Katz-Summercorn, Sarah Killcoyne, Calvin Cheah, Barbara Nutzinger, Nicola Grehan, Nuria Lopez-Bigas, OCCAMS Consortium, Rebecca C. Fitzgerald*, Maria Secrier* 


# Copyright
This code is free and is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY. See the GNU General Public License for more details.
