
## ReCognitION proteomics analyses are based on the DIA_pooled_DIA_NN datasets: DIA_pooled_DIA_NN
## Phenotype data is based on the OPTIMISTIC study: OPTIMISTIC_Phenotype_Data.xlsx
## External proteomics data: summary-results.proteinpgmaxlfq.parquet
## External ELISA based results of ITIH3: 20231024_ITIH3_ DM1_individual values.xlsx
## External phenotype data: Phenotype_ITIH3_DM1total.xlsx

## Scripts should be run in the order as shown below (so that generated datasets can be used in the next script).

##############################################
## Pre-processing and filtering of datasets ##
##############################################

## Filtering of phenotype data 
Script: Phenotype_filter.R
Input: OPTIMISTIC_Phenotype_Data.xlsx
Output: Philtered_Phenotypedata.RDATA

## Pre-processing of peptide data 
Script: peptide_analysis_p_DIA-NN.R
Input: merged_speclib.pr_matrix.tsv, All samples_LOG_randomize_runorder.xlsx, Filtered_Phenotypedata.RData
Output: normalized_peptide_intensities_p_DIA_NN.RDATA

## Pre-processing of protein data 
Script: protein_analysis_p_DIA-NN.R
Input: normalized_peptide_intensities_p_DIA_NN.RDATA, merged_speclib.pg_matrix.tsv, All samples_LOG_randomize_runorder.xlsx, Filtered_Phenotypedata.RDATA
Output: normalized_protein_intensities_p_DIA_NN.RDATA

## Pre-processing of external validation protein dataset, Canadian DM1 cohort only 
Script: external_validation_analysis.R
Input: summary-results.proteinpgmaxlfq.parquet
Output: ext_val_normalized_proteins.RDATA

## Pre-processing of external validation protein dataset without minimum threshold filtering, including controls 
# Results of this script exclusively used for plot generation
Script: external_validation_analysis_no_intensity_minimum.R
Input: summary-results.proteinpgmaxlfq.parquet
Output: ext_val_normalized_proteins_no_minimum.RDATA


##########################
## Statistical analyses ##
##########################

## Statistical analysis of the protein data with mixed effect models 
Script: Protein_statistics_p_DIA-NN.R
Input: normalized_protein_intensities_p_DIA_NN.RDATA 
Output: DIA-NN_p_Protein_statistics.RDATA, ReCognitION_CTG_hits.csv, ReCognitION_SMWT_hits.csv, ReCognitION_CBT_hits.csv, ReCognitION_GET_hits.csv, ReCognitION_OM_summary.csv

## External validation of statistical results, including ELISA-based ITIH3 validation 
Script: External_validation_statistics.R 
Input: 20231024_ITIH3_ DM1_individual values.xlsx, Phenotype_ITIH3_DM1total.xlsx, ext_val_normalized_proteins.RDATA
Output: all_p_ext_val_res.RDATA, ext_val_CTG_hits.csv, ext_val_SMWT_hits_sig.csv

## Bootstrap enhanced multivariate Elastic-Net modeling of CTG repeat and SMWT, including external validation | 
Script: Multivariate_Elastic_Net_p_DIA_NN.R
Input: all_p_ext_val_res.RDATA, normalized_protein_intensities_p_DIA_NN.RDATA, DIA-NN_p_Protein_statistics.RDATA
Output: p_DIA_NN_CTG&SMWT_prot_VIPs.RDATA

## Mediation analysis studying the impact of BMI on validated biomarkers associated with CTG-repeat and 6MWT
# Needs table 2 in addition to datasets generated above
Script: BMI_Complement_Mediation
Input: DIA-NN_p_Protein_statistics.RDATA, Table2.RDATA, OPTIMISTIC_Phenotype_Data.xlsx
