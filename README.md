
## Evolution of the nonsense mediated decay pathway is associated with decreased cytolytic immune infiltration.
Authors: Boyang Zhao and Justin Pritchard

### File info for folder `src`
* analyze.R: random forest model building and analyses
* analyze_indv.R: individual ad-hoc analyses
* analyze_mod.R: additional analyses of the random forest model results
* analyze_pan.R: aggregate pan-cancer analyses
* analyze_nullmodel.R: random forest null model building
* analyze_nullmodel_results.R: analyze null model results
* analyze_indv_globalNMD.R: individual ad-hoc analyses, for global NMD
* analyze_indv_survivalSKCM.R: individual ad-hoc analyses, for SKCM survival
* nmd_compr.R: NMD alterations analyses, association with NMD metrics
* nmd_compr_triple.R: NMD alterations analyses, broken to compare none, 1, or >1 variants
* nmd_compr_cyt.R: NMD alterations analyses, association with cytolytic activity
* utils_surv.R: survival methods used for analyses
* utils.R: methods used for analyses

### File info for folder `datasets`
The folders are organized by indication, and downloaded from Broad firehose. 
* *_clin contains clinical data
* *_cna contains CNA data
* *_mut contains MAF mutational data
* *_rnaseq contains RNA-seq data

### File info for folder `outputs`
* 18.0604: random forest model results
* 18.1212_survival: survival analyses. Parsed .rda files from 18.0604 were copied over to this folder, and survival analyses were then run
* 19.0108 cbio: genetic alterations of NMD and other pathways, results from cBioPortal
* 19.0108_cnacutoff_1_NMDpos_ampdel: cnv of genes in the NMD pathway across indications
* 19.0110_cnacutoff_1_NMDpos: association of NMD genes (any alterations) with NMD metrics, this is used for the subsequent analyses below
* 19.0110_cnacutoff_1_NMDpos_CooccurCyt: association of NMD co-alterations and cytolytic activity
* 19.0110_cnacutoff_1_NMDpos_CooccurVs1variantVsCtrl: trend analyses of NMD genetic alterations and NMD metrics
* 19.0110_cnacutoff_1_NMDpos_CooccurVsCtrl: association of NMD genes (co-alterations) with NMD metrics
* 19.0817_nullmodel: results from null model (randomized inputs)
* 19.0819_cox: residual results from Cox regression models
* 19.0820_NMDglobal: generate visualizations for NMD global differences plots
