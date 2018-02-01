
#####################################################
## start vignette
# devtools::use_vignette("systemsbio")
devtools::build()
devtools::build_vignettes()
devtools::clean_vignettes # to remove build tex/pdf files.

# output: rmarkdown::html_vignette
# output: pdf_document

## package dependencies
devtools::use_package("hexbin")
devtools::use_package("arrayQualityMetrics", "Suggests")




######################################################
## Create example data for systemsbio


# use_data_raw() # code for preparing examples not to be included in bundled package
# 
# * Creating `data-raw`.
# * Adding `data-raw` to `.Rbuildignore`.
# Next: 
#   * Add data creation scripts in data-raw
# * Use devtools::use_data() to add data to package


## used samples for example data
ctrl_sample_names <- c("BN0080", "BN0065", "BN0081") # c("BN0063", "BN0088", "BN0064")
case1_sample_names <- c("BN0319", "BN0236", "BN0168")
case2_sample_names <- c("BN0223", "BN0022", "BN0172") # CAD + MI

## projectfolder
projectfolder <- "W:/Projekte/BioNRW/BioNRW_Labor"

# prepare expression data
gex <- read_Illu_expr_array(
            dataFile = file.path(projectfolder, "Illumina GeneExpression/Analysis/SampleProbeProfile.txt"), 
            qcFile = file.path(projectfolder, "Illumina GeneExpression/Analysis/ControlProbeProfile.txt"), 
            sampleSheet = file.path(projectfolder, "Illumina GeneExpression/sample_sheet_BioNRW.csv"),
            ProbeID = "PROBE_ID", 
            skip = 0, 
            controlID= "ProbeID", 
            qc.skip = 0, 
            qc.columns = list(exprs = "AVG_Signal", Detection = "Detection Pval"),
            sampleColumn = "Sample_Name", 
            groupColumn  = "Sample_Group",  
            exprchip= "HumanHT-12 v4",
            org= "human",
            covarfile = file.path(projectfolder, "Probenselektion/BioNRW_cov.txt"), 
            covarsampleID = "ID",   
            matchvar = NULL, 
            method_norm = "none", 
            transform= "log2",  
            fields2Add= NULL 
            ) 
  
# process phenotype data
pData(gex)$Sample_Group <- as.character(pData(gex)$Sample_Group)
pData(gex)[pData(gex)$case_types == "case", "Sample_Group"] <- "case1"
pData(gex)[pData(gex)$case_types == "CAD_MI", "Sample_Group"] <- "case2"
pData(gex) <- pData(gex)[, colnames(pData(gex)) %in% c("Sample_Name", "Sample_Group", "Sentrix_ID", "Sentrix_Position",
                              "Age", "Gender")]


# modulate expression values with random numbers
set.seed(10)
randomValues <- matrix(data= runif(length(exprs(gex)), min=-0.3, max=0.3), nrow= nrow(exprs(gex)), ncol= ncol(exprs(gex)))
#randomValues <- matrix(data= rnorm(length(exprs(gex))), nrow= nrow(exprs(gex)), ncol= ncol(exprs(gex)))
exprs(gex) <- exprs(gex) + exprs(gex) * randomValues


# subset expression set
gex <- gex[, c(case1_sample_names, case2_sample_names, ctrl_sample_names)]
# gex <- gex[1:1000, ] # subset of genes


pData(gex)$Sample_Name <- sub("BN", "sample", pData(gex)$Sample_Name) # rename samples
sampleNames(gex) <- sub("BN", "sample", sampleNames(gex))

# save gex to data folder
use_data(gex, overwrite = TRUE, compress = "bzip2")









