

## Description
# Set parameters used in analysis pipeline 

## Usage
setparam <- function(workingdir,
                     outdir="Analysis",                    
                     projectname ="", 
                     org = "human", 
                     threshold_p = 0.05, 
                     threshold_fdr = 0.05, 
                     threshold_FC = 1.5, 
                     
                     # GEX parameter
                     exprchip=NULL, 
                     exprchip.manufacturer = "Illumina",
                     outdirGEX = file.path(outdir, "GEX"),
                     sampleGEX = "Sample_Name", 
                     groupGEX = "Sample_Group",
                     matchvarGEX = NULL,
                     covarGEXSampleID ="ID",
                     groupcomparisonsGEX = c("case-control"),

                     # Methylation parameter:
                     methPlatform = NULL, 
                     methPlatform.manufacturer = "Illumina",
                     methannotation.file = NULL,
                     outdirMT = file.path(outdir, "MT"),
                     threshold_beta = 0.1,
                     refgroupMeth = "control",
                     typeMethvar = "categorical",  
                     sampleMT = "Sample_Name", 
                     groupMT = "Sample_Group",
                     matchvarMeth = "none", 
                     covarMTSampleID = "IID",  
                     groupcomparsionsMT = groupcomparisonsGEX,

                     # GT Parameter
                     snpchip=NULL, 
                     snpchip.manufacturer = "Illumina",
                     outdirGT = file.path(outdir, "GT"),
                     covarGTSampleID ="IID",
                     
                     ### other parameter
                     updateBiocPackages = FALSE, # automatically updates all installed Bioconductor Packages
                     source.bioc="http://bioconductor.org/biocLite.R"
                     ) {

  
  ## Arguments
  # workingdir: character with working directory.
  # outdir: character with output directory within working directory for result files (will be created if not existing).                   
  # projectname: character used as suffix for output file names.
  # org: character with species name. Currently "human", "rat" or "mouse". Corresponding annotation packages will be loaded.
  # threshold_p: numeric p-value threshold. 
  # threshold_fdr: numeric threshold for FDR (false discovery rate) corrected p-value.
  # threshold_FC: numeric foldchange threshold. 
  # 
  # # GEX parameter (only relevant if exprchip != NULL)
  # exprchip: character with expression chip type, e.g. "HumanHT-12 v4", "MouseWG-6 v2", "MouseRef-8 v1".
  # exprchip.manufacturer: character with chip manufacturer, e.g. "Illumina" or "Affymetrix". Will be used 
  #                        to load corresponding annotation package. NULL if no supplier specified.
  # outdirGEX: character with output directory for gene expression results.
  # sampleGEX: character with column name containing sample names in sample data (e.g. in Illumina sample sheet). 
  # groupGEX: character with column name containing group names in sample data (e.g. in Illumina sample sheet).
  # matchvarGEX: In case of matched study design, character with column name indicating corresponding samples.
  #              If NULL, unmatched study design assumed.
  # covarGEXSampleID: character with column name of sample names in optional covariate file.
  # groupcomparisonsGEX: character vector with desired group comparisons to analyse in format "groupA-groupB" or 
  #                      nested comparisons in format "(groupA-groupB)-(groupC-groupD)". Group names must correspond
  #                      to data given by 'groupGEX'. Analysis can be performed either for any number of unpaired 
  #                      group comparisons or for one group comparison in paired design.
  # 
  # # Methylation parameter (only relevant if methPlatform != NULL)
  # methPlatform: character with methylation array platform, e.g. "450k-UCSC" for UCSC CpG Islands for 450k array probes, 
  #               "450k-HMM" for HMM CpG Islands for 450k array probes, "27k" for UCSC CpG Islands for 27k array probes.
  # methPlatform.manufacturer: character with chip manufacturer, e.g. "Illumina".
  # methannotation.file: (if no methPlatform specified) character with path to custom tab-text file for annotation 
  #                       with header "SiteID  Chr  Loc  Gene  Island".
  # outdirMT: character with output directory for gene methylation results
  # threshold_beta: numeric threshold for beta value. beta is defined as: methylated signal / (methylated + unmethylated signal).
  # refgroupMeth: character with reference group of methylation analysis.
  # typeMethvar: character with type of methylation variable ("categorical" or "continuous")
  # sampleMT: character with column name containing sample names in sample data. 
  # groupMT: character with column name containing group names in sample data.
  # matchvarMeth: In case of matched study design, character with column name indicating corresponding samples. "none" if unpaired Samples
  # covarMTSampleID: character with column name of sample names in optional covariate file.
  # groupcomparsionsMT = character vector with desired group comparisons to analyse.
  # 
  # # GT Parameter (only relevant if snpchip != NULL)
  # snpchip: character with genotyping chip type, e.g."HumanOmniExpress-24", "HumanOmni5-Quad", "HumanOmni5Exome", "HumanCore-24", "HumanCoreExome-24".
  # snpchip.manufacturer: character with chip manufacturer, e.g. "Illumina".
  # outdirGT: character with output directory for genotyping results.
  # covarGTSampleID: character with column name of sample names in optional covariate file.
  # 
  # ### other parameter
  # updateBiocPackages: boolean. if TRUE automatically updates all installed Bioconductor Packages.
  # source.bioc: character with bioconductor source script.
  
  
  ## Details
  # This functions optionally organizes parameter for analysis settings frequently used within analysis pipeline and
  # returns a list with relevant parameters, which may be used as standard input for further functions. 
  # The denoted working directory is set and required analysis and annotation packages are installed and/or loaded.
 
  
  ## Value
  # List with parameters given or derived from input arguments. A log-file showing these parameter is stored 
  # in the output directory as side-effect.
  
  
  ## Author(s) 
  # Frank Rühle 
  
  
 
  # Creating directory if not yet existing 
  if (!file.exists(file.path(workingdir))) {dir.create(file.path(workingdir), recursive=T) }
  setwd(workingdir)
  if (!file.exists(file.path(outdir))) {dir.create(file.path(outdir), recursive=T) }
  if (!is.null(exprchip) & !file.exists(file.path(outdirGEX))) {dir.create(file.path(outdirGEX), recursive=T) }
  if (!is.null(methPlatform) & !file.exists(file.path(outdirMT))) {dir.create(file.path(outdirMT), recursive=T) }
  if (!is.null(snpchip) & !file.exists(file.path(outdirGT))) {dir.create(file.path(outdirGT), recursive=T) }
  
  
  # shall all bioconductor packages be updated?
  if(updateBiocPackages) {
    cat("\nUpdating all installed Bioconductor Packages\n")
    source(source.bioc)
    biocLite("BiocUpgrade")
    biocLite(ask=FALSE) 
  }


  # pcksCRAN <- c("base64","flashClust", "matrixStats", "curl", "chron", "madsim", "survival", "plyr", "gtools", "digest", "colorspace", "gsl",
  #               "leaps", "Hmisc", "caTools", "lme4", "EMA", "WGCNA", "RBGL", "VennDiagram", "KernSmooth",
  #               "gplots", "ggplot2", "rJava", "MatrixEQTL", "sem", "doParallel", "reshape2",
  #               "devtools", "rstudioapi", "GenABEL", "stringi", "cluster", "ICC", "rgl")
  # pcksBioc <- c("Biobase", "biomaRt", "arrayQualityMetrics", "limma", "sva", "minet", "siggenes", "genefilter", 
  #               "Category", "GO.db","GOstats", "annotate", "impute", "beadarray","RDAVIDWebService", "GSEABase", 
  #               "pcaGoPromoter", "base64enc", "KEGG.db", "KEGGgraph", "BiocParallel", "minfi", "COHCAP", "S4Vectors", "quantsmooth", 
  #               "MotifDb", "MotIV", "PWMEnrich", "seqLogo", "motifStack", "Biostrings", "GenomicFeatures", "OmicCircos",
  #               "AnnotationHub", "IlluminaHumanMethylation450kmanifest", 
  #               "PWMEnrich.Hsapiens.background", "org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg19.knownGene", "BSgenome.Hsapiens.UCSC.hg19",
  #               "FDb.InfiniumMethylation.hg19",
  #               "affy", "affyQCReport", "hgu95av2.db")
  ##derzeit nicht verfügbar: "Rplinkseq"

  # Essential packages to load
  pcksCRAN <- c("gplots", "ggplot2", "rJava", "VennDiagram", "flashClust")
  pcksBioc <- c("Biobase", "limma")
  attach_package(pcksCRAN, pcksBioc, source.bioc=source.bioc)


  
  ### Load genomic annotation package (includes OrgDb, TxDb and GODb)
  annopkg <- switch(org, human = "Homo.sapiens", 
                         mouse = "Mus.musculus",
                         rat="Rattus.norvegicus")
  attach_package(pkg.bioc=annopkg, source.bioc=source.bioc)
    
  
  ### Load array annotation package:
  # Expression Chip
  if(!is.null(exprchip)) {
       if(tolower(exprchip.manufacturer) == "illumina") {  # for Illumina Arrays
          exprchip <- tolower(exprchip)
          exprchipBez     <- sub("-.*", "", exprchip)
          exprchipVersion <- sub(".*v", "v", exprchip)
          ArrayAnnotation.GEX <- switch(org, human = paste("Human", exprchipVersion, sep=""), 
                                        mouse = paste("Mouse", exprchipVersion, sep=""))
          # ArrayAnnotation.GEX: character string specifying the name of the annotation package, e.g. "Humanv4".
          }  else { # end Illumina

          if(tolower(exprchip.manufacturer) == "affymetrix") {  # for Affymetrix arrays
            ArrayAnnotation.GEX <- sub("cdf$", "", exprchip)
            } else { # end Affymetrix
        
              ArrayAnnotation.GEX <- NULL  # no supplier specified
              }
          }
      }

  # SNP-Chip
  if(!is.null(snpchip)){
      if(tolower(snpchip.manufacturer) == "illumina") { # for Illumina Arrays
          # not implemented yet
          } else { # end Illumina
        }
    }



  
  
  ############ return object with parameter settings
  pipepar <- list()
  
    pipepar[["workingdir"]] <- workingdir
    pipepar[["outdir"]] <- outdir
    pipepar[["projectname"]] <- projectname
    pipepar[["org"]] <- org
    pipepar[["annopkg"]] <- annopkg
    pipepar[["annopkgResources"]] <- names(get(annopkg)@resources)
    pipepar[["threshold_p"]] <- threshold_p
    pipepar[["threshold_fdr"]] <- threshold_fdr
    pipepar[["threshold_FC"]] <- threshold_FC
  
  
  # Expression Parameter
  if(!is.null(exprchip)) {
    pipepar[["exprchip"]] <- exprchip
    pipepar[["ArrayAnnotation.GEX"]] <- ArrayAnnotation.GEX
    pipepar[["outdirGEX"]] <- outdirGEX
    pipepar[["sampleGEX"]] <- sampleGEX
    pipepar[["groupGEX"]] <- groupGEX
    pipepar[["matchvarGEX"]] <- matchvarGEX  
    pipepar[["covarGEXSampleID"]] <- covarGEXSampleID   
    pipepar[["groupcomparisonsGEX"]] <- groupcomparisonsGEX
  }
  
  # Methylation Parameter
  if(!is.null(methPlatform)) {
    pipepar[["outdirMT"]] <- outdirMT
    pipepar[["threshold_beta"]] <- threshold_beta
    pipepar[["refgroupMeth"]] <- refgroupMeth
    pipepar[["typeMethvar"]] <- typeMethvar  
    pipepar[["methPlatform"]] <- methPlatform  
    pipepar[["methannotation.file"]] <- methannotation.file   
    pipepar[["sampleMT"]] <- sampleMT 
    pipepar[["groupMT"]] <- groupMT 
    pipepar[["matchvarMeth"]] <- matchvarMeth 
    pipepar[["covarMTSampleID"]] <- covarMTSampleID  
    pipepar[["groupcomparsionsMT"]] <- groupcomparsionsMT
  }
  
  # Genotyping Parameter
  if (!is.null(snpchip)){
    pipepar[["snpchip"]] <- snpchip
    pipepar[["outdirGT"]] <- outdirGT
    pipepar[["covarGTSampleID"]] <- covarGTSampleID
  }


  
  ############ write parameter settings to screen and file
    cat("\nWriting analysis parameter to", file.path(workingdir, outdir, paste0("Analysis_parameter_", Sys.Date(), ".txt")), "\n")
    sink(file=file.path(workingdir, outdir, paste0("Analysis_parameter_", Sys.Date(), ".txt")), split=T)
  
    # print provided parameter:
    cat(paste("Pipeline parameter set", Sys.Date()),
        "\nSet workingdir to: ", pipepar[["workingdir"]],
        "\noutdir: ", pipepar[["outdir"]],
        "\nprojectname: ", pipepar[["projectname"]],
        "\norg: ", pipepar[["org"]], 
        "\nannotation package: ", pipepar[["annopkg"]],
        "\nannotation package resources: ",  paste(pipepar[["annopkgResources"]], collapse=", "),
        "\nthreshold_p: ", pipepar[["threshold_p"]],
        "\nthreshold_fdr: ", pipepar[["threshold_fdr"]],
        "\nthreshold_FC: ", pipepar[["threshold_FC"]]
    )
    
    if(!is.null(exprchip)) {    
    cat("\n\n Gene Expression parameter:",
        "\noutdirGEX: ", pipepar[["outdirGEX"]],
        "\nExpression chip type: ",pipepar[["exprchip"]], 
        "\nArrayAnnotation.GEX: ", exprchip.manufacturer, pipepar[["ArrayAnnotation.GEX"]],
        "\noutdirGEX: ", pipepar[["outdirGEX"]] ,
        "\nmatchvarGEX: ", pipepar[["matchvarGEX"]]  ,
        "\ncovarGEXSampleID: ", pipepar[["covarGEXSampleID"]],   
        "\nsampleGEX: ", pipepar[["sampleGEX"]], 
        "\ngroupGEX: ", pipepar[["groupGEX"]], 
        "\ngroupcomparisonsGEX: ", pipepar[["groupcomparisonsGEX"]], "\n") 
    }
    
    if(!is.null(methPlatform)) {
    cat("\n\n Methylation Parameter:",
        "\noutdirMT: ", pipepar[["outdirMT"]],
        "\nthreshold_beta: ", pipepar[["threshold_beta"]],    
        "\nmethPlatform: ", pipepar[["methPlatform"]],  
        "\nmethannotation.file: ",  pipepar[["methannotation.file"]],
        "\nrefgroupMeth: ", pipepar[["refgroupMeth"]],
        "\ntypeMethvar: ", pipepar[["typeMethvar"]],  
        "\nmatchvarMeth: ", pipepar[["matchvarMeth"]],
        "\ncovarMTSampleID: ", pipepar[["covarMTSampleID"]],
        "\nsampleMT: ", pipepar[["sampleMT"]], 
        "\ngroupMT: ", pipepar[["groupMT"]],
        "\ngroupcomparsionsMT: ", pipepar[["groupcomparsionsMT"]], "\n") 
    }
    
    if (!is.null(snpchip)){
    cat("\n\n Genotyping Parameter:", 
        "\noutdirGT: ", pipepar[["outdirGT"]],
        "\nGenotyping chip type: ",pipepar[["snpchip"]], 
        "\ncovarGTSampleID: ", pipepar[["covarGTSampleID"]], "\n")
        }

  sink() # return to standard output

  return(pipepar)
}





