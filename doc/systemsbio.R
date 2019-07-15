## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
knitr::opts_chunk$set(fig.width=12, fig.height=8)

## ------------------------------------------------------------------------
library("systemsbio")
library("knitr")

## ------------------------------------------------------------------------
library(beadarray)
library(Biobase)
beadarray::exprs(gex)[1:5, 1:5]
Biobase::pData(gex)[,1:5]
Biobase::fData(gex)[1:5, ]

## ---- echo = TRUE, out.width="47%"---------------------------------------
deg_images<-c("figures/example_boxplot_controlprofile.png",
                         "figures/example_SampleDendrogram_noNorm_adjacency.png")
include_graphics(deg_images)

## ---- echo = TRUE, out.width="47%"---------------------------------------
deg_images<-c("figures/example_pcainfoplot_PC1_2.png",
              "figures/example_pcaplot3d_PC1_2_3.png")
include_graphics(deg_images)

## ---- echo = TRUE, out.width="47%"---------------------------------------
deg_images<-c("figures/module_boxplot_1_of_3.png",
              "figures/ModuleDendrogram_Block.png")
include_graphics(deg_images)

## ---- echo = TRUE, out.width="47%"---------------------------------------
deg_images<-c("figures/Heatmap_Module-Groupset_relationship.png",
              "figures/Intramodular_analysis_case1-control.png")
include_graphics(deg_images)

## ---- echo = TRUE, out.width="31%"---------------------------------------
deg_images<-c("figures/Heatmap_example_case1-control.png",
              "figures/Venn_Diagram_example_allcomp.png",
              "figures/Volcano_example_case1-control.png")
include_graphics(deg_images)

## ---- echo = TRUE, out.width="47%"---------------------------------------
Overrep_BP <- read.table("figures/example_case1-control_Overrep_BP_resulttable.txt", sep="\t", header=T)
head(Overrep_BP)

deg_images<-c("figures/example_case1-control_Overrep_BP_cnetPlot.png",
              "figures/example_case1-control_Overrep_BP_enrichmentMap.png")
include_graphics(deg_images)

## ---- echo = TRUE, out.width="47%"---------------------------------------
tf <- read.table("figures/example_case1-control_PWMEnrich_report.txt", sep="\t", header=T)
head(tf)

deg_images<-c("figures/example_case1-control_PWMEnrich_top.png",
              "figures/example_case2-control_PWMEnrich_top.png")
include_graphics(deg_images)

