#https://github.com/supatt-lab/SingCellaR-myelofibrosis
#https://supatt-lab.github.io/SingCellaR.Doc/SingCellaR_workflow_TARGET_Seq.html
setwd("/Library/Frameworks/R.framework/Versions/4.0/Resources/library/projects/")
install.packages("shiny")
BiocManager::install("SingleCellExperiment")
install.packages("Matrix")
install.packages("matrixStats")
install.packages("igraph")
library(devtools)
install_github("jokergoo/ComplexHeatmap")
install.packages("circlize")
install.packages("DT")

library(shiny)
library(SingleCellExperiment)
library(Matrix)
library(matrixStats)
library(ggplot2)
library(gridExtra)
library(igraph)
library(ComplexHeatmap)
library(circlize)
library(DT)

##Please extract the data and move the folder "projects" into the main package folder!.
find.package('shiny')

load(file = "Psaila_et_al.rdata")
library(shiny)
runApp()
