##########################################################################
### Pre-processing and differential analysis of TMT-labelled data ########
### Based on MSstatsTMT for MaxQuant outputs                      ########
### Miguel Cosenza v 0.1                                          ########
##########################################################################

## Required packages ----

### Load required packages or install them if necesary ----

packages <- c("dplyr", "here", "stringr", "tidyr", "MSstatsTMT")

if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
   install.packages(setdiff(packages, rownames(installed.packages())))  
} else {
   library(dplyr)
   library(stringr)
   library(tidyr)
   library(here)
   library(MSstatsTMT)
}

# Requirements: Data to process -----

proteingroups <- read.delim(here::here("proteinGroups_HeLa_UPS1_marked.txt"),
                            sep = "\t",
                            stringsAsFactors = FALSE)

evidence <- read.delim(here::here("evidence_HeLa_UPS1_marked.txt"),
                       sep = "\t",
                       stringsAsFactors = FALSE)

annotation <- read.delim(here::here("MaxQuant_annotation.csv"),
                         sep = ";",
                         stringsAsFactors = FALSE) %>%
            dplyr::mutate(Run = str_replace(Run, "161117", "161122"))

# Requirements: Your parameters -----

### Script execution start ------------  

# EXECUTE EVERY LINE CONSECUTIVELY 

## Create a folder to store big project-objects (to ) ----

if(dir.exists(here::here("objects")) == FALSE){
      dir.create(here::here("objects"))
}

## MaxQuant to MsstatsTMT formating ----

if (file.exists(here::here("objects/formated_data.Rda")) == FALSE){
      
      formated_data <- MaxQtoMSstatsTMTFormat(evidence = evidence,
                                              proteinGroups = proteingroups,
                                              annotation = annotation,
                                              which.proteinid = 'Leading.proteins',
                                              rmProt_Only.identified.by.site = TRUE)
      save(formated_data,
           file = here::here("objects/formated_data.Rda"))
      
} else {
      load(here::here("objects/formated_data.Rda"))
}

## MSstatsTMT protein summarization ----  

if (file.exists(here::here("objects/summarized_data.Rda")) == FALSE){
   
      summarized_data <- proteinSummarization(formated_data,
                                              method="msstats")

      save(summarized_data,
           file = here::here("objects/summarized_data.Rda"))

} else {
      load(here::here("objects/summarized_data.Rda"))
}

## Exploratory plots ----

## Profile plots ----

dataProcessPlotsTMT(data.peptide = tmt_mstsformat,
                    data.summarization = tmt_procs,
                    type = "ProfilePlot",
                    address = FALSE,
                    originalPlot = TRUE,
                    summaryPlot = FALSE,
                    which.Protein = "O00762ups",
                    width = 21,
                    height = 5)

## Quality control plots ----  

dataProcessPlotsTMT(data.peptide = tmt_mstsformat,
                    data.summarization = tmt_procs,
                    type = "QCPlot",
                    address = FALSE,
                    originalPlot = FALSE,
                    summaryPlot = TRUE,
                    which.Protein = "O00762ups",
                    width = 21,
                    height = 5)

## Differential expression analysis ----

## Creating a contrast matrix ----

comparison1<-matrix(c(-1,1,0,0),nrow=1)
comparison2<-matrix(c(-1,0,1,0),nrow=1)
comparison3<-matrix(c(-1,0,0,1),nrow=1)
comparison4<-matrix(c(0,-1,1,0),nrow=1)
comparison5<-matrix(c(0,-1,0,1),nrow=1)
comparison6<-matrix(c(0,0,-1,1),nrow=1)


comparison_all <- rbind(comparison1, 
                        comparison2, comparison3, 
                        comparison4, comparison5, comparison6)

colnames(comparison_all)<- c("0.125", "0.5", "0.667", "1")

row.names(comparison_all)<-c("0.5-0.125","0.667-0.125","1-0.125",
                             "0.667-0.5","1-0.5","1-0.667")

## Executing `groupComparisonTMT()` ----

if (file.exists(here::here("Data/TMT/groupCompar.Rds")) == FALSE){
   diffexpr <- groupComparisonTMT(data = tmt_procs,
                                  contrast.matrix = comparison_all)
   
   saveRDS(object = diffexpr, file = here::here("Data/TMT/groupCompar.Rds"))
} else {
   diffexpr <- readRDS(here::here("Data/TMT/groupCompar.Rds"))
   
## Volcano plots (MSstatsTMT) ----
   
groupComparisonPlots(data = diffexpr,
                        type = "VolcanoPlot",
                        address = FALSE,
                        which.Comparison = "0.5-0.125",
                        ProteinName = FALSE)
}





