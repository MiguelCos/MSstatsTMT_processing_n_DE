##########################################################################
### Pre-processing and differential analysis of TMT-labelled data ########
### Based on MSstatsTMT for MaxQuant outputs                      ########
### Miguel Cosenza v 0.1                                          ########
##########################################################################

## Required packages ----

### Load required packages or install them if necesary ----

packages <- c("dplyr", "here", "stringr", "tidyr")

if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
   install.packages(setdiff(packages, rownames(installed.packages())))  
} else {
   library(dplyr)
   library(stringr)
   library(tidyr)
   library(here)
}

bioc <- c("MSstats", "MSstatsTMT")

if (length(setdiff(bioc, rownames(installed.packages()))) > 0) {
   BiocManager::install(setdiff(bioc, rownames(installed.packages())))  
} else {

   library(MSstatsTMT)
   library(MSstats)
}

# Requirements: Data to process -----

proteingroups <- read.delim(here::here("proteinGroups.txt"),
                            sep = "\t",
                            stringsAsFactors = FALSE)

evidence <- read.delim(here::here("evidence.txt"),
                       sep = "\t",
                       stringsAsFactors = FALSE)

annotation <- read.delim(here::here("annotation.tsv"),
                         sep = "\t",
                         stringsAsFactors = FALSE) 


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


n_rep_per_cond <- summarized_data %>% group_by(Protein,Condition) %>%
   summarise(n_quantified_per_cond = n())

wider <- pivot_wider(n_rep_per_cond,
                     names_from = Condition,
                     values_from = n_quantified_per_cond)
## Exploratory plots ----

## Profile plots ----

dataProcessPlotsTMT(data.peptide = formated_data,
                    data.summarization = summarized_data,
                    type = "ProfilePlot",
                    address = FALSE,
                    originalPlot = TRUE,
                    summaryPlot = FALSE,
                    which.Protein = "Biognosys_pep-a",
                    width = 21,
                    height = 5)

## Quality control plots ----  

dataProcessPlotsTMT(data.peptide = formated_data,
                    data.summarization = summarized_data,
                    type = "QCPlot",
                    address = FALSE,
                    originalPlot = FALSE,
                    summaryPlot = TRUE,
                    which.Protein = "Biognosys_pep-a",
                    width = 21,
                    height = 5)

## Differential expression analysis ----

## Creating a contrast matrix ----
# this section shows a sample on how to create the comparison matrices
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

if (file.exists(here::here("objects/groupCompar.Rds")) == FALSE){
   diffexpr <- groupComparisonTMT(data = summarized_data,
                                  contrast.matrix = "pairwise",
                                  moderated = TRUE)
   
   saveRDS(object = diffexpr, file = here::here("objects/groupCompar.Rds"))
} else {
   diffexpr <- readRDS(here::here("objects/groupCompar.Rds"))}

# set gene names and descriptions associated to the uniprot IDs
unipr2genename = dplyr::select(evidence, Protein = Leading.razor.protein,
                               Gene.names, Protein.names)

diffexpr2 = left_join(diffexpr, unipr2genename,  by = "Protein") %>% 
   filter(duplicated(Protein) == FALSE)

diffexpr3 = left_join(diffexpr2, wider, by = "Protein")

write.csv(diffexpr3,"results.csv") # save file


## Volcano plots (MSstatsTMT) ----

library(MSstats)   

groupComparisonPlots(data = diffexpr,
                     type = "VolcanoPlot",
                     address = FALSE,
                     which.Comparison = "pairwise")




