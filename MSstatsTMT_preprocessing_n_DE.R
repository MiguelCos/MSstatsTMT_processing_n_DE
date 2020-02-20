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

