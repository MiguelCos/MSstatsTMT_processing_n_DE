#######################
# create annotation matrix for MSStatsTMT
#######################
# Klemens Froelich 14.5.2019
# Version 2: Includes Fraction, TechRepMixture
# KF 27.5.2019
# for MSstatsTMT Version 1.2.





############# README ##################

# Execute the code LINE BY LINE and read through the instructions.

# This annotation works for all TMT datasets that do not include technical remeasurements of the same high pH HPLC fraction!
# If you have technical replicates of this sort, please write the annotation table in excel

#######################################


#set working directory to current folder of R script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#create annotation dataframe for your first experiment:


###### Run ######
# What is the name of the Raw Files of your first TMT experiment?

# for example:    XY10, XY11, XY12
# is written as:  c(paste("XY", seq(from = 10, to = 12), sep = ""))

#alternative:     c("XY10", "XY11", "XY12")

#template:
Run1 <- c(paste("BM", seq(from = 4742, to = 4753), sep = ""))

###### Channel ######
# how many channels did you use in general? If you have used different numbers of channels in your experiments, please use a 
# MaxQuant Search which always includes all channels and label the empty Channels as "Empty" in conditions

# for example 1:      You have used only 5 channels in 1 TMT experiment
# Output should be: channel.0, channel.1, channel.2, channel.3, channel.4

#can be written as:  c(paste("channel.",-1 + seq(from = 1, to = 5), sep = ""))



# If you have used different channels in different TMT experiments, just use all 11 channels:
# output should be:  "channel.0"  "channel.1"  "channel.2"  "channel.3"  "channel.4"  "channel.5"  "channel.6"  "channel.7" 
#                    "channel.8"  "channel.9"  "channel.10"



#template (for example 1):
Channels <- c(paste("channel.",-1 + seq(from = 1, to = 11), sep = ""))


##### Conditions #####
# what conditions (e.g. Wildtype or Knockout) did you measure in the first experiment? Keep the same order as in the TMT
# channels which you used.

# For example: You used 5 channels. 1st channel: WT, 2nd channel: WT, 3rd channel: KO, 4th channel: KO, 5th channel: WT
# output should be: WT, WT, KO, KO, WT

# Can be written as: c("WT", "WT", "KO", "KO", "WT")
#Normlization channel is written as : "Norm"

#template:
Condition1 <- c("Empty", "Norm", "PatientAM1", "PatientAM2", "PatientAK1", "PatientAK2","cntrl11","cntrl12","cntrl21","cntrl22","Norm")

# If you have used more than 1 experiment and used different channels, add EMPTY as a condition in all non-used channels
# example: You have used the first 5 channels in the first experiment and 5 different channels in the second experiment:
# (always count up the channels: 126 is channel.0 , 127N is channel.1 , 127C is channel.2 etc)
# output should be for first experiment: 

#Condition1 <- c("WT", "WT", "KO", "KO", "WT", "Empty","Empty","Empty","Empty","Empty","Empty")
#Condition2 <- c("Empty","Empty","Empty","Empty", "WT","Empty", "WT", "KO","Empty", "KO", "WT")

###### BioReplicate ###### 
#BioReplicates allows MSStats to decide whether to combine multiple replicates (if they were technical)
# Or whether to use paired statistics, like patient matched settings

# Every biological replicate gets its own number. Technical replicates get the same number.

# example from Condition1 would consist of purely technial replicates and Condition2 would contain real biological replicates

# output should be: 
Bioreplicate1 <- c(100,200,1,2,3,4,5,6,7,8,200)

# Note that the empty channels get an arbitrary Bioreplicate Number, in this example 100. You also need to give an arbitrary replicate number
# to the Normalization channels you used e.g. 200



##### make annotations for each experiment #######
 
#annotation function takes in all vectors and combines them to a ready to use annotation file:
make_annotation <- function(Run, Channel, Condition, Mixture, Bioreplicate) 
{
  Run_vector      <- rep(Run, each = length(Channel))
  Fraction_Vector <- rep(seq(1:length(Run)), each = length(Channel))
  TechRepMixture_V<- Mixture
  Channel_vector  <- rep(Channel, length(Run))
  Condition_vector<- rep(Condition, length(Run))  
  Mixture_vector  <- rep(paste("Mixture", Mixture, sep = ""), length(Run) * length(Channel))
  Biorep_vector   <- rep(Bioreplicate, length(Run))
  df              <- data.frame(Run_vector, Fraction_Vector, TechRepMixture_V, Channel_vector, Condition_vector,
                                Mixture_vector, Biorep_vector)
  colnames(df)    <- c("Run", "Fraction", "TechRepMixture", "Channel", "Condition", "Mixture", "BioReplicate")
  return(df)
}


# for experiment 1 you can combine the information defined above:

#template:
annotation1 <- make_annotation(Run = Run1, Channel = Channels, Condition = Condition1,
                               Mixture = 1, Bioreplicate = Bioreplicate1)

# If you have a second experiment you want to combine in your analysis, repeat the steps above and summarize again:
# Note that you also have to change the Mixture to 2

#template:


#Combine all single experiment annotation files as follows:

#template:
annotation <- rbind(annotation1)

# make annotation output as csv or txt
output.file = "annotation.tsv"
write.table(annotation, output.file, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE, dec= ".")







