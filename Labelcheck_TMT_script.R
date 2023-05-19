# Labelcheck script
# 29.Jan 2020

# Note: This script uses msms.txt files as they are created by MQ search
# MQ Search Settings:
# Type: Standard
# Variable Modifications: Modi TMT6plex N-term, Modi TMT6plex K, Oxidation (M), Acetyl (N-term) 

########################################################################################################################################
########################################################################################################################################
# Note: The only thing that needs to be edited for new use of this script is 1) the working directory and 2) the msmsfile name
rm(list=ls())
setwd("D:/Labelcheck 3.0 Th1 January 2020/combined/txt")
list.files()

# read in data as a dataframe
msmsfile <- "msms.txt"                                                              
df <- read.delim(file=msmsfile,sep="\t", header=TRUE,stringsAsFactors=FALSE)
dim(df)

# filter out spectra with lowest 10% scores
summary(df$Score)
cutoff <- quantile(df$Score, probs=0.30)
cutoff
df <- df[df$Score > cutoff,]
dim(df)


# extract rawfile information and rawfile names
rawfiles <- df$Raw.file
rawfile_unique <- unique(rawfiles)
rawfile_unique


# initiate final output vector called lab_eff (labeling efficiency)
lab_eff <- numeric(length(rawfile_unique))
names(lab_eff) <- rawfile_unique


# calculate labeling efficiency for each raw file using a forward loop
for (i in rawfile_unique){
  
  df_i <- df[rawfiles == i,]
  bool_Nterm_noAc <- !grepl(df_i$Modifications, pattern= "Protein N-term")
  df_i <- df_i[bool_Nterm_noAc,]
  
  bool_Nterm_TMT <-  grepl(df_i$Modifications, pattern= "plex-Nter")
  lab_eff[i] <- mean(bool_Nterm_TMT)*100
  
  
  
}

print(lab_eff)

 ###################################################################################################################################
####################################################################################################################################

