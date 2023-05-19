# script reduce redundant entries in evidence table, change retention time to median of all evidences found for a unique feature

# set working directory:
setwd("//storage.imp.ac.at/groups/massspec/User/Madern/MQLive Setup_remove redundant entries")
list.files()


# read in evidence-file:
evidence <- read.delim(file="evidence.txt", sep= "\t", stringsAsFactors = FALSE, header=TRUE)
dim(evidence)


# create new column: paste "modified sequence" and "charge" together:
evidence$identifier <- paste(evidence$Modified.sequence, evidence$Charge, sep="_charged")


# get a feeling for the data
dim(evidence)
length(unique(evidence$Modified.sequence))
length(unique(evidence$identifier))


# get column of unique entries (corresponds to the number of rows in your output matrix):
unique_evidence <- unique(evidence$identifier)
length_unique_evidence <- length(unique_evidence)
length_unique_evidence


for (i in 1:length(unique_evidence)) {
  
  print(paste0(i, " of ", length(unique_evidence )))
  unique_i <- unique_evidence[i]                               # set unique identifier
  unique_df <- evidence[evidence$identifier %in% unique_i,]    # get all evidences with the same unique identifier
  ret_times_i <- unique_df$Retention.time                      # extract retention times for these evidences
  med_ret_times_i <- median(ret_times_i, na.rm=TRUE)           # calculate median retention time
  intensities_i <- unique_df$Intensity                         # extract intensities for these evidences
  intensities_i[is.na(intensities_i)] <- 0                     # replace NAs with zeros
  max_ind_i <- which.max(intensities_i)                        # which to keep? choose maximum intensity evidence
  keep_i <- unique_df[max_ind_i,]                              # extract the one evidence you want to keep 
  keep_i$Retention.time <- med_ret_times_i                     # replace retention time with median retention time of all evidences found for feature i
  
  
  if (i == 1){                                                 # generate new dataframe with truly unique rows   
    df_res <- keep_i
  } else{
    df_res <- rbind(df_res, keep_i)
  }
}  


# check if it worked:
dim(df_res)
length(unique(df_res$identifier))


# export final dataframe:
write.table(df_res, file="truly_unique_evidences_with_median_retention_time.txt", sep= "\t", row.names = FALSE)

