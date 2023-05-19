# calculate enrichment based on evidence.txt

# read in evidence.txt file
# you can use choose.files()
evidence_file <- "C:\\RAWFILES\\Master\\combined\\txt\\evidence.txt"
evidence <- read.delim(file=evidence_file,sep="\t", header=TRUE,stringsAsFactors=FALSE)

# get rid of CONs, and reverse
evidence <- evidence[!evidence$Reverse=="+",]
evidence <- evidence[!evidence$Potential.contaminant=="+",]
evidence <- evidence[!grepl(evidence$Proteins, pattern="iRT"),]

# get overview of raw files in evidence.txt
table(evidence$Raw.file)
length(table(evidence$Raw.file))

# only keep acetylome raw-files (only relevant when searched together with proteome files)
bool_acetylom <- grepl(evidence$Raw.file, pattern="acet")
evidence <- evidence[bool_acetylom,]
evidence$Modifications


# calculate enrichment
enrichment<- numeric(length(unique(evidence$Raw.file)))   # initialize
names(enrichment) <- unique(evidence$Raw.file)
number_acetyl_evidences <- numeric(length(unique(evidence$Raw.file)))
names(number_acetyl_evidences) <- unique(evidence$Raw.file)


for (i in 1:length(enrichment)){

  rawfile_i <- unique(evidence$Raw.file)[i]
  evidence_i <- evidence[evidence$Raw.file == rawfile_i,]
    
  modifications_i <- evidence_i$Modifications
  bool_acetylated_i<- grepl(modifications_i, pattern="Acetyl [(]K[)]")    #check if this pattern matches how acetyl-k peptides are written in the evidence.txt!
    
  enrichment_i <- sum(bool_acetylated_i)/length(bool_acetylated_i)
  number_acetyl_evidences[i] <-  sum(bool_acetylated_i)  
  enrichment[i]<- enrichment_i*100
  
}

number_acetyl_evidences
enrichment



