### R script for modification of proteinGroups.txt table

## Create file with nuclear proteins for PERSEUS analysis 
# Read datafiles, proteinGroups file (with added ID) and IDs to be selected
proteinGroups <- read.delim("proteinGroups_ID.txt", header = TRUE)
proteinIDs.tom<- read.delim("Step1_Proteomics_List.text", header = TRUE)
proteinIDs.nucl <- read.table("Step4_Nuclear_IDs.txt", header = FALSE)

# Subset tomato proteins only, remove Phytophthora and contaminants
proteinGroups.tom <- subset(proteinGroups, proteinGroups[,1] %in% proteinIDs.tom[,1])

# Subset and write tables.
proteinGroups.nucl <- subset(proteinGroups.tom, proteinGroups.tom[,1] %in% proteinIDs.nucl[,1])
proteinGroups.nonucl <- subset(proteinGroups.tom, !proteinGroups.tom[,1] %in% proteinIDs.nucl[,1])
write.table(proteinGroups.nucl, "Step6_proteinGroups_ID_nucl.txt", quote = FALSE, sep = "\t")
write.table(proteinGroups.nonucl, "Step6_proteinGroups_ID_notnucl.txt", quote = FALSE, sep = "\t")

# For PERSEUS, subset tables without 0 values
proteinGroups.LFQ8 <- proteinGroups.nucl[,111:116]
row_nozero8 = apply(proteinGroups.LFQ8, 1, function(row) all(row !=0 ))
proteinGroups.nozero8 <-  proteinGroups.nucl[row_nozero8,]
#write.table(proteinGroups.nozero8, "Step6_proteinGroups_nozero8.txt", quote = FALSE, sep = "\t")

proteinGroups.LFQ24 <- proteinGroups.nucl[,105:110]
row_nozero24 = apply(proteinGroups.LFQ24, 1, function(row) all(row !=0 ))
proteinGroups.nozero24 <-  proteinGroups.nucl[row_nozero24,]
#write.table(proteinGroups.nozero24, "Step6_proteinGroups_nozero24.txt", quote = FALSE, sep = "\t")

## Create files with intensity values for other analysis
# Subset LFQ values, other columns are not needed for this analysis
proteinLFQ <- as.data.frame(proteinGroups.nucl[,1])

proteinLFQ[,2:13] <- as.numeric(c(proteinGroups.nucl$LFQ.intensity.8hrsNI.Rep1, proteinGroups.nucl$LFQ.intensity.8hrsNI.Rep2, 
                      proteinGroups.nucl$LFQ.intensity.8hrsNI.Rep3, proteinGroups.nucl$LFQ.intensity.8hrsI.Rep1, 
                      proteinGroups.nucl$LFQ.intensity.8hrsI.Rep2, proteinGroups.nucl$LFQ.intensity.8hrsI.Rep3,
                      proteinGroups.nucl$LFQ.intensity.24hrsNI.Rep1, proteinGroups.nucl$LFQ.intensity.24hrsNI.Rep2, 
                      proteinGroups.nucl$LFQ.intensity.24hrsNI.Rep3, proteinGroups.nucl$LFQ.intensity.24hrsI.Rep1, 
                      proteinGroups.nucl$LFQ.intensity.24hrsI.Rep2, proteinGroups.nucl$LFQ.intensity.24hrsI.Rep3))
colnames(proteinLFQ) <- c( "Protein", 
                        "proteinGroups.nucl$LFQ.intensity.8hrsNI.Rep1", "proteinGroups.nucl$LFQ.intensity.8hrsNI.Rep2", 
                      "proteinGroups.nucl$LFQ.intensity.8hrsNI.Rep3", "proteinGroups.nucl$LFQ.intensity.8hrsI.Rep1", 
                      "proteinGroups.nucl$LFQ.intensity.8hrsI.Rep2", "proteinGroups.nucl$LFQ.intensity.8hrsI.Rep3",
                      "proteinGroups.nucl$LFQ.intensity.24hrsNI.Rep1", "proteinGroups.nucl$LFQ.intensity.24hrsNI.Rep2", 
                      "proteinGroups.nucl$LFQ.intensity.24hrsNI.Rep3", "proteinGroups.nucl$LFQ.intensity.24hrsI.Rep1", 
                      "proteinGroups.nucl$LFQ.intensity.24hrsI.Rep2", "proteinGroups.nucl$LFQ.intensity.24hrsI.Rep3")
rownames(proteinLFQ) <- proteinLFQ[,1]
proteinLFQ[,1] <-NULL

#write.table(proteinLFQ, "Step6_proteinLFQs.txt", quote = FALSE, sep = "\t")


### Repeat with proteins predicted to be not-nuclear


# For PERSEUS, subset tables without 0 values
proteinGroups.LFQ8n <- proteinGroups.nonucl[,111:116]
row_nozero8n = apply(proteinGroups.LFQ8n, 1, function(row) all(row !=0 ))
proteinGroups.nozero8n <-  proteinGroups.nonucl[row_nozero8n,]
write.table(proteinGroups.nozero8n, "Step6_proteinGroups_nozero8_nonucl.txt", quote = FALSE, sep = "\t")

proteinGroups.LFQ24n <- proteinGroups.nonucl[,105:110]
row_nozero24n = apply(proteinGroups.LFQ24n, 1, function(row) all(row !=0 ))
proteinGroups.nozero24n <-  proteinGroups.nonucl[row_nozero24n,]
write.table(proteinGroups.nozero24n, "Step6_proteinGroups_nozero24_nonucl.txt", quote = FALSE, sep = "\t")

## Create files with intensity values for other analysis
# Subset LFQ values, other columns are not needed for this analysis
proteinLFQ_nonucl <- as.data.frame(proteinGroups.nonucl[,1])

proteinLFQ_nonucl[,2:13] <- as.numeric(c(proteinGroups.nonucl$LFQ.intensity.8hrsNI.Rep1, proteinGroups.nonucl$LFQ.intensity.8hrsNI.Rep2, 
                                         proteinGroups.nonucl$LFQ.intensity.8hrsNI.Rep3, proteinGroups.nonucl$LFQ.intensity.8hrsI.Rep1, 
                                         proteinGroups.nonucl$LFQ.intensity.8hrsI.Rep2, proteinGroups.nonucl$LFQ.intensity.8hrsI.Rep3,
                                         proteinGroups.nonucl$LFQ.intensity.24hrsNI.Rep1, proteinGroups.nonucl$LFQ.intensity.24hrsNI.Rep2, 
                                         proteinGroups.nonucl$LFQ.intensity.24hrsNI.Rep3, proteinGroups.nonucl$LFQ.intensity.24hrsI.Rep1, 
                                         proteinGroups.nonucl$LFQ.intensity.24hrsI.Rep2, proteinGroups.nonucl$LFQ.intensity.24hrsI.Rep3))
colnames(proteinLFQ_nonucl) <- c( "Protein", 
                           "proteinGroups.nonucl$LFQ.intensity.8hrsNI.Rep1", "proteinGroups.nonucl$LFQ.intensity.8hrsNI.Rep2", 
                           "proteinGroups.nonucl$LFQ.intensity.8hrsNI.Rep3", "proteinGroups.nonucl$LFQ.intensity.8hrsI.Rep1", 
                           "proteinGroups.nonucl$LFQ.intensity.8hrsI.Rep2", "proteinGroups.nonucl$LFQ.intensity.8hrsI.Rep3",
                           "proteinGroups.nonucl$LFQ.intensity.24hrsNI.Rep1", "proteinGroups.nonucl$LFQ.intensity.24hrsNI.Rep2", 
                           "proteinGroups.nonucl$LFQ.intensity.24hrsNI.Rep3", "proteinGroups.nonucl$LFQ.intensity.24hrsI.Rep1", 
                           "proteinGroups.nonucl$LFQ.intensity.24hrsI.Rep2", "proteinGroups.nonucl$LFQ.intensity.24hrsI.Rep3")
rownames(proteinLFQ_nonucl) <- proteinLFQ_nonucl[,1]
proteinLFQ_nonucl[,1] <-NULL

write.table(proteinLFQ_nonucl, "Step6_proteinLFQs_nonucl.txt", quote = FALSE, sep = "\t")
