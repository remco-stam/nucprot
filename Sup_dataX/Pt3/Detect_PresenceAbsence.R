### R script for analysis of proteins with presence/absence between treatments

## Split time points for individual analysis
# Read data file and subset
proteinLFQ <- read.table("Step6_proteinLFQs.txt", sep = "\t")
df8 <- proteinLFQ[,1:6]
df24 <- proteinLFQ[,7:12]

## Find presence/absence. 
eightDOWN <- df8[df8[,1]>0 & df8[,2] >0 & df8[,3]>0 & df8[,4]==0 & df8[,5]==0 & df8[,6]==0,]
eightUP <- df8[df8[,1]==0 & df8[,2] ==0 & df8[,3]==0 & df8[,4]>0 & df8[,5]>0 & df8[,6]>0,]

twfrDOWN <- df24[df24[,1]>0 & df24[,2] >0 & df24[,3]>0 & df24[,4]==0 & df24[,5]==0 & df24[,6]==0,]
twfrUP <- df24[df24[,1]==0 & df24[,2] ==0 & df24[,3]==0 & df24[,4]>0 & df24[,5]>0 & df24[,6]>0,]

#Write tables
write.table(eightDOWN, "Step7_eightDOWN.txt", quote = FALSE, sep = "\t")
write.table(eightUP, "Step7_eightUP.txt", quote = FALSE, sep = "\t")
write.table(twfrDOWN, "Step7_twfrDOWN.txt", quote = FALSE, sep = "\t")
write.table(twfrUP, "Step7_twfrUP.txt", quote = FALSE, sep = "\t")
