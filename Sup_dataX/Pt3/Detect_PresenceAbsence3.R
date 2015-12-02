### R script for analysis of proteins with presence/absence between treatments

## Split time points for individual analysis
# Read data file and subset
proteinLFQ <- read.table("Step6_proteinLFQs2.txt", sep = "\t", row.names =1)
colnames(proteinLFQ) <- c("24I1","24I2","24I3","24NI1","24NI2","24NI3","8I1","8I2","8I3","8NI1","8NI2","8NI3")
proteinLFQ <- as.matrix(proteinLFQ[2:3572,])
df24 <- proteinLFQ[,1:6]
df8 <- proteinLFQ[,7:12]
dfI <- proteinLFQ[,c(1,2,3,7,8,9)]
dfNI <- proteinLFQ[,c(4,5,6,10,11,12)]


## Find presence/absence. 
eightUP <- df8[df8[,1]>0 & df8[,2] >0 & df8[,3]>0 & df8[,4]==0 & df8[,5]==0 & df8[,6]==0,]
eightDOWN <- df8[df8[,1]==0 & df8[,2] ==0 & df8[,3]==0 & df8[,4]>0 & df8[,5]>0 & df8[,6]>0,]

twfrUP <- df24[df24[,1]>0 & df24[,2] >0 & df24[,3]>0 & df24[,4]==0 & df24[,5]==0 & df24[,6]==0,]
twfrDOWN <- df24[df24[,1]==0 & df24[,2] ==0 & df24[,3]==0 & df24[,4]>0 & df24[,5]>0 & df24[,6]>0,]

NI_UP <- dfNI[dfNI[,1]>0 & dfNI[,2] >0 & dfNI[,3]>0 & dfNI[,4]==0 & dfNI[,5]==0 & dfNI[,6]==0,]
NI_DOWN <- dfNI[dfNI[,1]==0 & dfNI[,2] ==0 & dfNI[,3]==0 & dfNI[,4]>0 & dfNI[,5]>0 & dfNI[,6]>0,]

I_UP <- dfI[dfI[,1]>0 & dfI[,2] >0 & dfI[,3]>0 & dfI[,4]==0 & dfI[,5]==0 & dfI[,6]==0,]
I_DOWN <- dfI[dfI[,1]==0 & dfI[,2] ==0 & dfI[,3]==0 & dfI[,4]>0 & dfI[,5]>0 & dfI[,6]>0,]


#Write tables
write.table(eightDOWN, "eightDOWN.txt", quote = FALSE, sep = "\t")
write.table(eightUP, "eightUP.txt", quote = FALSE, sep = "\t")
write.table(twfrDOWN, "twfrDOWN.txt", quote = FALSE, sep = "\t")
write.table(twfrUP, "twfrUP.txt", quote = FALSE, sep = "\t")
write.table(I_DOWN, "IDOWN.txt", quote = FALSE, sep = "\t")
write.table(I_UP, "IUP.txt", quote = FALSE, sep = "\t")
write.table(NI_DOWN, "NIDOWN.txt", quote = FALSE, sep = "\t")
write.table(NI_UP, "NIUP.txt", quote = FALSE, sep = "\t")