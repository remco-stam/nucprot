### R script for analysis of proteins with presence/absence between treatments

## Split time points for individual analysis
# Read data file and subset
proteinLFQ <- read.table("Step6_proteinLFQs.txt", sep = "\t")
proteinLFQn <- read.table("Step6_proteinLFQs_nonucl.txt", sep = "\t")

dfNI <- proteinLFQ[,c(1,2,3,7,8,9)]
dfI <- proteinLFQ[,c(4,5,6,10,11,12)]

dfNIn <- proteinLFQn[,c(1,2,3,7,8,9)]
dfIn <- proteinLFQn[,c(4,5,6,10,11,12)]

## Find presence/absence. (nuclear proteins)
NI_DOWN <- dfNI[dfNI[,1]>0 & dfNI[,2] >0 & dfNI[,3]>0 & dfNI[,4]==0 & dfNI[,5]==0 & dfNI[,6]==0,]
NI_UP <- dfNI[dfNI[,1]==0 & dfNI[,2] ==0 & dfNI[,3]==0 & dfNI[,4]>0 & dfNI[,5]>0 & dfNI[,6]>0,]

I_DOWN <- dfI[dfI[,1]>0 & dfI[,2] >0 & dfI[,3]>0 & dfI[,4]==0 & dfI[,5]==0 & dfI[,6]==0,]
I_UP <- dfI[dfI[,1]==0 & dfI[,2] ==0 & dfI[,3]==0 & dfI[,4]>0 & dfI[,5]>0 & dfI[,6]>0,]

#Write tables
write.table(NI_DOWN, "Step11_NI_DOWN.txt", quote = FALSE, sep = "\t")
write.table(NI_UP, "Step11_NI_UP.txt", quote = FALSE, sep = "\t")
write.table(I_DOWN, "Step11_I_DOWN.txt", quote = FALSE, sep = "\t")
write.table(I_UP, "Step11_I_UP.txt", quote = FALSE, sep = "\t")



## Find presence/absence. (notnuclear proteins)
NIn_DOWN <- dfNIn[dfNIn[,1]>0 & dfNIn[,2] >0 & dfNIn[,3]>0 & dfNIn[,4]==0 & dfNIn[,5]==0 & dfNIn[,6]==0,]
NIn_UP <- dfNIn[dfNIn[,1]==0 & dfNIn[,2] ==0 & dfNIn[,3]==0 & dfNIn[,4]>0 & dfNIn[,5]>0 & dfNIn[,6]>0,]

In_DOWN <- dfIn[dfIn[,1]>0 & dfIn[,2] >0 & dfIn[,3]>0 & dfIn[,4]==0 & dfIn[,5]==0 & dfIn[,6]==0,]
In_UP <- dfIn[dfIn[,1]==0 & dfIn[,2] ==0 & dfIn[,3]==0 & dfIn[,4]>0 & dfIn[,5]>0 & dfIn[,6]>0,]

#Write tables
write.table(NIn_DOWN, "Step11_NI_DOWN_nonucl.txt", quote = FALSE, sep = "\t")
write.table(NIn_UP, "Step11_NI_UP_nonucl.txt", quote = FALSE, sep = "\t")
write.table(In_DOWN, "Step11_I_DOWN_nonucl.txt", quote = FALSE, sep = "\t")
write.table(In_UP, "Step11_I_UP_nonucl.txt", quote = FALSE, sep = "\t")

