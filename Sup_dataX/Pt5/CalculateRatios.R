### R script for analysis of proteins with presence/absence between treatments

# Subset I and NI hours
proteinLFQ <- read.table("Step6_proteinLFQs.txt", sep = "\t")
proteinLFQn <- read.table("Step6_proteinLFQs_nonucl.txt", sep = "\t")

dfNI <- proteinLFQ[,c(1,2,3,7,8,9)]
row_nozeroNI = apply(dfNI, 1, function(row) all(row !=0 )) # remove complete cases, they were used in PERSEUS
dfNI <-  dfNI[!row_nozeroNI,]

dfI <- proteinLFQ[,c(4,5,6,10,11,12)]
row_nozeroI = apply(dfI, 1, function(row) all(row !=0 )) # remove complete cases, they were used in PERSEUS
dfI <-  dfI[!row_nozeroI,]

dfNIn <- proteinLFQn[,c(1,2,3,7,8,9)]
row_nozeroNIn = apply(dfNIn, 1, function(row) all(row !=0 )) # remove complete cases, they were used in PERSEUS
dfNIn <-  dfNIn[!row_nozeroNIn,]

dfIn <- proteinLFQn[,c(4,5,6,10,11,12)]
row_nozeroIn = apply(dfIn, 1, function(row) all(row !=0 )) # remove complete cases, they were used in PERSEUS
dfIn <-  dfIn[!row_nozeroIn,]




## Calculate ratios NI
RepNI_1 <- dfNI$proteinGroups.nucl.LFQ.intensity.8hrsNI.Rep1/dfNI$proteinGroups.nucl.LFQ.intensity.24hrsNI.Rep1
RepNI_2 <- dfNI$proteinGroups.nucl.LFQ.intensity.8hrsNI.Rep2/dfNI$proteinGroups.nucl.LFQ.intensity.24hrsNI.Rep2
RepNI_3 <- dfNI$proteinGroups.nucl.LFQ.intensity.8hrsNI.Rep3/dfNI$proteinGroups.nucl.LFQ.intensity.24hrsNI.Rep3
dfNI.ratios <- as.data.frame(cbind(RepNI_1,RepNI_2,RepNI_3)) # bind as dataframe
rownames(dfNI.ratios) <- rownames(dfNI)
dfNI.ratios <- dfNI.ratios[complete.cases(dfNI.ratios),] # remove all NaNs, these are 0/0 cases so not informative
#Proteins up
dfNI.twofold <- dfNI.ratios[dfNI.ratios[,1]>2 & dfNI.ratios[,2]>2 & dfNI.ratios[,3]>2,] # subset all genes with twofold change
row_allInfNI = apply(dfNI.twofold, 1, function(row) all(row ==Inf )) # remove presence/absence genes number/0 = Inf
dfNI.poi <- dfNI.twofold[!row_allInfNI,]
#Proteins up
dfNI.twofoldD <- dfNI.ratios[dfNI.ratios[,1]<0.5 & dfNI.ratios[,2]<0.5 & dfNI.ratios[,3]<0.5,] # subset all genes with twofold change
row_allInfNID = apply(dfNI.twofoldD, 1, function(row) all(row ==0 )) # remove presence/absence genes 0/number = 0
dfNI.poiD <- dfNI.twofoldD[!row_allInfNID,]

## Calculate ratios NI non-nuclear
RepNIn_1 <- dfNIn$proteinGroups.nonucl.LFQ.intensity.8hrsNI.Rep1/dfNIn$proteinGroups.nonucl.LFQ.intensity.24hrsNI.Rep1
RepNIn_2 <- dfNIn$proteinGroups.nonucl.LFQ.intensity.8hrsNI.Rep2/dfNIn$proteinGroups.nonucl.LFQ.intensity.24hrsNI.Rep2
RepNIn_3 <- dfNIn$proteinGroups.nonucl.LFQ.intensity.8hrsNI.Rep3/dfNIn$proteinGroups.nonucl.LFQ.intensity.24hrsNI.Rep3
dfNIn.ratios <- as.data.frame(cbind(RepNIn_1,RepNIn_2,RepNIn_3)) # bind as dataframe
rownames(dfNIn.ratios) <- rownames(dfNIn)
dfNIn.ratios <- dfNIn.ratios[complete.cases(dfNIn.ratios),] # remove all NaNs, these are 0/0 cases so not informative
#Proteins up
dfNIn.twofold <- dfNIn.ratios[dfNIn.ratios[,1]>2 & dfNIn.ratios[,2]>2 & dfNIn.ratios[,3]>2,] # subset all genes with twofold change
row_allInfNIn = apply(dfNIn.twofold, 1, function(row) all(row ==Inf )) # remove presence/absence genes number/0 = Inf
dfNIn.poi <- dfNIn.twofold[!row_allInfNIn,]
#Proteins up
dfNIn.twofoldD <- dfNIn.ratios[dfNIn.ratios[,1]<0.5 & dfNIn.ratios[,2]<0.5 & dfNIn.ratios[,3]<0.5,] # subset all genes with twofold change
row_allInfNInD = apply(dfNIn.twofoldD, 1, function(row) all(row ==0 )) # remove presence/absence genes 0/number = 0
dfNIn.poiD <- dfNIn.twofoldD[!row_allInfNInD,]


## Calculate ratios I
RepI_1 <- dfI$proteinGroups.nucl.LFQ.intensity.8hrsI.Rep1/dfI$proteinGroups.nucl.LFQ.intensity.24hrsI.Rep1
RepI_2 <- dfI$proteinGroups.nucl.LFQ.intensity.8hrsI.Rep2/dfI$proteinGroups.nucl.LFQ.intensity.24hrsI.Rep2
RepI_3 <- dfI$proteinGroups.nucl.LFQ.intensity.8hrsI.Rep3/dfI$proteinGroups.nucl.LFQ.intensity.24hrsI.Rep3
dfI.ratios <- as.data.frame(cbind(RepI_1,RepI_2,RepI_3)) # bind as dataframe
rownames(dfI.ratios) <- rownames(dfI)
dfI.ratios <- dfI.ratios[complete.cases(dfI.ratios),] # remove all NaNs, these are 0/0 cases so not informative
#Proteins up
dfI.twofold <- dfI.ratios[dfI.ratios[,1]>2 & dfI.ratios[,2]>2 & dfI.ratios[,3]>2,] # subset all genes with twofold change
row_allInfI = apply(dfI.twofold, 1, function(row) all(row ==Inf )) # remove presence/absence genes number/0 = Inf
dfI.poi <- dfI.twofold[!row_allInfI,]
#Proteins up
dfI.twofoldD <- dfI.ratios[dfI.ratios[,1]<0.5 & dfI.ratios[,2]<0.5 & dfI.ratios[,3]<0.5,] # subset all genes with twofold change
row_allInfID = apply(dfI.twofoldD, 1, function(row) all(row ==0 )) # remove presence/absence genes 0/number = 0
dfI.poiD <- dfI.twofoldD[!row_allInfID,]

## Calculate ratios I non-nuclear
RepIn_1 <- dfIn$proteinGroups.nonucl.LFQ.intensity.8hrsI.Rep1/dfIn$proteinGroups.nonucl.LFQ.intensity.24hrsI.Rep1
RepIn_2 <- dfIn$proteinGroups.nonucl.LFQ.intensity.8hrsI.Rep2/dfIn$proteinGroups.nonucl.LFQ.intensity.24hrsI.Rep2
RepIn_3 <- dfIn$proteinGroups.nonucl.LFQ.intensity.8hrsI.Rep3/dfIn$proteinGroups.nonucl.LFQ.intensity.24hrsI.Rep3
dfIn.ratios <- as.data.frame(cbind(RepIn_1,RepIn_2,RepIn_3)) # bind as dataframe
rownames(dfIn.ratios) <- rownames(dfIn)
dfIn.ratios <- dfIn.ratios[complete.cases(dfIn.ratios),] # remove all NaNs, these are 0/0 cases so not informative
#Proteins up
dfIn.twofold <- dfIn.ratios[dfIn.ratios[,1]>2 & dfIn.ratios[,2]>2 & dfIn.ratios[,3]>2,] # subset all genes with twofold change
row_allInfIn = apply(dfIn.twofold, 1, function(row) all(row ==Inf )) # remove presence/absence genes number/0 = Inf
dfIn.poi <- dfIn.twofold[!row_allInfIn,]
#Proteins up
dfIn.twofoldD <- dfIn.ratios[dfIn.ratios[,1]<0.5 & dfIn.ratios[,2]<0.5 & dfIn.ratios[,3]<0.5,] # subset all genes with twofold change
row_allInfInD = apply(dfIn.twofoldD, 1, function(row) all(row ==0 )) # remove presence/absence genes 0/number = 0
dfIn.poiD <- dfIn.twofoldD[!row_allInfInD,]

# Write tables
write.table(dfNI.poi, "Step12_NI_UP.txt", quote = FALSE, sep = "\t", row.names = TRUE)
write.table(dfNI.poiD, "Step12_NI_DOWN.txt", quote = FALSE, sep = "\t", row.names = TRUE)
write.table(dfI.poi, "Step12_I_UP.txt", quote = FALSE, sep = "\t", row.names = TRUE)
write.table(dfI.poiD, "Step12_I_DOWN.txt", quote = FALSE, sep = "\t", row.names = TRUE)
write.table(dfNIn.poi, "Step12_NI_UP_nonucl.txt", quote = FALSE, sep = "\t", row.names = TRUE)
write.table(dfNIn.poiD, "Step12_NI_DOWN_nonucl.txt", quote = FALSE, sep = "\t", row.names = TRUE)
write.table(dfIn.poi, "Step12_I_UP_nonucl.txt", quote = FALSE, sep = "\t", row.names = TRUE)
write.table(dfIn.poiD, "Step12_I_DOWN_nonucl.txt", quote = FALSE, sep = "\t", row.names = TRUE)
