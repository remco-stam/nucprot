### R script for analysis of proteins with presence/absence between treatments

# Subset 8 and 24 hours
proteinLFQ <- read.table("Step6_proteinLFQs_nonucl.txt", sep = "\t")
df8 <- proteinLFQ[,1:6]
row_nozero8 = apply(df8, 1, function(row) all(row !=0 )) # remove complete cases, they were used in PERSEUS
df8 <-  df8[!row_nozero8,]

df24 <- proteinLFQ[,7:12]
row_nozero24 = apply(df24, 1, function(row) all(row !=0 )) # remove complete cases, they were used in PERSEUS
df24 <-  df24[!row_nozero24,]

## Calculate ratios 8 hr
Rep8_1 <- df8$proteinGroups.nonucl.LFQ.intensity.8hrsI.Rep1/df8$proteinGroups.nonucl.LFQ.intensity.8hrsNI.Rep1
Rep8_2 <- df8$proteinGroups.nonucl.LFQ.intensity.8hrsI.Rep2/df8$proteinGroups.nonucl.LFQ.intensity.8hrsNI.Rep2
Rep8_3 <- df8$proteinGroups.nonucl.LFQ.intensity.8hrsI.Rep3/df8$proteinGroups.nonucl.LFQ.intensity.8hrsNI.Rep3
df8.ratios <- as.data.frame(cbind(Rep8_1,Rep8_2,Rep8_3)) # bind as dataframe
rownames(df8.ratios) <- rownames(df8)
df8.ratios <- df8.ratios[complete.cases(df8.ratios),] # remove all NaNs, these are 0/0 cases so not informative
#Proteins up
df8.twofold <- df8.ratios[df8.ratios[,1]>2 & df8.ratios[,2]>2 & df8.ratios[,3]>2,] # subset all genes with twofold change
row_allInf8 = apply(df8.twofold, 1, function(row) all(row ==Inf )) # remove presence/absence genes number/0 = Inf
df8.poi <- df8.twofold[!row_allInf8,]
#Proteins up
df8.twofoldD <- df8.ratios[df8.ratios[,1]<0.5 & df8.ratios[,2]<0.5 & df8.ratios[,3]<0.5,] # subset all genes with twofold change
row_allInf8D = apply(df8.twofoldD, 1, function(row) all(row ==0 )) # remove presence/absence genes 0/number = 0
df8.poiD <- df8.twofoldD[!row_allInf8D,]

## Calculate ratios 24 hr
Rep24_1 <- df24$proteinGroups.nonucl.LFQ.intensity.24hrsI.Rep1/df24$proteinGroups.nonucl.LFQ.intensity.24hrsNI.Rep1
Rep24_2 <- df24$proteinGroups.nonucl.LFQ.intensity.24hrsI.Rep2/df24$proteinGroups.nonucl.LFQ.intensity.24hrsNI.Rep2
Rep24_3 <- df24$proteinGroups.nonucl.LFQ.intensity.24hrsI.Rep3/df24$proteinGroups.nonucl.LFQ.intensity.24hrsNI.Rep3
df24.ratios <- as.data.frame(cbind(Rep24_1,Rep24_2,Rep24_3)) # bind as dataframe
rownames(df24.ratios) <- rownames(df24)
df24.ratios <- df24.ratios[complete.cases(df24.ratios),] # remove all NaNs, these are 0/0 cases so not informative
#Proteins up
df24.twofold <- df24.ratios[df24.ratios[,1]>2 & df24.ratios[,2]>2 & df24.ratios[,3]>2,] # subset all genes with twofold change
row_allInf24 = apply(df24.twofold, 1, function(row) all(row ==Inf )) # remove presence/absence genes number/0 = Inf
df24.poi <- df24.twofold[!row_allInf24,]
#Proteins down
df24.twofoldD <- df24.ratios[df24.ratios[,1]<0.5 & df24.ratios[,2]<0.5 & df24.ratios[,3]<0.5,] # subset all genes with twofold change
row_allInf24D = apply(df24.twofoldD, 1, function(row) all(row ==0 )) # remove presence/absence genes 0/number = 0
df24.poiD <- df24.twofoldD[!row_allInf24D,]

# Write tables
write.table(df8.poi, "Step10_eightUP.txt", quote = FALSE, sep = "\t", row.names = TRUE)
write.table(df8.poiD, "Step10_eightDOWN.txt", quote = FALSE, sep = "\t", row.names = TRUE)
write.table(df24.poi, "Step10_twfrUP.txt", quote = FALSE, sep = "\t", row.names = TRUE)
write.table(df24.poiD, "Step10_twfrDOWN.txt", quote = FALSE, sep = "\t", row.names = TRUE)
