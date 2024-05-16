library(readr)
library(dplyr)
library(decoupleR)
library(reshape2)

GDSC_rnaseq_tpm_20220624 <- as.data.frame(
  read_csv("data/rnaseq_tpm_20220624.csv"))[,-1]
names(GDSC_rnaseq_tpm_20220624)[1] <- "GENE_SYMBOLS"

GDSC_rnaseq_tpm_20220624 <- GDSC_rnaseq_tpm_20220624 %>% group_by(GENE_SYMBOLS) %>% summarise_each(mean)
GDSC_rnaseq_tpm_20220624 <- as.data.frame(GDSC_rnaseq_tpm_20220624)
GDSC_rnaseq_tpm_20220624 <- GDSC_rnaseq_tpm_20220624[!is.na(GDSC_rnaseq_tpm_20220624$GENE_SYMBOLS),]

row.names(GDSC_rnaseq_tpm_20220624) <- GDSC_rnaseq_tpm_20220624$GENE_SYMBOLS

GDSC_rnaseq_tpm_20220624 <- GDSC_rnaseq_tpm_20220624[,-1]
GDSC_rnaseq_tpm_20220624 <- GDSC_rnaseq_tpm_20220624[complete.cases(GDSC_rnaseq_tpm_20220624),]

#to match the CCLE dataset
GDSC_rnaseq_tpm_20220624_logp1 <- log2(GDSC_rnaseq_tpm_20220624+1)

#extremelly basic feature selection to have a reasonable sized model
plot(density(rowSums(GDSC_rnaseq_tpm_20220624_logp1)))
GDSC_rnaseq_tpm_20220624_logp1 <- GDSC_rnaseq_tpm_20220624_logp1[rowSums(GDSC_rnaseq_tpm_20220624_logp1) > 3500,]


to_write <- as.data.frame(cbind(row.names(GDSC_rnaseq_tpm_20220624_logp1),GDSC_rnaseq_tpm_20220624_logp1))
names(to_write)[1] <- "GENE_SYMBOLS"


write_csv(to_write, file = "data/rnaseq_tpm_20220624_logp1.csv")


### collectrI
# collectrI <- decoupleR::get_collectri()
# write_csv(collectrI, file = "support/collectrI_052024.csv")
collectrI <- as.data.frame(read_csv("support/collectrI_052024.csv"))

TF_activities <- run_ulm(GDSC_rnaseq_tpm_20220624_logp1, collectrI)

TF_activities_df <- reshape2::dcast(TF_activities[,c(2,3,4)], source~condition)

write_csv(TF_activities_df, file = "results/GDSC_TF_activities.csv")
