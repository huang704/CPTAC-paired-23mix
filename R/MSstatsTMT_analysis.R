library(MSstatsTMT) # version 2.2.7


######################################################################
# estimate protein abundances
MSstatsTMT_annotation <- read.csv("MSstatsTMT_annotation.csv")
MSstatsTMT_annotation$Run <- paste0(MSstatsTMT_annotation$Run, ".raw")

processed_psm <- PhilosophertoMSstatsTMTFormat(path = "Philosopher", 
                                               annotation = MSstatsTMT_annotation, 
                                               protein_id_col = "ProteinAccessions", 
                                               peptide_id_col = "PeptideSequence", Purity_cutoff = 0.6,
                                               PeptideProphet_prob_cutoff = 0.7, useUniquePeptide = TRUE,
                                               rmPSM_withfewMea_withinRun = TRUE, rmPeptide_OxidationM = TRUE,
                                               rmProtein_with1Feature = FALSE, summaryforMultipleRows = sum)

save(processed_psm, file = "processed_psm.rda")

## count the number of proteins
length(unique(processed_psm$ProteinName))  # 11984
length(unique(processed_psm$PSM))  # 321110
length(unique(processed_psm$PeptideSequence))  # 214864

# use MSstats for protein summarization
MSstatsTMT.abun.norm <- proteinSummarization(processed_psm,
                                             method="msstats",
                                             global_norm=TRUE,
                                             reference_norm=TRUE)

head(MSstatsTMT.abun.norm$ProteinLevelData)
save(MSstatsTMT.abun.norm, file='MSstatsTMT.abun.norm.rda')
unique(MSstatsTMT.abun.norm$ProteinLevelData$Run)

######################################################################
# variance components
library(tidyverse)
library(lmerTest)
library(lme4)
load("MSstatsTMT.abun.norm.rda")
data <- MSstatsTMT.abun.norm$ProteinLevelData
data <- data %>% filter(!Condition %in% c("Norm", "Empty") & !Mixture %in% c("19", "20")
                        & !grepl("contam_sp", Protein))
length(unique(data$Protein)) # 11894

filtered_data <- data %>% filter(!is.na(Abundance)) %>%
  dplyr::group_by(Protein) %>%
  dplyr::summarise(n = n()) %>%
  right_join(data) %>%
  filter(n == max(n))

proteins <- as.character(unique(data$Protein)) #7215
varscomp <- matrix(rep(0, 2*length(proteins)), ncol = 2)
terms <- c("BioReplicate", "Residual")
for(i in 1:length(proteins)) {
  message("Protein: ", i)
  sub_data <- data[data$Protein == proteins[i],]
  ## linear mixed model
  fit.mixed <- try(lmer(Abundance ~ 1 + Condition + (1|BioReplicate), data=sub_data), TRUE)
  if(!inherits(fit.mixed, "try-error")){
    vc <- as.data.frame(VarCorr(fit.mixed, comp="Variance"))
    rownames(vc) <- vc$grp
    varscomp[i,] <- vc[terms, "vcov"]
  } else{
    varscomp[i,] <- NA
  }
}
colnames(varscomp) <- terms
rownames(varscomp) <- proteins
save(varscomp, file = 'CS-variance-comp.rda')

mean(varscomp[,"BioReplicate"], na.rm=TRUE)

load('CS-variance-comp.rda')
library(plyr)
plot.data.raw <-na.omit(varscomp)
plot.data.raw <- as.data.frame(plot.data.raw) # 2 proteins removed
plot.data.raw$ProteinName <- rownames(plot.data.raw)
plot.data <- plot.data.raw %>% gather(Source, variance, -ProteinName)
perc.plot.data = ddply(plot.data, "ProteinName", mutate, percent_variance = variance/sum(variance) * 100)

perc.plot.data%>%
  group_by(Source)%>%
  dplyr::summarise(mean = mean(percent_variance), median = median(percent_variance))

perc.plot.data%>%
  group_by(Source)%>%
  dplyr::summarise(mean = mean(variance), median = median(variance))
