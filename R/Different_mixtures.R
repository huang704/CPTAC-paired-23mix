library(dplyr)
library(tidyr)
library(MSstatsTMT)
library(invgamma)

setwd("~/Dropbox/Ting-thesis/Complex_design/results/TMT-CCRCC/github")
source('R/modelEvaluation.R')
source("R/groupComparisonTMT.R")
source("R/utils_group_comparison.R")

load("MSstatsTMT.abun.norm.rda")
data <- MSstatsTMT.abun.norm$ProteinLevelData
data <- data %>% filter(!Condition %in% c("Norm", "Empty") & !Mixture %in% c("19", "20")
                        & !grepl("contam_sp", Protein))

length(unique(data$Protein))

########################################################################
########################################################################
# multiple mixtures
fx1 <- formula("Abundance ~ 1 + Group + (1|Mixture) + (1|Mixture:Subject)")
fx2 <- formula("Abundance ~ 1 + Group + (1|Subject)")
fx3 <- formula("Abundance ~ 1 + Group + (1|Mixture)")
fx4 <- formula("Abundance ~ 1 + Group")

methods <- c("CMS", "CS", "CM", "C")
methods_2 <- c("CS", "C")
model_list <- c(fx1, fx2, fx3, fx4)
model_list_2 <- c(fx2, fx4)
num_mixtures <- c(1, 3, 5, 10, 25)
# statistical testing
for(j in seq_along(num_mixtures)){
  set.seed(0623+j)
  
  selected_mixtures <- sample(unique(data$Mixture), num_mixtures[j])
  
  if(num_mixtures[j] == 1){ # single mixture
    for(k in seq_along(model_list_2)){
      # without moderation
      data.res <- groupComparisonTMT.v2(data = data %>% filter(Mixture %in% selected_mixtures),
                                        formula = model_list_2[[k]],
                                        contrast.matrix = "pairwise",
                                        moderated = FALSE)
      filename1 <- paste0('diff-mixtures/', num_mixtures[j], 'mixtures-by-moderated', methods_2[k], ".testing.rda")
      save(data.res, file=filename1)
      rm(data.res)
      rm(filename1)
      
      # with moderation
      data.res <- groupComparisonTMT.v2(data = data %>% filter(Mixture %in% selected_mixtures),
                                        formula = model_list_2[[k]],
                                        contrast.matrix = "pairwise",
                                        moderated = TRUE)
      filename2 <- paste0('diff-mixtures/', num_mixtures[j], 'mixtures-by-moderated', methods_2[k], ".testing-EB.rda")
      save(data.res, file=filename2)
      rm(data.res)
      rm(filename2)
    }
  } else{ # multiple mixtures
    for(k in seq_along(model_list)){
      # without moderation
      data.res <- groupComparisonTMT.v2(data = data %>% filter(Mixture %in% selected_mixtures),
                                        formula = model_list[[k]],
                                        contrast.matrix = "pairwise",
                                        moderated = FALSE)
      filename1 <- paste0('diff-mixtures/', num_mixtures[j], 'mixtures-by-moderated', methods[k], ".testing.rda")
      save(data.res, file=filename1)
      rm(data.res)
      rm(filename1)
      
      # with moderation
      data.res <- groupComparisonTMT.v2(data = data %>% filter(Mixture %in% selected_mixtures),
                                        formula = model_list[[k]],
                                        contrast.matrix = "pairwise",
                                        moderated = TRUE)
      
      filename2 <- paste0('diff-mixtures/', num_mixtures[j], 'mixtures-by-moderated', methods[k], ".testing-EB.rda")
      save(data.res, file=filename2)
      rm(data.res)
      rm(filename2)
    }
    
    # run two-way limma
    data.res <- groupComparison(data = data %>% filter(Mixture %in% selected_mixtures),
                                method = "limma+timecourse+twoway",
                                contrast.matrix = "pairwise",
                                moderated = TRUE, 
                                adj.method = "BH")
    filename4 <- paste0('diff-mixtures/', num_mixtures[j], 'mixtures-by-limmaCMS.testing.rda')
    save(data.res, file=filename4)
    rm(data.res)
    rm(filename4)
  }
  
  data.res <- groupComparison(data = data %>% filter(Mixture %in% selected_mixtures),
                              method = "limma+timecourse",
                              contrast.matrix = "pairwise",
                              moderated = TRUE, 
                              adj.method = "BH")
  filename3 <- paste0('diff-mixtures/', num_mixtures[j], 'mixtures-by-limmaCS.testing.rda')
  save(data.res, file=filename3)
  rm(data.res)
  rm(filename3)
}
