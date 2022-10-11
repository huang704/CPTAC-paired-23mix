library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pROC)
library(RColorBrewer)
library(MSstats)

num_mixtures <- c(1, 3, 5, 10, 21)
methods <- c("moderatedCMS", "moderatedCS", "moderatedCM", "moderatedC", "limmaCS", "limmaCMS")
methods_2 <- c("moderatedCS", "moderatedC", "limmaCS")

new_methods <- c("moderatedTMS", "moderatedTS", "moderatedTM", "moderatedT", "limmaTS", "limmaTMS")
new_methods_2 <- c("moderatedTS", "moderatedT", "limmaTS")

res <- list()
count = 0

for(j in seq_along(num_mixtures)){
  if(num_mixtures[j] == 1){
    for(k in seq_along(methods_2)){
      count <- count + 1
      
      if(k <= 2){
        load(paste0('diff-mixtures/', num_mixtures[j], 'mixtures-by-', methods_2[k], ".testing-EB.rda"))
        data.res <-data.res$ComparisonResult
      }
      else{
        load(paste0('diff-mixtures/', num_mixtures[j], 'mixtures-by-', methods_2[k], ".testing.rda"))
      }
      
      data.res <- data.res[, c("Protein", "Label", "log2FC", "SE", "DF", "pvalue",  "adj.pvalue")]
      data.res <- data.res %>% filter(!is.na(adj.pvalue))
      data.res <- as.data.frame(data.res)
      data.res$Label <- as.character(data.res$Label)
      data.res$Protein <- as.character(data.res$Protein)
      data.res$Method <- new_methods_2[k]
      data.res$Num_mixtures <- num_mixtures[j]
      res[[count]] <- data.res
    }
  } else{
    for(k in seq_along(methods)){
      count <- count + 1
      
      if(k <= 4){
        load(paste0('diff-mixtures/', num_mixtures[j], 'mixtures-by-', methods[k], ".testing-EB.rda"))
        data.res <-data.res$ComparisonResult
      }
      else{
        load(paste0('diff-mixtures/', num_mixtures[j], 'mixtures-by-', methods[k], ".testing.rda"))
      }
      
      data.res <- data.res[, c("Protein", "Label", "log2FC", "SE", "DF", "pvalue",  "adj.pvalue")]
      data.res <- data.res %>% filter(!is.na(adj.pvalue))
      data.res <- as.data.frame(data.res)
      data.res$Label <- as.character(data.res$Label)
      data.res$Protein <- as.character(data.res$Protein)
      data.res$Method <- new_methods[k]
      data.res$Num_mixtures <- num_mixtures[j]
      res[[count]] <- data.res
    }
  }
  
}


res1 <- data.table::rbindlist(res)
res1$log2FC <- as.numeric(as.character(res1$log2FC))
res1$adj.pvalue <- as.numeric(as.character(res1$adj.pvalue))

res1 <- within(res1, Method <- factor(Method, levels = new_methods))

res1$pred <- ifelse(res1$adj.pvalue <= 0.05, 1, 0)
res1[res1$Label == "Tumor-Normal", "Label"] <- "Tumor vs Normal"
res1[res1$Label == "Tumor-Normal", "Label"] <- "Tumor vs Normal"
res1[res1$Label == "Normal-Tumor", "log2FC"] <- -res1[res1$Label == "Normal-Tumor", "log2FC"]
res1[res1$Label == "Normal-Tumor", "Label"] <- "Tumor vs Normal"
res1[res1$Label == "Normal vs Tumor", "log2FC"] <- -res1[res1$Label == "Normal vs Tumor", "log2FC"]
res1[res1$Label == "Normal vs Tumor", "Label"] <- "Tumor vs Normal"
unique(res1$Label)
save(res1, file = "diff-mixtures/Results_diff_mixtures.rda")

res1 <- res1 %>% filter(Method != "moderatedTMS")
num_testable_prots <- res1 %>% group_by(Num_mixtures, Method) %>%
  dplyr::summarise(n= n_distinct(Protein)) %>% 
  spread(Num_mixtures, n)
write.csv(num_testable_prots, file = "num_testable_prots_multi_mixtures.csv")

sig <- res1 %>% filter(adj.pvalue <= 0.05)
perf <- sig %>% group_by(Num_mixtures, Method) %>%
  dplyr::summarise(n= n_distinct(Protein)) %>% 
  spread(Num_mixtures, n)
perf
write.csv(perf, file = "diff-mixtures/performance_multi_mixtures.csv")

blues_colors <- brewer.pal(n = 9, name = "Blues")
purple_colors <- brewer.pal(n = 9, name = "PuRd")
colors <- c(purple_colors[9], blues_colors[6], blues_colors[8], purple_colors[7], purple_colors[4])

pdf("Compare_different_mixtures.pdf",width=10,height=3)
sig %>% group_by(Num_mixtures, Method) %>%
  dplyr::summarise(n= n_distinct(Protein)) %>% 
  spread(Num_mixtures, n) %>% 
  gather(Num_mixtures, n, -Method) %>%
  filter(Num_mixtures != 1) %>%
  ggplot(aes(x=factor(Num_mixtures, levels = c(3, 5, 10, 21)), y=n, fill=Method)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_text(aes(label=n), vjust=1.6, color="white", position = position_dodge(0.9), size=3)+
  scale_fill_manual(values=colors) +
  labs(x="# of mixtures", y = "# of significant proteins", size=12)+ theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(angle=0, size=12), #, vjust = 0.43
        axis.text.x=element_text(angle=0, size=12),
        axis.title=element_text(size=12,face="bold"),
        strip.text.x = element_text(size=12), # angle=75
        strip.text.y = element_text(size=12),
        strip.background = element_rect(colour=NA, fill=NA),
        legend.position = 'right',
        legend.title = element_text(size=11, face="bold"),
        legend.text = element_text(size=11))
dev.off()

########################################################################
# Distribution of SE
blues_colors <- brewer.pal(n = 9, name = "Blues")
purple_colors <- brewer.pal(n = 9, name = "PuRd")
colors <- c(purple_colors[9], blues_colors[6], blues_colors[8], purple_colors[7], purple_colors[4])

pdf("Compare_different_methods_fc.pdf",width=6,height=3.5)
g4 <- res1 %>% filter(Method != "moderatedTMS" & Num_mixtures != 1) %>% 
  ggplot(aes(x=as.factor(Num_mixtures), y= log2FC, color= Method))+
  #geom_jitter(position=position_jitter(0.1))+ # v1 : width=0.5, v2 n v3 : width=1
  geom_boxplot(lwd=0.6, outlier.alpha = 0.02, notch = TRUE)+ 
  scale_color_manual(values=colors) + 
  ylim(-2,2)+
  labs(x=("# of mixtures"), y = ("Estimated log2-fold change"), size=12)+ 
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(angle=0, size=12), #, vjust = 0.43
        axis.text.x=element_text(angle=0, size=12),
        axis.title=element_text(size=12,face="bold"),
        strip.text.x = element_text(size=12), # angle=75
        strip.text.y = element_text(size=12),
        strip.background = element_rect(colour=NA, fill=NA),
        legend.position = 'bottom')

g4
dev.off()


pdf("Compare_different_methods_se.pdf",width=6,height=3.5)
g4 <- res1 %>% filter(Method != "moderatedTMS" & Num_mixtures != 1) %>% 
  ggplot(aes(x=as.factor(Num_mixtures), y= SE, color= Method))+
  #geom_jitter(position=position_jitter(0.1))+ # v1 : width=0.5, v2 n v3 : width=1
  geom_boxplot(lwd=0.6, outlier.alpha = 0.02, notch = TRUE)+ 
  scale_color_manual(values=colors) + 
  ylim(0,0.3)+
  labs(x=("# of mixtures"), y = ("Standard error"), size=12)+ 
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(angle=0, size=12), #, vjust = 0.43
        axis.text.x=element_text(angle=0, size=12),
        axis.title=element_text(size=12,face="bold"),
        strip.text.x = element_text(size=12), # angle=75
        strip.text.y = element_text(size=12),
        strip.background = element_rect(colour=NA, fill=NA),
        legend.position = 'bottom')

g4
dev.off()

pdf("Compare_different_methods_df.pdf",width=6,height=3.5)
g4 <- res1 %>% filter(Method != "moderatedTMS" & Num_mixtures != 1) %>% 
  ggplot(aes(x=as.factor(Num_mixtures), y= DF, fill= Method))+
  #geom_jitter(position=position_jitter(0.1))+ # v1 : width=0.5, v2 n v3 : width=1
  geom_boxplot(lwd=0.6, outlier.alpha = 0.02)+ 
  scale_fill_manual(values=colors) + 
  labs(x=("# of mixtures"), y = ("Degree of freedoms"), size=12)+ 
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(angle=0, size=12), #, vjust = 0.43
        axis.text.x=element_text(angle=0, size=12),
        axis.title=element_text(size=12,face="bold"),
        strip.text.x = element_text(size=12), # angle=75
        strip.text.y = element_text(size=12),
        strip.background = element_rect(colour=NA, fill=NA),
        legend.position = 'bottom')

g4
dev.off()

pdf("Compare_different_methods_legend.pdf",width=8,height=3.5)
g4 <- res1 %>% filter(Method != "moderatedTMS" & Num_mixtures != 1) %>% 
  ggplot(aes(x=as.factor(Num_mixtures), y= DF, fill= Method))+
  #geom_jitter(position=position_jitter(0.1))+ # v1 : width=0.5, v2 n v3 : width=1
  geom_boxplot(lwd=0.6, outlier.alpha = 0.02)+ 
  scale_fill_manual(values=colors) + 
  labs(x=("# of mixtures"), y = ("Degree of freedoms"), size=12)+ 
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(angle=0, size=12), #, vjust = 0.43
        axis.text.x=element_text(angle=0, size=12),
        axis.title=element_text(size=12,face="bold"),
        strip.text.x = element_text(size=12), # angle=75
        strip.text.y = element_text(size=12),
        strip.background = element_rect(colour=NA, fill=NA),
        legend.position = 'bottom')

g4
dev.off()

