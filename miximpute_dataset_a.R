setwd("C:/Users/ys25.stu/OneDrive - UBC/BeeCSI/Proteomics/imputation")

#package
library(readr)
library(rlang)
library(tibble)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(matrixStats)
library(factoextra) #pca
library(imputeLCMD) #knn, svd, MLEimputation
library(pcaMethods) #bpca, lls imputation
library(missForest) #random forest
library(purrr)
library(svglite)
library(matrixStats)
library(xlsx)
library(pROC)
library(ROCR)
library(doParallel)
library(data.table)
library(cowplot)
library(limma)
library(ggbreak)
library(preprocessCore)
library(MKinfer)
library(pheatmap)
library(RColorBrewer)
library(scales)
library(sfsmisc)


#import peptide intensity data
hela <- read.delim("combined_modified_peptide_MBR.tsv")
hela <- hela[,-c(1, 3:7,9:10, 12:16)]

#extract Max LFQ values
hela_maxLFQ_intensity <- hela[, c(1:3, grep("MaxLFQ", colnames(hela)))]

#calculate peptide average intensity and missing rate
hela_maxLFQ_intensity <- add_column(hela_maxLFQ_intensity, mean_intensity = NA, .after = "Protein.ID")
hela_maxLFQ_intensity <- add_column(hela_maxLFQ_intensity, NA_count = NA, .after = "mean_intensity")
hela_maxLFQ_intensity <- add_column(hela_maxLFQ_intensity, NA_rate = NA, .after = "NA_count")
hela_maxLFQ_intensity$mean_intensity <- rowMeans(hela_maxLFQ_intensity[, -c(1:6)])
hela_maxLFQ_intensity$NA_count <- rowSums(hela_maxLFQ_intensity[, -c(1:6)] == 0)
hela_maxLFQ_intensity$NA_rate <- hela_maxLFQ_intensity$NA_count/(ncol(hela_maxLFQ_intensity) - 6)

#calculate 25% and 75% log2 intensity quantile
int_low <- log2(quantile(hela_maxLFQ_intensity$mean_intensity)[2])
int_high <- log2(quantile(hela_maxLFQ_intensity$mean_intensity)[4])

#log intensity vs NA rate correlation
ggplot(data = hela_maxLFQ_intensity[-which(hela_maxLFQ_intensity$NA_rate == 1),], aes(x = log2(mean_intensity), y = NA_rate)) + 
  geom_point(alpha = 0.1, shape = 20) +
  labs(title = paste("A. HeLa cell lysate QC dataset", "Peptide mean intensity vs. missing rate", sep = "\n"),
       x = bquote(paste("Mean log"["2"]*"(Intensity)")),
       y = "Peptide Missing Rate") +
  lims(x = c(3,25), y = c(0,1)) +
  theme_classic() +
  theme(axis.text = element_text(size = 15, color = "black"),
        axis.title = element_text(size = 15),
        title = element_text(size = 15))
last_plot() +
  stat_density2d(aes(fill = ..density..), geom = "tile", contour = FALSE, n=200, alpha = 0.8) +
  scale_fill_continuous(low = "white", high = "dodgerblue4") +
  labs(fill = "Density") +
  stat_smooth(method = "lm", geom = "smooth", formula = y~x, se = TRUE, color = "black", size = 0.5) +
  geom_line(aes(x = int_low), color = "red", size = 0.5) +
  geom_line(aes(x = int_high), color = "red", size = 0.5) +
  geom_line(aes(y = 0.75), color = "red", size = 0.5) +
  geom_line(aes(y = 0.25), color = "red", size = 0.5) +
  theme(legend.title = element_text(size = 13),
        legend.text = element_text(size = 12))
ggsave("figure1_int_mr_hela.svg",
       plot = last_plot(),
       path = "figures/manuscript/",
       device = "svg",
       width = 2244,
       height = 1496,
       units = "px",
       dpi = 300)
#check linear regression result
summary(lm(NA_rate~log2(mean_intensity), data = hela_maxLFQ_intensity[-which(hela_maxLFQ_intensity$NA_rate == 1),]))


#sepearte data by intensity and missing rate
#high intensity, high missing peptides
hihm <- hela_maxLFQ_intensity[which(hela_maxLFQ_intensity$mean_intensity >= int_high & 
                                      hela_maxLFQ_intensity$NA_rate >= 0.75 &
                                      hela_maxLFQ_intensity$NA_rate < 1),]
#high intensity, medium missing peptides
himm <- hela_maxLFQ_intensity[which(hela_maxLFQ_intensity$mean_intensity >= int_high & 
                                      hela_maxLFQ_intensity$NA_rate >= 0.25 &
                                      hela_maxLFQ_intensity$NA_rate < 0.75),]
#high intensity, low missing peptides
hilm <- hela_maxLFQ_intensity[which(hela_maxLFQ_intensity$mean_intensity >= int_high & 
                                      hela_maxLFQ_intensity$NA_rate < 0.25),]
#medium intensity, high missing peptides
mihm <- hela_maxLFQ_intensity[which(hela_maxLFQ_intensity$mean_intensity < int_high & 
                                      hela_maxLFQ_intensity$mean_intensity >= int_low &
                                      hela_maxLFQ_intensity$NA_rate >= 0.75 &
                                      hela_maxLFQ_intensity$NA_rate < 1),]
#medium intensity, medium missing peptides
mimm <- hela_maxLFQ_intensity[which(hela_maxLFQ_intensity$mean_intensity < int_high & 
                                      hela_maxLFQ_intensity$mean_intensity >= int_low &
                                      hela_maxLFQ_intensity$NA_rate >= 0.25 &
                                      hela_maxLFQ_intensity$NA_rate < 0.75),]
#medium intensity, low missing peptides
milm <- hela_maxLFQ_intensity[which(hela_maxLFQ_intensity$mean_intensity < int_high & 
                                      hela_maxLFQ_intensity$mean_intensity >= int_low &
                                      hela_maxLFQ_intensity$NA_rate < 0.25),]
#low intensity, high missing peptides
lihm <- hela_maxLFQ_intensity[which(hela_maxLFQ_intensity$mean_intensity < int_low & 
                                      hela_maxLFQ_intensity$NA_rate >= 0.75 &
                                      hela_maxLFQ_intensity$NA_rate < 1),]
#low intensity, medium missing peptides
limm <- hela_maxLFQ_intensity[which(hela_maxLFQ_intensity$mean_intensity < int_low & 
                                      hela_maxLFQ_intensity$NA_rate >= 0.25 &
                                      hela_maxLFQ_intensity$NA_rate < 0.75),]
#low intensity, low missing peptides
lilm <- hela_maxLFQ_intensity[which(hela_maxLFQ_intensity$mean_intensity < int_low & 
                                      hela_maxLFQ_intensity$NA_rate < 0.25),]
#lilm no data

#keep intensity
hihm <- hihm[,-c(2:6)]
himm <- himm[,-c(2:6)]
hilm <- hilm[,-c(2:6)]
mihm <- mihm[,-c(2:6)]
mimm <- mimm[,-c(2:6)]
milm <- milm[,-c(2:6)]
lihm <- lihm[,-c(2:6)]
limm <- limm[,-c(2:6)]

#replace 0 with NA
hihm[hihm == 0] <- NA
himm[himm == 0] <- NA
hilm[hilm == 0] <- NA
mihm[mihm == 0] <- NA
mimm[mimm == 0] <- NA
milm[milm == 0] <- NA
lihm[lihm == 0] <- NA
limm[limm == 0] <- NA

#peptide name as rownames
rownames(hihm) <- hihm[,1]
rownames(himm) <- himm[,1]
rownames(hilm) <- hilm[,1]
rownames(mihm) <- mihm[,1]
rownames(mimm) <- mimm[,1]
rownames(milm) <- milm[,1]
rownames(lihm) <- lihm[,1]
rownames(limm) <- limm[,1]
#remove first column
hihm <- hihm[,-1]
himm <- himm[,-1]
hilm <- hilm[,-1]
mihm <- mihm[,-1]
mimm <- mimm[,-1]
milm <- milm[,-1]
lihm <- lihm[,-1]
limm <- limm[,-1]


#prepare datasets for imputation, remove sample and pepetide that are all NA
hihm_na <- hihm[which(rowSums(is.na(hihm)) < ncol(hihm)),which(colSums(is.na(hihm)) < nrow(hihm))]
himm_na <- himm[which(rowSums(is.na(himm)) < ncol(himm)),which(colSums(is.na(himm)) < nrow(himm))]
hilm_na <- hilm[which(rowSums(is.na(hilm)) < ncol(hilm)),which(colSums(is.na(hilm)) < nrow(hilm))]
mihm_na <- mihm[which(rowSums(is.na(mihm)) < ncol(mihm)),which(colSums(is.na(mihm)) < nrow(mihm))]
mimm_na <- mimm[which(rowSums(is.na(mimm)) < ncol(mimm)),which(colSums(is.na(mimm)) < nrow(mimm))]
milm_na <- milm[which(rowSums(is.na(milm)) < ncol(milm)),which(colSums(is.na(milm)) < nrow(milm))]
lihm_na <- lihm[which(rowSums(is.na(lihm)) < ncol(lihm)),which(colSums(is.na(lihm)) < nrow(lihm))]
limm_na <- limm[which(rowSums(is.na(limm)) < ncol(limm)),which(colSums(is.na(limm)) < nrow(limm))]

colnames(limm)[which(!colnames(limm) %in% colnames(limm_na))]
#11837, 11333 is missing 100% peptide in limm and milm, so remove them from all
hihm_na <- hihm_na[,-grep(paste(c("11837", "11333"), collapse = "|"), colnames(hihm_na))]
himm_na <- himm_na[,-grep(paste(c("11837", "11333"), collapse = "|"), colnames(himm_na))]
hilm_na <- hilm_na[,-grep(paste(c("11837", "11333"), collapse = "|"), colnames(hilm_na))]
mihm_na <- mihm_na[,-grep(paste(c("11837", "11333"), collapse = "|"), colnames(mihm_na))]
mimm_na <- mimm_na[,-grep(paste(c("11837", "11333"), collapse = "|"), colnames(mimm_na))]
milm_na <- milm_na[,-grep(paste(c("11837", "11333"), collapse = "|"), colnames(milm_na))]
lihm_na <- lihm_na[,-grep(paste(c("11837", "11333"), collapse = "|"), colnames(lihm_na))]

#log2 transformation
hihm_na <- hihm_na %>% log2()
himm_na <- himm_na %>% log2()
hilm_na <- hilm_na %>% log2()
mihm_na <- mihm_na %>% log2()
mimm_na <- mimm_na %>% log2()
milm_na <- milm_na %>% log2()
lihm_na <- lihm_na %>% log2()
limm_na <- limm_na %>% log2()

#count non NA
sum(!is.na(hihm_real))
sum(!is.na(himm_real))
sum(!is.na(hilm_real))
sum(!is.na(mihm_real))
sum(!is.na(mimm_real))
sum(!is.na(milm_real))
sum(!is.na(lihm_real))
sum(!is.na(limm_real))


#randomly divide data, 60% as training set, 40% as testing set
set.seed(100)
trs <- sample(1:80, 48, replace = FALSE)

#parallel computing for random forest
registerDoParallel(cores = 18)

for (i in 1:length(regions)) {
  temp <- get(paste(regions[i], "na", sep = "_"))
  #randomly divide data, 60% as training set
  temp1 <- temp[, trs]
  #40% as testing set
  temp2 <- temp[, -trs]
  
  #keep real dataset
  assign(paste(regions[i], "trs", "real", sep = "_"), temp1) #training set
  assign(paste(regions[i], "tes", "real", sep = "_"), temp2) #testing set
  print(paste(regions[i], "spliting complete", sep = " "))
  
  #masking 10% of the non-NA value
  #training set
  na_pos <- sample(length(!temp1), 0.1*length(!temp1))
  temp1 <- as.matrix(temp1)
  temp1[na_pos] <- NA
  temp1 <- as.data.frame(temp1)
  #testing set
  na_pos <- sample(length(!temp2), 0.1*length(!temp2))
  temp2 <- as.matrix(temp2)
  temp2[na_pos] <- NA
  temp2 <- as.data.frame(temp2)
  
  #count na percentage, remove peptides that are all missing after masking
  if (length(which(rowSums(is.na(temp1)) == ncol(temp1))) > 0) {
    temp1 <- temp1[-which(rowSums(is.na(temp1)) == ncol(temp1)),]
  }
  if (length(which(rowSums(is.na(temp2)) == ncol(temp2))) > 0) {
    temp2 <- temp2[-which(rowSums(is.na(temp2)) == ncol(temp2)),]
  }
  
  #keep the na datasets
  assign(paste(regions[i], "trs", "na", sep = "_"), temp1)
  assign(paste(regions[i], "tes", "na", sep = "_"), temp2)
  print(paste(regions[i], "masking complete", sep = " "))
  
  #bpca
  #training set
  temp11 <- pca(temp1, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
  temp11 <- as.data.frame(temp11@completeObs)
  assign(paste(regions[i], "trs", "bpca", sep = "_"), temp11)
  #testing set
  temp21 <- pca(temp2, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
  temp21 <- as.data.frame(temp21@completeObs)
  assign(paste(regions[i], "tes", "bpca", sep = "_"), temp21)
  print(paste(regions[i], "bpca complete", sep = " "))
  
  #svd
  #training set
  temp11 <- as.data.frame(impute.wrapper.SVD(temp1, K = 5))
  assign(paste(regions[i], "trs", "svd", sep = "_"), temp11)
  #testing set
  temp21 <- as.data.frame(impute.wrapper.SVD(temp2, K = 5))
  assign(paste(regions[i], "tes", "svd", sep = "_"), temp21)
  print(paste(regions[i], "svd complete", sep = " "))
  
  #knn
  #trainging set
  temp11 <- as.data.frame(impute.wrapper.KNN(as.matrix(temp1), K = 5))
  assign(paste(regions[i], "trs", "knn", sep = "_"), temp11)
  #testing set
  temp21 <- as.data.frame(impute.wrapper.KNN(as.matrix(temp2), K = 5))
  assign(paste(regions[i], "tes", "knn", sep = "_"), temp21)
  print(paste(regions[i], "knn complete", sep = " "))
  
  #mle
  temp11 <- as.data.frame(impute.wrapper.MLE(as.matrix(temp1)))
  assign(paste(regions[i], "trs", "mle", sep = "_"), temp11)
  temp21 <- as.data.frame(impute.wrapper.MLE(as.matrix(temp2)))
  assign(paste(regions[i], "tes", "mle", sep = "_"), temp21)
  print(paste(regions[i], "mle complete", sep = " "))
  
  #lls
  temp11 <- llsImpute(temp1, k = 5, center = TRUE, completeObs = TRUE)
  temp11 <- as.data.frame(temp11@completeObs)
  assign(paste(regions[i], "trs", "lls", sep = "_"), temp11)
  temp21 <- llsImpute(temp2, k = 5, center = TRUE, completeObs = TRUE)
  temp21 <- as.data.frame(temp21@completeObs)
  assign(paste(regions[i], "tes", "lls", sep = "_"), temp21)
  print(paste(regions[i], "lls complete", sep = " "))
  
  #rf
  temp11 <- missForest(temp1, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
  temp11 <- temp11$ximp
  assign(paste(regions[i], "trs", "rf", sep = "_"), temp11)
  temp21 <- missForest(temp2, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
  temp21 <- temp21$ximp
  assign(paste(regions[i], "tes", "rf", sep = "_"), temp21)
  print(paste(regions[i], "rf complete", sep = " "))
  
  #transpose data frame for pimms
  write.csv(t(temp1), file = paste(regions[i], "trs_na.csv", sep = "_"))
  write.csv(t(temp2), file = paste(regions[i], "tes_na.csv", sep = "_"))
}

#import imputed result from pimms
#cf, vae, dae
datatype <- c("trs", "tes")
regions <- c("hihm", "himm", "hilm", "mihm", "mimm", "milm", "lihm", "limm", "whole")

for (i in 1:length(regions)) {
  for (j in 1:length(datatype)) {
    #import cf
    if (length(grep(paste(regions[i], datatype[j], "cf", sep = "_"), dir())) > 0) {
      temp1 <- read.csv(paste(regions[i], datatype[j], "cf.csv", sep = "_"))
      rownames(temp1) <- temp1$protein.group
      temp1 <- temp1[,-1]
      assign(paste(regions[i], datatype[j], "cf", sep = "_"), temp1)
    }
    #import dae
    if (length(grep(paste(regions[i], datatype[j], "dae", sep = "_"), dir())) > 0) {
      temp2 <- read.csv(paste(regions[i], datatype[j], "dae.csv", sep = "_"))
      rownames(temp2) <- temp2$protein.group
      temp2 <- temp2[,-c(1:2)]
      assign(paste(regions[i], datatype[j], "dae", sep = "_"), temp2)
    }
    #import vae
    if (length(grep(paste(regions[i], datatype[j], "vae", sep = "_"), dir())) > 0) {
      temp3 <- read.csv(paste(regions[i], datatype[j], "vae.csv", sep = "_"))
      rownames(temp3) <- temp3$protein.group
      temp3 <- temp3[,-c(1:2)]
      assign(paste(regions[i], datatype[j], "vae", sep = "_"), temp3)
    }
    rm(temp1)
    rm(temp2)
    rm(temp3)
  }
}

#NRMSE evaluation by region
for (i in 1:length(regions)) {
  for (j in 1:length(datatype)) {
    #get na data
    comp_na <- get(paste(regions[i], datatype[j], "na", sep = "_"))
    #get real data
    comp_real <- get(paste(regions[i], datatype[j], "real", sep = "_"))
    comp_real <- comp_real[which(rownames(comp_real) %in% rownames(comp_na)),]
    
    #NRMSE
    nrmse_sum <- data.frame(sample = colnames(comp_na),
                            nrmse_bpca = NA,
                            nrmse_knn = NA,
                            nrmse_lls = NA,
                            nrmse_mle = NA,
                            nrmse_svd = NA, 
                            nrmse_rf = NA,
                            nrmse_cf = NA,
                            nrmse_dae = NA,
                            nrmse_vae = NA)
    
    #get data for different methods
    for (k in 1:length(methods)) {
      if (paste(regions[i], datatype[j], methods[k], sep = "_") %in% ls() == TRUE) {
        temp <- get(paste(regions[i], datatype[j], methods[k], sep = "_"))
        for (l in 1:length(nrmse_sum$sample)) {
          nrmse_sum[l, grep(methods[k], colnames(nrmse_sum))] <- sqrt(mean((temp[,l] - comp_real[,l])^2, na.rm = TRUE))/sd(comp_real[,l], na.rm = TRUE)
        }
      } else {
        next
      }
      #if highest value > 4 fold of median, remove values higher than 90% values in a method
      if (max(nrmse_sum[,grep(methods[k], colnames(nrmse_sum))]) > 4*quantile(nrmse_sum[,grep(methods[k], colnames(nrmse_sum))], probs = c(0,0.5,0.9,0.95,1))[2]) {
        print(paste(regions[i], datatype[j], methods[k], sep = "_"))
        print(quantile(nrmse_sum[,grep(methods[k], colnames(nrmse_sum))], probs = c(0,0.5,0.9,0.95,1)))
        nrmse_sum[which(nrmse_sum[,grep(methods[k], colnames(nrmse_sum))] > 5*quantile(nrmse_sum[,grep(methods[k], colnames(nrmse_sum))], probs = c(0,0.5,0.9,0.95,1))[2]),
                  grep(methods[k], colnames(nrmse_sum))] <- NA
      }
    }
    nrmse_long <- nrmse_sum %>% pivot_longer(cols = 2:length(.), values_to = "NRMSE", names_to = "method")
    nrmse_long$method <- factor(nrmse_long$method, levels = c("nrmse_bpca", "nrmse_knn", "nrmse_lls", "nrmse_mle", 
                                                              "nrmse_svd", "nrmse_rf", "nrmse_cf", "nrmse_dae", "nrmse_vae"))
    #barplot
    #calculate mean and standard error
    nrmse_stat <- data.frame(method = colnames(nrmse_sum)[-c(1)],
                             mean = colMeans(nrmse_sum[,-c(1)], na.rm = TRUE),
                             se = colSds(as.matrix(nrmse_sum[,-c(1)]), na.rm = TRUE)/sqrt(nrow(nrmse_sum)))
    nrmse_stat$method <- factor(nrmse_stat$method, levels = c("nrmse_rf", "nrmse_bpca", "nrmse_knn", "nrmse_lls",  "nrmse_mle",
                                                              "nrmse_svd", "nrmse_cf", "nrmse_dae", "nrmse_vae"))
    nrmse_anova <- anova(lm(NRMSE ~ method, data = nrmse_long))
    TukeyHSD(aov(lm(NRMSE ~ method, data = nrmse_long)))
    
    #barplot of NRMSE
    if (max(nrmse_stat$mean, na.rm = TRUE) > 10) {
      ggplot(data = nrmse_stat, aes(x = method, y = mean, fill = method)) +
        geom_bar(stat = "identity") +
        geom_errorbar(aes(ymin = mean - 1.96*se, ymax = mean + 1.96*se, width = 0.2)) + #95% confidence interval
        geom_text(aes(label = round(mean, digits = 2)), vjust = -0.5, size = 5) +
        scale_fill_brewer(palette = "Set3", labels = c("RF", "BPCA", "kNN", "LLS", "MLE", "SVD", "CF", "DAE", "VAE")) +
        scale_y_break(c(1.5, 200), scales = 0.5) +
        labs(x = toupper(regions[i]),
             y = "NRMSE") +
        theme_bw() +
        theme(axis.title.x = element_text(face = "bold", size = 20),
              axis.title.y = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_text(face = "bold", size = 20),
              panel.grid.major.x = element_blank(),
              panel.grid.major = element_line(linewidth = 1),
              legend.position = "none")
    } else {
      ggplot(data = nrmse_stat, aes(x = method, y = mean, fill = method)) +
        geom_bar(stat = "identity") +
        geom_errorbar(aes(ymin = mean - 1.96*se, ymax = mean + 1.96*se, width = 0.2)) + #95% confidence interval
        geom_text(aes(label = round(mean, digits = 2)), vjust = -0.8, size = 5) +
        scale_fill_brewer(palette = "Set3", labels = c("RF", "BPCA", "kNN", "LLS", "MLE", "SVD", "CF", "DAE", "VAE")) +
        ylim(0,1.5) +
        labs(x = toupper(regions[i]),
             y = "NRMSE") +
        theme_bw() +
        theme(axis.title.x = element_text(face = "bold", size = 20),
              axis.title.y = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_text(face = "bold", size = 20),
              panel.grid.major.x = element_blank(),
              panel.grid.major = element_line(linewidth = 1),
              legend.position = "none")
    }
    ggsave(paste("figures/manuscript/Hela/NRMSE_barplot_", paste(regions[i], datatype[j], sep = "_"), ".jpeg", sep = ""), device = "jpeg", width = 1400, height = 1400, units = "px", dpi = 300)
    print(regions[i])
    print(colMeans(nrmse_sum[,2:length(nrmse_sum)], na.rm = TRUE))
    quants <- c(0, 0.25, 0.50, 0.75, 1)
    print(apply(nrmse_sum[,2:length(nrmse_sum)], 2 , quantile , probs = quants , na.rm = TRUE))
    print(paste(gsub("nrmse_", "", nrmse_stat$method[which(nrmse_stat$mean==min(nrmse_stat$mean, na.rm = TRUE))]),
                "is the imputation method with lowest NRMSE in",
                regions[i],
                datatype[j],
                sep = " "))
  }
}

#ANOVA of NRMSE
aov_sum <- data.frame(region = toupper(rep(regions, length(datatype))),
                      datatype = c(rep("trs", length(regions)), rep("tes", length(regions))),
                      Df = NA,
                      sum_sq = NA,
                      mean_sq = NA,
                      F_value = NA,
                      p_value = NA)
for (i in 1:length(regions)) {
  for (j in 1:length(datatype)) {
    #get na data
    comp_na <- get(paste(regions[i], datatype[j], "na", sep = "_"))
    #get real data
    comp_real <- get(paste(regions[i], datatype[j], "real", sep = "_"))
    comp_real <- comp_real[which(rownames(comp_real) %in% rownames(comp_na)),]
    
    #NRMSE
    nrmse_sum <- data.frame(sample = colnames(comp_na),
                            nrmse_bpca = NA,
                            nrmse_knn = NA,
                            nrmse_lls = NA,
                            nrmse_mle = NA,
                            nrmse_svd = NA, 
                            nrmse_rf = NA,
                            nrmse_cf = NA,
                            nrmse_dae = NA,
                            nrmse_vae = NA)
    
    #get data for different methods
    for (k in 1:length(methods)) {
      if (paste(regions[i], datatype[j], methods[k], sep = "_") %in% ls() == TRUE) {
        temp <- get(paste(regions[i], datatype[j], methods[k], sep = "_"))
        for (l in 1:length(nrmse_sum$sample)) {
          nrmse_sum[l, grep(methods[k], colnames(nrmse_sum))] <- sqrt(mean((temp[,l] - comp_real[,l])^2, na.rm = TRUE))/sd(comp_real[,l], na.rm = TRUE)
        }
      } else {
        next
      }
      #if highest value > 4 fold of median, remove values higher than 90% values in a method
      if (max(nrmse_sum[,grep(methods[k], colnames(nrmse_sum))]) > 4*quantile(nrmse_sum[,grep(methods[k], colnames(nrmse_sum))], probs = c(0,0.5,0.9,0.95,1))[2]) {
        print(paste(regions[i], datatype[j], methods[k], sep = "_"))
        print(quantile(nrmse_sum[,grep(methods[k], colnames(nrmse_sum))], probs = c(0,0.5,0.9,0.95,1)))
        nrmse_sum[which(nrmse_sum[,grep(methods[k], colnames(nrmse_sum))] > 5*quantile(nrmse_sum[,grep(methods[k], colnames(nrmse_sum))], probs = c(0,0.5,0.9,0.95,1))[2]),
                  grep(methods[k], colnames(nrmse_sum))] <- NA
      }
    }
    nrmse_long <- nrmse_sum %>% pivot_longer(cols = 2:length(.), values_to = "NRMSE", names_to = "method")
    nrmse_long$method <- factor(nrmse_long$method, levels = c("nrmse_bpca", "nrmse_knn", "nrmse_lls", "nrmse_mle", 
                                                              "nrmse_svd", "nrmse_rf", "nrmse_cf", "nrmse_dae", "nrmse_vae"))
    #barplot
    #calculate mean and standard error
    nrmse_stat <- data.frame(method = colnames(nrmse_sum)[-c(1)],
                             mean = colMeans(nrmse_sum[,-c(1)], na.rm = TRUE),
                             se = colSds(as.matrix(nrmse_sum[,-c(1)]), na.rm = TRUE)/sqrt(nrow(nrmse_sum)))
    nrmse_stat$method <- factor(nrmse_stat$method, levels = c("nrmse_rf", "nrmse_bpca", "nrmse_knn", "nrmse_lls",  "nrmse_mle",
                                                              "nrmse_svd", "nrmse_cf", "nrmse_dae", "nrmse_vae"))
    
    nrmse_anova <- anova(lm(NRMSE ~ method, data = nrmse_long))
    aov_sum$Df[intersect(which(aov_sum$region == regions[i]), which(aov_sum$datatype == datatype[j]))] <- nrmse_anova$Df[1]
    aov_sum$sum_sq[intersect(which(aov_sum$region == regions[i]), which(aov_sum$datatype == datatype[j]))] <- round(nrmse_anova$`Sum Sq`[1], digits = 3)
    aov_sum$mean_sq[intersect(which(aov_sum$region == regions[i]), which(aov_sum$datatype == datatype[j]))] <- round(nrmse_anova$`Mean Sq`[1], digits = 3)  
    aov_sum$F_value[intersect(which(aov_sum$region == regions[i]), which(aov_sum$datatype == datatype[j]))] <- round(nrmse_anova$`F value`[1], digits = 3)  
    aov_sum$p_value[intersect(which(aov_sum$region == regions[i]), which(aov_sum$datatype == datatype[j]))] <- nrmse_anova$`Pr(>F)`[1]
  }
}

write.csv(aov_sum, "aov_sum_hela_region.csv")


#assemble best practice from each region
#training set
whole_trs_mix <- rbind(hihm_trs_rf, himm_trs_rf, hilm_trs_rf, mihm_trs_bpca, mimm_trs_rf, milm_trs_rf, lihm_trs_bpca, limm_trs_rf)
#testing set
whole_tes_mix <- rbind(hihm_tes_rf, himm_tes_rf, hilm_tes_rf, mihm_tes_bpca, mimm_tes_rf, milm_tes_rf, lihm_tes_bpca, limm_tes_rf)

#mpute as whole
#training set
whole_trs_na <- rbind(hihm_trs_na, himm_trs_na, hilm_trs_na,
                      mihm_trs_na, mimm_trs_na, milm_trs_na,
                      lihm_trs_na, limm_trs_na)
#testing set
whole_tes_na <- rbind(hihm_tes_na, himm_tes_na, hilm_tes_na,
                      mihm_tes_na, mimm_tes_na, milm_tes_na,
                      lihm_tes_na, limm_tes_na)

colnames(hihm_trs_na)%in% colnames(himm_trs_na)

whole_tes_na <- rbind(hihm_tes_na, himm_tes_na, hilm_tes_na,
                      mihm_tes_na, mimm_tes_na, milm_tes_na,
                      lihm_tes_na, limm_tes_na)
#export for pimms
write.csv(t(whole_trs_na), file = "whole_trs_na.csv")
write.csv(t(whole_tes_na), file = "whole_tes_na.csv")

#real
#training set
whole_trs_real <- rbind(hihm_trs_real[rownames(hihm_trs_real) %in% rownames(hihm_trs_na),],
                        himm_trs_real[rownames(himm_trs_real) %in% rownames(himm_trs_na),],
                        hilm_trs_real[rownames(hilm_trs_real) %in% rownames(hilm_trs_na),],
                        mihm_trs_real[rownames(mihm_trs_real) %in% rownames(mihm_trs_na),],
                        mimm_trs_real[rownames(mimm_trs_real) %in% rownames(mimm_trs_na),],
                        milm_trs_real[rownames(milm_trs_real) %in% rownames(milm_trs_na),],
                        lihm_trs_real[rownames(lihm_trs_real) %in% rownames(lihm_trs_na),],
                        limm_trs_real[rownames(limm_trs_real) %in% rownames(limm_trs_na),])

whole_tes_real <- rbind(hihm_tes_real[rownames(hihm_tes_real) %in% rownames(hihm_tes_na),],
                        himm_tes_real[rownames(himm_tes_real) %in% rownames(himm_tes_na),],
                        hilm_tes_real[rownames(hilm_tes_real) %in% rownames(hilm_tes_na),],
                        mihm_tes_real[rownames(mihm_tes_real) %in% rownames(mihm_tes_na),],
                        mimm_tes_real[rownames(mimm_tes_real) %in% rownames(mimm_tes_na),],
                        milm_tes_real[rownames(milm_tes_real) %in% rownames(milm_tes_na),],
                        lihm_tes_real[rownames(lihm_tes_real) %in% rownames(lihm_tes_na),],
                        limm_tes_real[rownames(limm_tes_real) %in% rownames(limm_tes_na),])

#training set
#bpca
whole_trs_bpca <- pca(whole_trs_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
whole_trs_bpca <- as.data.frame(whole_trs_bpca@completeObs)
#svd
whole_trs_svd <- as.data.frame(impute.wrapper.SVD(whole_trs_na, K = 5))
#knn
whole_trs_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(whole_trs_na), K = 5))
#mle
whole_trs_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(whole_trs_na)))
#lls
whole_trs_lls <- llsImpute(whole_trs_na, k = 5, center = TRUE, completeObs = TRUE)
whole_trs_lls <- as.data.frame(whole_trs_lls@completeObs)
#rf
whole_trs_rf <- missForest(whole_trs_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
whole_trs_rf <- whole_trs_rf$ximp

#testing set
#bpca
whole_tes_bpca <- pca(whole_tes_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
whole_tes_bpca <- as.data.frame(whole_tes_bpca@completeObs)
#svd
whole_tes_svd <- as.data.frame(impute.wrapper.SVD(whole_tes_na, K = 5))
#knn
whole_tes_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(whole_tes_na), K = 5))
#mle
whole_tes_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(whole_tes_na)))
#lls
whole_tes_lls <- llsImpute(whole_tes_na, k = 5, center = TRUE, completeObs = TRUE)
whole_tes_lls <- as.data.frame(whole_tes_lls@completeObs)
#rf
whole_tes_rf <- missForest(whole_tes_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
whole_tes_rf <- whole_tes_rf$ximp

#pimms result imported above

#NRMSE
regions <- c("whole")
datatype <- c("trs", "tes")
methods <- c("bpca", "knn", "lls", "mle", "svd", "rf", "cf", "dae", "vae", "mix")

#evaluation by region
for (i in 1:length(regions)) {
  for (j in 1:length(datatype)) {
    #get na data
    comp_na <- get(paste(regions[i], datatype[j], "na", sep = "_"))
    #get real data
    comp_real <- get(paste(regions[i], datatype[j], "real", sep = "_"))
    comp_real <- comp_real[which(rownames(comp_real) %in% rownames(comp_na)),]
    
    #NRMSE
    nrmse_sum <- data.frame(sample = colnames(comp_na),
                            nrmse_bpca = NA,
                            nrmse_knn = NA,
                            nrmse_lls = NA,
                            nrmse_mle = NA,
                            nrmse_svd = NA, 
                            nrmse_rf = NA,
                            nrmse_cf = NA,
                            nrmse_dae = NA,
                            nrmse_vae = NA,
                            nrmse_mix = NA)
    
    #get data for different methods
    for (k in 1:length(methods)) {
      if (paste(regions[i], datatype[j], methods[k], sep = "_") %in% ls() == TRUE) {
        temp <- get(paste(regions[i], datatype[j], methods[k], sep = "_"))
        for (l in 1:length(nrmse_sum$sample)) {
          nrmse_sum[l, grep(methods[k], colnames(nrmse_sum))] <- sqrt(mean((temp[,l] - comp_real[,l])^2, na.rm = TRUE))/sd(comp_real[,l], na.rm = TRUE)
        }
      } else {
        next
      }
      #if highest value > 4 fold of median, remove values higher than 90% values in a method
      if (max(nrmse_sum[,grep(methods[k], colnames(nrmse_sum))]) > 4*quantile(nrmse_sum[,grep(methods[k], colnames(nrmse_sum))], probs = c(0,0.5,0.9,0.95,1))[2]) {
        print(paste(regions[i], datatype[j], methods[k], sep = "_"))
        print(quantile(nrmse_sum[,grep(methods[k], colnames(nrmse_sum))], probs = c(0,0.5,0.9,0.95,1)))
        nrmse_sum[which(nrmse_sum[,grep(methods[k], colnames(nrmse_sum))] > 5*quantile(nrmse_sum[,grep(methods[k], colnames(nrmse_sum))], probs = c(0,0.5,0.9,0.95,1))[2]),
                  grep(methods[k], colnames(nrmse_sum))] <- NA
      }
    }
    nrmse_long <- nrmse_sum %>% pivot_longer(cols = 2:length(.), values_to = "NRMSE", names_to = "method")
    nrmse_long$method <- factor(nrmse_long$method, levels = c("nrmse_bpca", "nrmse_knn", "nrmse_lls", "nrmse_mle", 
                                                              "nrmse_svd", "nrmse_rf", "nrmse_cf", "nrmse_dae", "nrmse_vae", "nrmse_mix"))
    nrmse_anova <- anova(lm(NRMSE ~ method, data = nrmse_long))
    nrmse_tukey <- TukeyHSD(aov(lm(NRMSE ~ method, data = nrmse_long)))$method
    write.csv(nrmse_anova, paste("nrmse_anova_", datatype[j], ".csv", sep = ""))
    write.csv(nrmse_tukey, paste("nrmse_tukey_", datatype[j], ".csv", sep = ""))
    
    #barplot
    #calculate mean and standard error
    nrmse_stat <- data.frame(method = colnames(nrmse_sum)[-1],
                             mean = colMeans(nrmse_sum[,-1], na.rm = TRUE),
                             se = colSds(as.matrix(nrmse_sum[,-1]), na.rm = TRUE)/sqrt(nrow(nrmse_sum)))
    nrmse_stat$method <- factor(nrmse_stat$method, levels = c("nrmse_rf", "nrmse_bpca", "nrmse_knn", "nrmse_lls", "nrmse_mle",
                                                              "nrmse_svd", "nrmse_cf", "nrmse_dae", "nrmse_vae", "nrmse_mix"))
    
    ggplot(data = nrmse_stat, aes(x = method, y = mean, fill = method)) +
      geom_bar(stat = "identity") +
      geom_errorbar(aes(ymin = mean - 1.96*se, ymax = mean + 1.96*se, width = 0.2)) + #95% confidence interval
      geom_text(aes(label = round(mean, digits = 2)), vjust = -1, size = 5) +
      scale_x_discrete(labels = c("RF", "BPCA", "kNN", "LLS", "MLE", "SVD", "CF", "DAE", "VAE", "Mix")) +
      labs(x = "Method",
           y = "NRMSE") +
      ylim(0,0.7) +
      scale_fill_brewer(labels = c("RF", "BPCA", "kNN", "LLS", "MLE", "SVD", "CF", "DAE", "VAE", "Mix"), palette = "Set3") +
      theme_bw() +
      theme(axis.title = element_text(size = 15),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(), 
            panel.grid.major.x = element_blank(),
            panel.grid.major = element_line(linewidth = 1),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 15))
    
    ggsave(paste("figures/manuscript/Hela/NRMSE_barplot_", paste(regions[i], datatype[j], sep = "_"), ".jpeg", sep = ""), device = "jpeg", width = 1800, height = 1000, units = "px", dpi = 300)
    print(regions[i])
    print(colMeans(nrmse_sum[,2:length(nrmse_sum)], na.rm = TRUE))
    quants <- c(0, 0.25, 0.50, 0.75, 1)
    print(apply(nrmse_sum[,2:length(nrmse_sum)], 2 , quantile , probs = quants , na.rm = TRUE))
    print(paste(gsub("nrmse_", "", nrmse_stat$method[which(nrmse_stat$mean==min(nrmse_stat$mean, na.rm = TRUE))]),
                "is the imputation method with lowest NRMSE in",
                regions[i],
                datatype[j],
                sep = " "))
  }
}

#ANOVA
comp_real <- whole_trs_real[which(rownames(whole_trs_real) %in% rownames(whole_trs_na)),] %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "real")
comp_na <- whole_trs_na %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "na")
comp_knn <- whole_trs_knn %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "knn")
comp_mle <- whole_trs_mle %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "mle")
comp_svd <- whole_trs_svd %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "svd")
comp_bpca <- whole_trs_bpca %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "bpca")
comp_lls <- whole_trs_lls %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "lls")
comp_rf <- whole_trs_rf %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "rf")
comp_cf <- whole_trs_cf %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "cf")
comp_vae <- whole_trs_vae %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "vae")
comp_dae <- whole_trs_dae %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "dae")
comp_mix <- whole_trs_mix %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "mix")

comp <- data.frame(sample = comp_real$sample,
                   real = comp_real$real,
                   na = comp_na$na,
                   knn = comp_knn$knn,
                   mle = comp_mle$mle,
                   svd = comp_svd$svd,
                   bpca = comp_bpca$bpca,
                   lls = comp_lls$lls,
                   rf = comp_rf$rf,
                   cf = comp_cf$cf,
                   vae = comp_vae$vae,
                   dae = comp_dae$dae,
                   mix = comp_mix$mix)

#anova
comp_long <- comp[!is.na(comp$na),-3] %>% pivot_longer(., c(2, 3:12), names_to = "method", values_to = "value")
comp_long$method <- as.factor(comp_long$method)

#anova
result_anova <- aov(value~method, data = comp_long)
summary(result_anova)
#pairwise comparison
result_tukey <- TukeyHSD(result_anova)


result_tukey$method[grep("real", rownames(result_tukey$method)),]


res_cor <- cor(as.matrix(comp[!is.na(comp$real),-c(1,3)]), method = "spearman")
summary(res_cor)

cor.test(as.matrix(comp[!is.na(comp$real),-c(1,3)]), method = "spearman")

res_cor_p <- data.frame(method = colnames(comp)[4:13],
                        r = NA,
                        p = NA)


for (i in 1:nrow(res_cor_p)) {
  res_cor_p$r[i] <- cor.test(comp$real[!is.na(comp$real)], comp[!is.na(comp$real),i+3], method = "pearson")
  res_cor_p$r[i] <- cor.test(!is.na(comp$real), comp[!is.na(comp$real),i+3], method = "spearman")
  
}

