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
library(plotROC)
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

#validation using dataset b
PXD000279 <- read_delim("~/datasets/PXD000279_maxlfq_benchmark_human_ecoli_mixture/peptides.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
PXD000279_peptide <- as.data.frame(PXD000279[,c(1,33,74:79)])
rownames(PXD000279_peptide) <- PXD000279_peptide[,1]
PXD000279_peptide <- PXD000279_peptide[,-1]
PXD000279_peptide[,-1] <- log2(PXD000279_peptide[,-1])
PXD000279_peptide[PXD000279_peptide == -Inf] <- NA

PXD000279_proteinGroups <- read_delim("reference/Pietz.2024/data/datasets/PXD000279_maxlfq_benchmark_human_ecoli_mixture/proteinGroups.txt", 
                                      +     delim = "\t", escape_double = FALSE, 
                                      +     trim_ws = TRUE)

#calculate NA rate, mean intensity
PXD000279_peptide <- add_column(PXD000279_peptide, mean_intensity = NA, .before = "Intensity H1")
PXD000279_peptide <- add_column(PXD000279_peptide, NA_rate = NA, .after = "mean_intensity")

PXD000279_peptide$mean_intensity <- rowMeans(PXD000279_peptide[, -c(1,2,3)], na.rm = TRUE)
PXD000279_peptide$NA_rate <- rowSums(is.na(PXD000279_peptide[,-c(1,2,3)]))/(ncol(PXD000279_peptide) - 3)

PXD000279_int_low <- quantile(PXD000279_peptide$mean_intensity, na.rm = TRUE)[2]
PXD000279_int_high <- quantile(PXD000279_peptide$mean_intensity, na.rm = TRUE)[4]

#NA rate vs intensity
ggplot(data = PXD000279_peptide, aes(x = mean_intensity, y = NA_rate)) + 
  geom_point(alpha = 0.1, shape = 20) +
  labs(title = paste("B. PXD000279 Human E. coli mixture", "Peptide mean intensity and missing rate", sep = "\n"),
       x = bquote(paste("mean log"["2"]*"(intensity)")),
       y = "Missing rate") +
  lims(y = c(0,1)) +
  theme_classic() +
  theme(axis.text = element_text(size = 15, color = "black"),
        axis.title = element_text(size = 15),
        title = element_text(size = 15)) 
last_plot() +
  stat_density2d(aes(fill = ..density..), geom = "tile", contour = FALSE, n=200, alpha = 0.8) +
  scale_fill_continuous(low = "white", high = "dodgerblue4") +
  stat_smooth(method = "lm", geom = "smooth", formula = y~x, se = TRUE, color = "black", size = 0.5) +
  labs(fill = "Density") +
  geom_line(aes(x = PXD000279_int_high), color = "red", size = 0.5) +
  geom_line(aes(x = PXD000279_int_low), color = "red", size = 0.5) +
  geom_line(aes(y = 0.65), color = "red", size = 0.5) +
  geom_line(aes(y = 0.35), color = "red", size = 0.5) +
  theme(legend.title = element_text(size = 13),
        legend.text = element_text(size = 12))
ggsave("figure3_int_mr_PXD000279.tiff",
       plot = last_plot(),
       path = "figures/manuscript/",
       device = "tiff",
       width = 2244,
       height = 1496,
       units = "px",
       dpi = 300)
ggsave("figure3_int_mr_PXD000279.jpeg",
       plot = last_plot(),
       path = "figures/manuscript/",
       device = "jpeg",
       width = 2244,
       height = 1496,
       units = "px",
       dpi = 300)

summary(lm(NA_rate~mean_intensity, data = PXD000279_peptide))

#keep mapping relationship between peptide and protein
ppmap <- data.frame(peptide = rownames(PXD000279_peptide),
                    protein = PXD000279_peptide$`Leading razor protein`)

#high intensity, high missing peptides
PXD000279_hihm <- PXD000279_peptide[which(PXD000279_peptide$mean_intensity >= PXD000279_int_high & 
                                            PXD000279_peptide$NA_rate >= 0.65 &
                                            PXD000279_peptide$NA_rate < 1),]
#high intensity, medium missing peptides
PXD000279_himm <- PXD000279_peptide[which(PXD000279_peptide$mean_intensity >= PXD000279_int_high & 
                                            PXD000279_peptide$NA_rate >= 0.35 &
                                            PXD000279_peptide$NA_rate < 0.65),]
#high intensity, low missing peptides
PXD000279_hilm <- PXD000279_peptide[which(PXD000279_peptide$mean_intensity >= PXD000279_int_high & 
                                            PXD000279_peptide$NA_rate < 0.35),]
#medium intensity, high missing peptides
PXD000279_mihm <- PXD000279_peptide[which(PXD000279_peptide$mean_intensity < PXD000279_int_high & 
                                            PXD000279_peptide$mean_intensity >= PXD000279_int_low &
                                            PXD000279_peptide$NA_rate >= 0.65 &
                                            PXD000279_peptide$NA_rate < 1),]
#medium intensity, medium missing peptides
PXD000279_mimm <- PXD000279_peptide[which(PXD000279_peptide$mean_intensity < PXD000279_int_high & 
                                            PXD000279_peptide$mean_intensity >= PXD000279_int_low &
                                            PXD000279_peptide$NA_rate >= 0.35 &
                                            PXD000279_peptide$NA_rate < 0.65),]
#medium intensity, low missing peptides
PXD000279_milm <- PXD000279_peptide[which(PXD000279_peptide$mean_intensity < PXD000279_int_high & 
                                            PXD000279_peptide$mean_intensity >= PXD000279_int_low &
                                            PXD000279_peptide$NA_rate < 0.35),]
#low intensity, high missing peptides
PXD000279_lihm <- PXD000279_peptide[which(PXD000279_peptide$mean_intensity < PXD000279_int_low & 
                                            PXD000279_peptide$NA_rate >= 0.65 &
                                            PXD000279_peptide$NA_rate < 1),]
#low intensity, medium missing peptides
PXD000279_limm <- PXD000279_peptide[which(PXD000279_peptide$mean_intensity < PXD000279_int_low & 
                                            PXD000279_peptide$NA_rate >= 0.35 &
                                            PXD000279_peptide$NA_rate < 0.65),]
#low intensity, low missing peptides
PXD000279_lilm <- PXD000279_peptide[which(PXD000279_peptide$mean_intensity < PXD000279_int_low & 
                                            PXD000279_peptide$NA_rate < 0.35),]

#keep intensity
PXD000279_hihm <- PXD000279_hihm[,-c(1,2,3)]
PXD000279_himm <- PXD000279_himm[,-c(1,2,3)]
PXD000279_hilm <- PXD000279_hilm[,-c(1,2,3)]
PXD000279_mihm <- PXD000279_mihm[,-c(1,2,3)]
PXD000279_mimm <- PXD000279_mimm[,-c(1,2,3)]
PXD000279_milm <- PXD000279_milm[,-c(1,2,3)]
PXD000279_lihm <- PXD000279_lihm[,-c(1,2,3)]
PXD000279_limm <- PXD000279_limm[,-c(1,2,3)]
PXD000279_lilm <- PXD000279_lilm[,-c(1,2,3)]

#prepare datasets for imputation, remove sample and pepetide that are all NA
PXD000279_hihm_na <- PXD000279_hihm[which(rowSums(is.na(PXD000279_hihm)) < ncol(PXD000279_hihm)),
                                    which(colSums(is.na(PXD000279_hihm)) < nrow(PXD000279_hihm))]
PXD000279_himm_na <- PXD000279_himm[which(rowSums(is.na(PXD000279_himm)) < ncol(PXD000279_himm)),
                                    which(colSums(is.na(PXD000279_himm)) < nrow(PXD000279_himm))]
PXD000279_hilm_na <- PXD000279_hilm[which(rowSums(is.na(PXD000279_hilm)) < ncol(PXD000279_hilm)),
                                    which(colSums(is.na(PXD000279_hilm)) < nrow(PXD000279_hilm))]
PXD000279_mihm_na <- PXD000279_mihm[which(rowSums(is.na(PXD000279_mihm)) < ncol(PXD000279_mihm)),
                                    which(colSums(is.na(PXD000279_mihm)) < nrow(PXD000279_mihm))]
PXD000279_mimm_na <- PXD000279_mimm[which(rowSums(is.na(PXD000279_mimm)) < ncol(PXD000279_mimm)),
                                    which(colSums(is.na(PXD000279_mimm)) < nrow(PXD000279_mimm))]
PXD000279_milm_na <- PXD000279_milm[which(rowSums(is.na(PXD000279_milm)) < ncol(PXD000279_milm)),
                                    which(colSums(is.na(PXD000279_milm)) < nrow(PXD000279_milm))]
PXD000279_lihm_na <- PXD000279_lihm[which(rowSums(is.na(PXD000279_lihm)) < ncol(PXD000279_lihm)),
                                    which(colSums(is.na(PXD000279_lihm)) < nrow(PXD000279_lihm))]
PXD000279_limm_na <- PXD000279_limm[which(rowSums(is.na(PXD000279_limm)) < ncol(PXD000279_limm)),
                                    which(colSums(is.na(PXD000279_limm)) < nrow(PXD000279_limm))]
PXD000279_lilm_na <- PXD000279_lilm[which(rowSums(is.na(PXD000279_lilm)) < ncol(PXD000279_lilm)),
                                    which(colSums(is.na(PXD000279_lilm)) < nrow(PXD000279_lilm))]

#real datasets
PXD000279_hihm_real <- PXD000279_hihm_na
PXD000279_himm_real <- PXD000279_himm_na
PXD000279_hilm_real <- PXD000279_hilm_na
PXD000279_mihm_real <- PXD000279_mihm_na
PXD000279_mimm_real <- PXD000279_mimm_na
PXD000279_milm_real <- PXD000279_milm_na
PXD000279_lihm_real <- PXD000279_lihm_na
PXD000279_limm_real <- PXD000279_limm_na
PXD000279_lilm_real <- PXD000279_lilm_na

#generate random sample points, 10%
which(is.na(PXD000279_hihm_na) == FALSE)

na_pos_PXD000279_hihm <- sample(which(is.na(PXD000279_hihm_na) == FALSE), 
                                0.1*length(which(is.na(PXD000279_hihm_na) == FALSE)))
na_pos_PXD000279_himm <- sample(which(is.na(PXD000279_himm_na) == FALSE), 
                                0.1*length(which(is.na(PXD000279_himm_na) == FALSE)))
na_pos_PXD000279_hilm <- sample(which(is.na(PXD000279_hilm_na) == FALSE), 
                                0.1*length(which(is.na(PXD000279_hilm_na) == FALSE)))
na_pos_PXD000279_mihm <- sample(which(is.na(PXD000279_mihm_na) == FALSE), 
                                0.1*length(which(is.na(PXD000279_mihm_na) == FALSE)))
na_pos_PXD000279_mimm <- sample(which(is.na(PXD000279_mimm_na) == FALSE), 
                                0.1*length(which(is.na(PXD000279_mimm_na) == FALSE)))
na_pos_PXD000279_milm <- sample(which(is.na(PXD000279_milm_na) == FALSE), 
                                0.1*length(which(is.na(PXD000279_milm_na) == FALSE)))
na_pos_PXD000279_lihm <- sample(which(is.na(PXD000279_lihm_na) == FALSE), 
                                0.1*length(which(is.na(PXD000279_lihm_na) == FALSE)))
na_pos_PXD000279_limm <- sample(which(is.na(PXD000279_limm_na) == FALSE), 
                                0.1*length(which(is.na(PXD000279_limm_na) == FALSE)))
na_pos_PXD000279_lilm <- sample(which(is.na(PXD000279_lilm_na) == FALSE), 
                                0.1*length(which(is.na(PXD000279_lilm_na) == FALSE)))

#replace true values with NA
PXD000279_hihm_na <- as.matrix(PXD000279_hihm_na)
PXD000279_hihm_na[na_pos_PXD000279_hihm] <- NA
PXD000279_hihm_na <- as.data.frame(PXD000279_hihm_na)
PXD000279_himm_na <- as.matrix(PXD000279_himm_na)
PXD000279_himm_na[na_pos_PXD000279_himm] <- NA
PXD000279_himm_na <- as.data.frame(PXD000279_himm_na)
PXD000279_hilm_na <- as.matrix(PXD000279_hilm_na)
PXD000279_hilm_na[na_pos_PXD000279_hilm] <- NA
PXD000279_hilm_na <- as.data.frame(PXD000279_hilm_na)
PXD000279_mihm_na <- as.matrix(PXD000279_mihm_na)
PXD000279_mihm_na[na_pos_PXD000279_mihm] <- NA
PXD000279_mihm_na <- as.data.frame(PXD000279_mihm_na)
PXD000279_mimm_na <- as.matrix(PXD000279_mimm_na)
PXD000279_mimm_na[na_pos_PXD000279_mimm] <- NA
PXD000279_mimm_na <- as.data.frame(PXD000279_mimm_na)
PXD000279_milm_na <- as.matrix(PXD000279_milm_na)
PXD000279_milm_na[na_pos_PXD000279_milm] <- NA
PXD000279_milm_na <- as.data.frame(PXD000279_milm_na)
PXD000279_lihm_na <- as.matrix(PXD000279_lihm_na)
PXD000279_lihm_na[na_pos_PXD000279_lihm] <- NA
PXD000279_lihm_na <- as.data.frame(PXD000279_lihm_na)
PXD000279_limm_na <- as.matrix(PXD000279_limm_na)
PXD000279_limm_na[na_pos_PXD000279_limm] <- NA
PXD000279_limm_na <- as.data.frame(PXD000279_limm_na)
PXD000279_lilm_na <- as.matrix(PXD000279_lilm_na)
PXD000279_lilm_na[na_pos_PXD000279_lilm] <- NA
PXD000279_lilm_na <- as.data.frame(PXD000279_lilm_na)


#change na pos: hihm, mihm, lihm
#count na percentage
which(rowSums(is.na(PXD000279_hihm_na)) == ncol(PXD000279_hihm_na))
which(rowSums(is.na(PXD000279_himm_na)) == ncol(PXD000279_himm_na))
which(rowSums(is.na(PXD000279_hilm_na)) == ncol(PXD000279_hilm_na))
which(rowSums(is.na(PXD000279_mihm_na)) == ncol(PXD000279_mihm_na))
which(rowSums(is.na(PXD000279_mimm_na)) == ncol(PXD000279_mimm_na))
which(rowSums(is.na(PXD000279_milm_na)) == ncol(PXD000279_milm_na))
which(rowSums(is.na(PXD000279_lihm_na)) == ncol(PXD000279_lihm_na))
which(rowSums(is.na(PXD000279_limm_na)) == ncol(PXD000279_limm_na))
which(rowSums(is.na(PXD000279_lilm_na)) == ncol(PXD000279_lilm_na))


#impute
id <- "PXD000279"
for (i in 1:length(regions)) {
  temp <- get(paste("id", regions[i], "na", sep = "_"))
  
  #bpca
  #training set
  temp <- pca(temp, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
  temp <- as.data.frame(temp@completeObs)
  assign(paste(regions[i], "bpca", sep = "_"), temp)
  print(paste("id", regions[i], "bpca complete", sep = " "))
  
  #svd
  #training set
  temp <- as.data.frame(impute.wrapper.SVD(temp, K = 5))
  assign(paste(regions[i], "svd", sep = "_"), temp)
  print(paste("id", regions[i], "svd complete", sep = " "))
  
  #knn
  #trainging set
  temp <- as.data.frame(impute.wrapper.KNN(as.matrix(temp), K = 5))
  assign(paste(regions[i], "knn", sep = "_"), temp)
  print(paste("id", regions[i], "knn complete", sep = " "))
  
  #mle
  temp <- as.data.frame(impute.wrapper.MLE(as.matrix(temp)))
  assign(paste(regions[i], "mle", sep = "_"), temp)
  print(paste("id", regions[i], "mle complete", sep = " "))
  
  #lls
  temp <- llsImpute(temp, k = 5, center = TRUE, completeObs = TRUE)
  temp <- as.data.frame(temp@completeObs)
  assign(paste(regions[i], "lls", sep = "_"), temp)
  print(paste("id", regions[i], "lls complete", sep = " "))
  
  #rf
  temp <- missForest(temp, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
  temp <- temp$ximp
  assign(paste(regions[i], "rf", sep = "_"), temp)
  print(paste("id", regions[i], "rf complete", sep = " "))
}


#evaluation NRMSE
id <- "PXD000279"
regions <- c("hihm", "himm", "hilm", "mihm", "mimm", "milm", "lihm", "limm", "lilm")
methods <- c("bpca", "knn", "lls", "mle", "svd", "rf", "cf", "dae", "vae")

for (i in 1:length(regions)) {
  #get na data
  comp_na <- get(paste(id, regions[i], "na", sep = "_"))
  #get real data
  comp_real <- get(paste(id, regions[i], "real", sep = "_"))
  comp_real <- comp_real[which(rownames(comp_real) %in% rownames(comp_na)),]
  
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
  for (j in 1:length(methods)) {
    if (paste(id, regions[i], methods[j], sep = "_") %in% ls() == TRUE) {
      temp <- get(paste(id, regions[i], methods[j], sep = "_"))
      for (k in 1:length(nrmse_sum$sample)) {
        nrmse_sum[k, grep(methods[j], colnames(nrmse_sum))] <- sqrt(mean((temp[,k] - comp_real[,k])^2, na.rm = TRUE))/sd(comp_real[,k], na.rm = TRUE)
      }
    } else {
      next
    }
    #if highest value > 4 fold of median, remove values higher than 90% values in a method
    if (max(nrmse_sum[,grep(methods[j], colnames(nrmse_sum))]) > 4*quantile(nrmse_sum[,grep(methods[j], colnames(nrmse_sum))], probs = c(0,0.5,0.9,0.95,1))[2]) {
      print(paste(regions[i], methods[j]))
      print(quantile(nrmse_sum[,grep(methods[j], colnames(nrmse_sum))], probs = c(0,0.5,0.9,0.95,1)))
      nrmse_sum[which(nrmse_sum[,grep(methods[j], colnames(nrmse_sum))] > 5*quantile(nrmse_sum[,grep(methods[j], colnames(nrmse_sum))], probs = c(0,0.5,0.9,0.95,1))[2]),
                grep(methods[j], colnames(nrmse_sum))] <- NA
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
  
  if (max(nrmse_stat$mean, na.rm = TRUE) > 10) {
    ggplot(data = nrmse_stat, aes(x = method, y = mean, fill = method)) +
      geom_bar(stat = "identity") +
      geom_errorbar(aes(ymin = mean - 1.96*se, ymax = mean + 1.96*se, width = 0.2)) + #95% confidence interval
      geom_text(aes(label = round(mean, digits = 2)), vjust = -2.5, size = 5) +
      scale_fill_brewer(palette = "Set3", labels = c("RF", "BPCA", "kNN", "LLS", "MLE", "SVD", "CF", "DAE", "VAE")) +
      scale_y_break(c(1.6, 200), scales = 0.5) +
      labs(x =  toupper(regions[i]),
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
      geom_text(aes(label = round(mean, digits = 2)), vjust = -3, size = 5) +
      scale_fill_brewer(palette = "Set3", labels = c("RF", "BPCA", "kNN", "LLS", "MLE", "SVD", "CF", "DAE", "VAE")) +
      ylim(0,1.6) +
      scale_y_continuous(labels = label_number(accuracy = 0.1)) +
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
  ggsave(paste("figures/manuscript/PXD000279/NRMSE_barplot_", regions[i], ".jpeg", sep = ""), device = "jpeg", width = 1400, height = 1400, units = "px", dpi = 300)
  print(regions[i])
  print(colMeans(nrmse_sum[,2:length(nrmse_sum)], na.rm = TRUE))
  quants <- c(0, 0.25, 0.50, 0.75, 1)
  print(apply(nrmse_sum[,2:length(nrmse_sum)], 2 , quantile , probs = quants , na.rm = TRUE))
  print(paste(gsub("nrmse_", "", nrmse_stat$method[which(nrmse_stat$mean==min(nrmse_stat$mean, na.rm = TRUE))]),
              "is the imputation method with lowest NRMSE in",
              regions[i],
              sep = " "))
}


##evaluation differential expression analysis
PXD000279_meta <- data.frame("sample" = colnames(PXD000279_whole_real),
                             "dose" = substr(str_split_i(colnames(PXD000279_whole_real), " ", 2), 1,1))

DE_result <- data.frame("df" = ls()[grep("PXD000279_whole_", ls())],
                        "method" = NA,
                        "sig_DEP" = NA,
                        "sig_DF" = NA,
                        "DF" = NA)
DE_result$method <- str_split_i(DE_result$df, "_", 3)

#limma
design <- model.matrix(~ 0 + dose, data = PXD000279_meta)
for (i in 1:length(DE_result$df)) {
  data <- get(DE_result$df[i])
  rid <- rownames(data)
  cid <- colnames(data)
  data <- normalize.quantiles(as.matrix(data))
  rownames(data) <- rid
  colnames(data) <- cid
  fit <- lmFit(data, design)
  con <- makeContrasts(doseH - doseL, levels = design)
  fit.con <- contrasts.fit(fit, con)
  bayes <- eBayes(fit.con)
  result <- topTable(bayes, sort.by="none", n=Inf, adjust.method = "BH")
  sig <- result[result$adj.P.Val < 0.05 & abs(result$logFC) > 1.5,] #1 or 1.5
  sig <- sig[!is.na(sig$logFC),]
  dfname <- paste(DE_result$df[i], "sigDEP", sep = "_")
  assign(dfname, sig)
  DE_result$sig_DEP[i] <- length(rownames(sig))
  DE_result$sig_DF[i] <- dfname
  dfres <- paste(DE_result$df[i], "result", sep = "_")
  assign(dfres, result)
  DE_result$DF[i] <- dfres
  print(paste(DE_result$df[i], length(sig$logFC), "sigDEP", sep = " "))
}

DE_result <- DE_result[DE_result$method != "na",]
DE_result$method <- factor(DE_result$method, levels = c("real", "rf", "bpca", "knn", "lls", "mle", 
                                                        "svd", "cf", "dae", "vae", "mix"))

#map species to peptide
ppmap$species <- NA
for (i in 1:length(ppmap$peptide)) {
  pos1 <- grep(ppmap$protein[i], PXD000279_proteinGroups$`Protein IDs`)
  ppmap$species[i] <- ifelse(grepl("HUMAN", PXD000279_proteinGroups$`Fasta headers`[pos1]) == TRUE,
                             "Human",
                             "E.coli")
}
table(ppmap$species)

#ROC
roc_plot <- data.frame(peptide = rownames(PXD000279_whole_real_result),
                       truth = NA,
                       real = NA,
                       knn = NA,
                       mle = NA,
                       svd = NA,
                       bpca = NA,
                       lls = NA,
                       rf = NA,
                       cf = NA,
                       vae = NA,
                       dae = NA,
                       mix = NA)

for (i in 1:length(roc_plot$peptide)) {
  
  roc_plot$truth[i] <- ifelse(roc_plot$peptide[i] %in% ppmap$peptide[ppmap$species == "E.coli"], 1, 0)
  
}
id <- "PXD000279"
for (i in 3:length(colnames(roc_plot))) {
  temp <- get(paste(id, "whole", colnames(roc_plot)[i], "result", sep = "_"))
  for (j in 1:nrow(roc_plot)) {
    roc_plot[j,i] <- ifelse(temp$adj.P.Val[j] < 0.05, 1, 0)
  }
}

roc_plot_long <- melt_roc(roc_plot, "truth", colnames(roc_plot)[-c(1:2)])
ggplot(roc_plot_long, aes(d = D, m = M, color = name)) + 
  geom_roc() + 
  style_roc() 

roc_data <- roc_plot_long %>%
  group_by(name) %>%
  do({
    roc_obj <- roc(.$D, .$M)
    data.frame(
      FPR = 1 - roc_obj$specificities,
      TPR = roc_obj$sensitivities,
      name = unique(.$name)
    )
  })
roc_data <- roc_data[-which(roc_data$name=="real"),]
roc_data$name <- factor(roc_data$name, levels = c("rf", "bpca", "knn", "lls", "mle", 
                                                  "svd", "cf", "dae", "vae", "mix"))

ggplot(roc_data, aes(x = FPR, y = TPR, color = name)) +
  #geom_line(size = 1) +
  geom_smooth(method = "loess", se = FALSE, size = 1) +
  labs(title = "ROC Curve for Each Method", x = "False Positive Rate", y = "True Positive Rate", color = "Method") +
  theme_minimal() +
  theme(legend.title = element_blank())+  
  style_roc() +
  scale_color_manual(values = c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462","#b3de69","#fccde5","#d9d9d9","#bc80bd"),
                     labels =  c("RF", "BPCA", "kNN", "LLS", "MLE", "SVD", "CF", "DAE", "VAE", "Mix")) +
  theme(axis.text = element_text(size = 15, color = "black"),
        axis.title = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 13, color = "black"),
        lengend.title = element_text(size = 18, color = "black"))

#get auc
auc(roc_plot$truth, roc_plot$mix)
auc(roc_plot$truth, roc_plot$rf)
auc(roc_plot$truth, roc_plot$bpca)
auc(roc_plot$truth, roc_plot$knn)
auc(roc_plot$truth, roc_plot$mle)
auc(roc_plot$truth, roc_plot$svd)
auc(roc_plot$truth, roc_plot$lls)
auc(roc_plot$truth, roc_plot$cf)
auc(roc_plot$truth, roc_plot$vae)
auc(roc_plot$truth, roc_plot$dae)
