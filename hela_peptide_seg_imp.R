setwd("C:/Users/ys25.stu/OneDrive - UBC/BeeCSI/Proteomics/imputation")

library(readr)
library(tibble)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(matrixStats)
library(factoextra) #pca
library(imputeLCMD) #knn, svd, MLE(not working properly?) imputation
library(pcaMethods) #bpca, lls imputation
library(purrr)
library(svglite)
library(gprofiler2)
library(matrixStats)
library(xlsx)
library(VIM) #KNN
library(pROC)
library(ROCR)
library(missForest)
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
library(rlang)

#import data
hela <- read.delim("combined_modified_peptide_MBR.tsv")
hela <- hela[,-c(1, 3:7,9:10, 12:16)]

hela_maxLFQ_intensity <- hela[, c(1:3, grep("MaxLFQ", colnames(hela)))]
hela_maxLFQ_intensity <- add_column(hela_maxLFQ_intensity, mean_intensity = NA, .after = "Protein.ID")
hela_maxLFQ_intensity <- add_column(hela_maxLFQ_intensity, NA_count = NA, .after = "mean_intensity")
hela_maxLFQ_intensity <- add_column(hela_maxLFQ_intensity, NA_rate = NA, .after = "NA_count")

hela_maxLFQ_intensity$mean_intensity <- rowMeans(hela_maxLFQ_intensity[, -c(1:6)])
hela_maxLFQ_intensity$NA_count <- rowSums(hela_maxLFQ_intensity[, -c(1:6)] == 0)
hela_maxLFQ_intensity$NA_rate <- hela_maxLFQ_intensity$NA_count/(ncol(hela_maxLFQ_intensity) - 6)

int_low <- log2(quantile(hela_maxLFQ_intensity$mean_intensity)[2])
int_high <- log2(quantile(hela_maxLFQ_intensity$mean_intensity)[4])

#log intensity vs NA rate
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
ggsave("figure1_int_mr_hela.tiff",
       plot = last_plot(),
       path = "figures/manuscript/",
       device = "tiff",
       width = 2244,
       height = 1496,
       units = "px",
       dpi = 300)
summary(lm(NA_rate~log2(mean_intensity), data = hela_maxLFQ_intensity[-which(hela_maxLFQ_intensity$NA_rate == 1),]))

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



#donut plot
group_stats_plot_donut_new <- read.csv("donut_plot_new.csv")
group_stats_plot_donut_new[1,2:10] <- c(269, 1594, 25704, 11943, 23957, 19234, 20942, 577, 0)

list_im <- c("hihm", "himm", "hilm", "mihm", "mimm", "milm", "lihm", "limm", "lilm")
for (i in 1:length(list_im)) {
  temp_data <- group_stats_plot_donut_new[c(11, grep(list_im[i], colnames(group_stats_plot_donut_new)))]
  temp_name <- paste("donut_plot", list_im[i], sep = "_")
  temp_plot <- ggplot(temp_data[2:4,],
                      aes(ymax = temp_data[2:4,3], ymin = temp_data[2:4,4],
                          xmax = 4, xmin = 3, fill = temp_data[2:4,1])) +
    geom_rect() +
    geom_text(x = 2, 
              aes(y = temp_data[2:4,5], 
                  label = paste(round((temp_data[2:4,2] * 100), 2), "%", sep = "") ), 
              size = 10) +
    geom_text(x = 0, y = 0, hjust = 0.5, vjust = 0.5,
              label = paste(temp_data[1,2], toupper(list_im[i]), sep = "\n"), size = 15) +
    coord_polar(theta = "y") +
    xlim(c(-1, 4)) +
    theme_void() +
    scale_fill_brewer(palette = "Set1") +
    labs(fill = "Spectrum Match Type")
  #theme(panel.background = element_rect(fill='transparent'),
  #     plot.background = element_rect(fill='transparent', color = NA),
  #    panel.grid.major = element_blank(),
  #   panel.grid.minor = element_blank(),
  #  legend.background = element_rect(fill='transparent'),
  # legend.box.background = element_rect(fill='transparent'))
  print(temp_plot)
  ggsave(paste("piechart_", list_im[i], ".svg", sep = ""),
         plot = last_plot(),
         path = "figures/manuscript/piechart/",
         device = "svg",
         width = 3000,
         height = 3000,
         units = "px",
         dpi = 300)
  cat("plot", temp_name, "succesfully saved \n")
}



##get retention time, ion mobility, observed precursor m/z from psm
#read in psm files
for (i in 1:length(dir("./psm"))) {
  temp_name <- gsub(".tsv", "", dir("./psm")[i])
  temp <- read_delim(paste("./psm", dir("./psm")[i], sep = "/"), delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  assign(temp_name, temp)
}

psm_list_rt <- gsub(".tsv", "", dir("./psm")) %>% paste(., "RT", sep = "_")
psm_list_omz <- gsub(".tsv", "", dir("./psm")) %>% paste(., "ob.m.z", sep = "_")
psm_list_id <- gsub(".tsv", "", dir("./psm")) %>% gsub("psm", "", .)

hihm[c(psm_list_rt, psm_list_omz)] <- NA
himm[c(psm_list_rt, psm_list_omz)] <- NA
hilm[c(psm_list_rt, psm_list_omz)] <- NA
mihm[c(psm_list_rt, psm_list_omz)] <- NA
mimm[c(psm_list_rt, psm_list_omz)] <- NA
milm[c(psm_list_rt, psm_list_omz)] <- NA
lihm[c(psm_list_rt, psm_list_omz)] <- NA
limm[c(psm_list_rt, psm_list_omz)] <- NA
lilm[c(psm_list_rt, psm_list_omz)] <- NA

grep(paste(psm_list_id, collapse = "|"), colnames(hihm))

hihm <- hihm[,c(1:8, grep(paste(psm_list_id, collapse = "|"), colnames(hihm)))]
himm <- himm[,c(1:8, grep(paste(psm_list_id, collapse = "|"), colnames(himm)))]
hilm <- hilm[,c(1:8, grep(paste(psm_list_id, collapse = "|"), colnames(hilm)))]
mihm <- mihm[,c(1:8, grep(paste(psm_list_id, collapse = "|"), colnames(mihm)))]
mimm <- mimm[,c(1:8, grep(paste(psm_list_id, collapse = "|"), colnames(mimm)))]
milm <- milm[,c(1:8, grep(paste(psm_list_id, collapse = "|"), colnames(milm)))]
lihm <- lihm[,c(1:8, grep(paste(psm_list_id, collapse = "|"), colnames(lihm)))]
limm <- limm[,c(1:8, grep(paste(psm_list_id, collapse = "|"), colnames(limm)))]
lilm <- lilm[,c(1:8, grep(paste(psm_list_id, collapse = "|"), colnames(lilm)))]

#extract retention time, observed m/z
for (i in 1:length(list_im)) {
  temp <- get(list_im[i])
  print(list_im[i])
  for (j in 1:length(psm_list_id)) {
    #load psm
    temp_psm <- get(paste("psm", psm_list_id[j], sep = ""))
    #intensity
    pos1 <- grep(paste(psm_list_id[j], "MaxLFQ", sep = "."), colnames(temp))
    #retention time
    pos2 <- grep(paste(psm_list_id[j], "RT", sep = "_"), colnames(temp))
    #observed m/z
    pos3 <- grep(paste(psm_list_id[j], "ob.m.z", sep = "_"), colnames(temp))
    for (k in 1:length(temp$Modified.Sequence)) {
      if (temp[k, pos1] > 0) {
        pos4 <- which(temp_psm$Peptide == temp$Modified.Sequence[k])
        if (length(pos4) == 0) {
          temp[k, pos2] <- "MBR"
          temp[k, pos3] <- "MBR"
        } else if (length(pos4) == 1) {
          #retention time
          temp[k, pos2] <- temp_psm$Retention[pos4]
          #observed m/z
          temp[k, pos3] <- temp_psm$`Observed M/Z`[pos4]
        } else {
          #combined retention time
          temp[k, pos2] <- paste(temp_psm$Retention[pos4], collapse = ", ")
          #combined observed m/z
          temp[k, pos3] <- paste(temp_psm$`Observed M/Z`[pos4], collapse = ", ")
        }
      } else {
        temp[k, pos2] <- "unmatched"
        temp[k, pos3] <- "unmatched"
      } 
    }
  }
  assign(list_im[i], temp)
}

#remove 11333, 11837
hihm <- hihm[,-grep(paste(c("11333", "11837"), collapse = "|"), colnames(hihm))]




#MBR percentage in each group
group_stats <- data.frame(group = list_im,
                          peptide = NA,
                          MSMS_rate_mean = NA,
                          MSMS_rate_sd = NA,
                          MBR_rate_mean = NA,
                          MBR_rate_sd = NA,
                          unmatched_rate_mean = NA,
                          unmatched_rate_sd = NA)
group_stats[paste("MSMS", psm_list_id, sep = "_")] <- NA
group_stats[paste("MBR", psm_list_id, sep = "_")] <- NA
group_stats[paste("unmatched", psm_list_id, sep = "_")] <- NA

stats <- c("MSMS", "MBR", "unmatched")

for (i in 1:length(group_stats$group)) {
  temp <- get(group_stats$group[i])
  group_stats$peptide[i] <- length(temp$Modified.Sequence)
  pos1 <- which(group_stats$group == group_stats$group[i])
  for (j in 1:length(psm_list_id)) {
    pos2 <- grep(paste(psm_list_id[j], "RT", sep = "_"), colnames(temp))
    for (k in 1:length(stats)) {
      if (stats[k] == "MSMS") {
        group_stats[pos1, colnames(group_stats) == paste(stats[k], psm_list_id[j], sep = "_")] <- length(which(temp[,pos2] != "MBR" & temp[,pos2] != "unmatched"))/group_stats$peptide[i]
      } else {
        pos3 <- grep(stats[k], temp[,pos2])
        group_stats[pos1, colnames(group_stats) == paste(stats[k], psm_list_id[j], sep = "_")] <- length(pos3)/group_stats$peptide[i]
      }
    }
  }
}

group_stats$MSMS_rate_mean <- rowMeans(group_stats[,grep("MSMS", colnames(group_stats))][,-c(1:2)])
group_stats$MSMS_rate_sd <- rowSds(as.matrix(group_stats[,grep("MSMS", colnames(group_stats))][,-c(1:2)]))
group_stats$MBR_rate_mean <- rowMeans(group_stats[,grep("MBR", colnames(group_stats))][,-c(1:2)])
group_stats$MBR_rate_sd <- rowSds(as.matrix(group_stats[,grep("MBR", colnames(group_stats))][,-c(1:2)]))
group_stats$unmatched_rate_mean <- rowMeans(group_stats[,grep("unmatched", colnames(group_stats))][,-c(1:2)])
group_stats$unmatched_rate_sd <- rowSds(as.matrix(group_stats[,grep("unmatched", colnames(group_stats))][,-c(1:2)]))

group_stats_plot <- group_stats[,-c(2:8)] %>% pivot_longer(cols = 2:58,
                                                           names_to = "name",
                                                           values_to = "value")
group_stats_plot <- separate_wider_delim(group_stats_plot, cols = name, delim = "_", names = c("stats", "sample"))

#x axis order
group_stats_plot$group <- factor(group_stats_plot$group, levels = list_im)
group_stats_plot$stats <- factor(group_stats_plot$stats, levels = c("MSMS", "MBR", "unmatched"))

#boxplot
ggplot(data = group_stats_plot, aes(x = group, y = value, fill = stats)) +
  geom_boxplot() +
  geom_jitter(color="black", size = 0.6, alpha = 0.9) +
  facet_wrap(~stats, scales = "fixed") +
  xlab("Group") +
  ylab("Rate") +
  labs(fill = "Spectrum Match Type") +
  theme(axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 14),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(face = "bold", size = 13),
        strip.text = element_text(face = "bold", size = 13),
        panel.grid = element_blank(),
        panel.background = element_blank())

#donut chart
group_stats_plot_donut <- as.data.frame(t(group_stats[,1:8]))
colnames(group_stats_plot_donut) <- group_stats_plot_donut[1,]
group_stats_plot_donut <- group_stats_plot_donut[-c(1,4,6,8),]
group_stats_plot_donut$type <- rownames(group_stats_plot_donut)
group_stats_plot_donut[,1:9] <- lapply(group_stats_plot_donut[,1:9], as.numeric)
group_stats_plot_donut[,paste("ymax", list_im, sep = "_")] <- NA
group_stats_plot_donut[,paste("ymin", list_im, sep = "_")] <- NA

group_stats_plot_donut[2:4,11:19] <- lapply(group_stats_plot_donut[2:4,1:9], cumsum)
group_stats_plot_donut[2:4,20:28] <- lapply(group_stats_plot_donut[2:4,1:9], cumsum)
for (i in 20:28) {
  group_stats_plot_donut[2:4, i] <- c(0, group_stats_plot_donut[2:3, i-9])
}

group_stats_plot_donut$type <- c("peptide amount", "MS/MS", "MBR", "Unmatched")

group_stats_plot_donut[,paste("labelposition", list_im, sep = "_")] <- NA
for (i in 29:37) {
  group_stats_plot_donut[2:4, i] <- (group_stats_plot_donut[2:4, i-18] + group_stats_plot_donut[2:4, i-9])/2
}

ggplot(group_stats_plot_donut[2:4, c(1, 11, 10, 20, 29)],
       aes(ymax = ymax_hihm, ymin = ymin_hihm,
           xmax = 4, xmin = 3, fill = type)) +
  geom_rect() +
  geom_text(x = 2.4, 
            aes(y = labelposition_hihm, 
                label = paste(round((hihm * 100), 2), "%", sep = "") ), 
            size = 6) +
  geom_text(x = 0, y = 0, label = "1573", size = 10) +
  coord_polar(theta = "y") +
  xlim(c(-1, 4)) +
  theme_void() +
  scale_fill_brewer(palette = "Set1") +
  labs(fill = "Spectrum Match Type")

list_im <- list_im[-9]

group_stats_plot_donut$type <- factor(group_stats_plot_donut$type, levels = c("MS/MS", "MBR", "Unmatched"))

for (i in 1:length(list_im)) {
  temp_data <- group_stats_plot_donut[c(10, grep(list_im[i], colnames(group_stats_plot_donut)))]
  temp_name <- paste("donut_plot", list_im[i], sep = "_")
  svg(paste("./figures/", temp_name, ".svg", sep = ""))
  temp_plot <- ggplot(temp_data[2:4,],
                      aes(ymax = temp_data[2:4,3], ymin = temp_data[2:4,4],
                          xmax = 4, xmin = 3, fill = temp_data[2:4,1])) +
    geom_rect() +
    geom_text(x = 2, 
              aes(y = temp_data[2:4,5], 
                  label = paste(round((temp_data[2:4,2] * 100), 2), "%", sep = "") ), 
              size = 6) +
    geom_text(x = 0, y = 0, hjust = 0.5, #vjust = 0.5,
              label = paste(temp_data[1,2], toupper(list_im[i]), sep = "\n"), size = 10) +
    coord_polar(theta = "y") +
    
    xlim(c(-1, 4)) +
    theme_void() +
    scale_fill_brewer(palette = "Set1") +
    labs(fill = "Spectrum Match Type")
  #theme(panel.background = element_rect(fill='transparent'),
  #     plot.background = element_rect(fill='transparent', color = NA),
  #    panel.grid.major = element_blank(),
  #   panel.grid.minor = element_blank(),
  #  legend.background = element_rect(fill='transparent'),
  # legend.box.background = element_rect(fill='transparent'))
  print(temp_plot)
  dev.off()
  cat("plot", temp_name, "succesfully saved \n")
}



regions <- c("hihm", "himm", "hilm", "mihm", "mimm", "milm", "lihm", "limm")
methods <- c("bpca", "knn", "lls", "mle", "svd", "rf", "cf", "dae", "vae")

for (i in 1:length(regions)) {
  #get na data
  comp_na <- get(paste(regions[i], "na", sep = "_"))
  #get real data
  comp_real <- get(paste(regions[i], "real", sep = "_"))
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
    if (paste(regions[i], methods[j], sep = "_") %in% ls() == TRUE) {
      temp <- get(paste(regions[i], methods[j], sep = "_"))
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
      geom_text(aes(label = round(mean, digits = 2)), vjust = -3, size = 5) +
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
      geom_text(aes(label = round(mean, digits = 2)), vjust = -3, size = 5) +
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
  ggsave(paste("figures/manuscript/Hela/NRMSE_barplot_", regions[i], ".jpeg", sep = ""), device = "jpeg", width = 1200, height = 1200, units = "px", dpi = 300)
  print(regions[i])
  print(colMeans(nrmse_sum[,2:length(nrmse_sum)], na.rm = TRUE))
  quants <- c(0, 0.25, 0.50, 0.75, 1)
  print(apply(nrmse_sum[,2:length(nrmse_sum)], 2 , quantile , probs = quants , na.rm = TRUE))
  print(paste(gsub("nrmse_", "", nrmse_stat$method[which(nrmse_stat$mean==min(nrmse_stat$mean, na.rm = TRUE))]),
              "is the imputation method with lowest NRMSE in",
              regions[i],
              sep = " "))
}












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

#pimms
#cf, vae, dae
himm_trs_cf <- read.csv("himm_trs_cf.csv")

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

regions <- c("hihm", "himm", "hilm", "mihm", "mimm", "milm", "lihm", "limm")
datatype <- c("trs", "tes")
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

expr_print(pretty10exp(ax, drop.1=TRUE))

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



#wilxocon rank test

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

library(stamats)
library(exactRankTests)
wilcox.exact(comp$real, comp$mix, paired = TRUE, conf.int = TRUE)
mi.wilcox.test(comp, x = "bpca", y = "mix", paired = TRUE)

wilcox.exact(comp$real, comp$rf, paired = TRUE, conf.int = TRUE)

res_wilcox <- data.frame(method1 = rep(colnames(comp)[c(2,4:13)], ncol(comp)-1),
                         method2 = rep(colnames(comp)[c(2,4:13)], each = ncol(comp)-1),
                         v = NA,
                         pvalue = NA)

for (i in 1:nrow(res_wilcox)) {
  temp1 <- comp[,grep(res_wilcox$method1[i], colnames(comp))]
  temp2 <- comp[,grep(res_wilcox$method2[i], colnames(comp))]
  
  if (res_wilcox$method1[i] != res_wilcox$method2[i]) {
    res <- wilcox.exact(temp1, temp2, paired = TRUE, conf.int = TRUE)
    res_wilcox$v[i] <- res$statistic
    res_wilcox$pvalue[i] <- res$p.value
  } else {
    next
  }
  if (i%%4 == 0) {
    print(i)
  }
}

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


pheatmap(res_cor[-c(3,9,10),-c(3,9,10)],
         brewer.pal(n = 8, name = "PuRd"))

set.seed(123)
x <- rnorm(25, mean = 1)
x[sample(1:25, 5)] <- NA
y <- rnorm(20, mean = -1)
y[sample(1:20, 4)] <- NA
pair <- c(rnorm(25, mean = 1), rnorm(20, mean = -1))
g <- factor(c(rep("yes", 25), rep("no", 20)))
D <- data.frame(ID = 1:45, response = c(x, y), pair = pair, group = g)




#test all possible combination
methods <- c("bpca", "knn", "lls", "mle", "svd", "rf", "cf", "dae", "vae")
regions <- c("hihm", "himm", "hilm", "mihm", "mimm", "milm", "lihm", "limm")

# Create all combinations of methods
imp_comb <- expand.grid(rep(list(methods), length(regions)))
# Assign column names from regions
colnames(imp_comb) <- regions

#add NRMSE
imp_comb$nrmse_mean <- NA
imp_comb$nrmse_sd <- NA

system.time(
  for (i in 1:500) {
    #assemble testing dataset
    for (j in 1:(ncol(imp_comb)-2)) {
      temp_name <- paste(colnames(imp_comb)[j], "trs", imp_comb[i,j], sep = "_")
      if (temp_name %in% ls() == TRUE) {
        temp1 <- get(temp_name)
      } else {
        next
      }
      if (j == 1) {
        temp <- temp1
      } else {
        temp <- rbind(temp, temp1)
      }
    }
    
    temp_nrmse <- data.frame(sample = colnames(temp),
                             nrmse = NA)
    for (k in 1:ncol(temp)) {
      temp_nrmse$nrmse[k] <-  sqrt(mean((temp[,k] - whole_trs_real[,k])^2, na.rm = TRUE))/sd(whole_trs_real[,k], na.rm = TRUE)
    }
    
    imp_comb$nrmse_mean[i] <- mean(temp_nrmse$nrmse)
    imp_comb$nrmse_sd[i] <- sd(temp_nrmse$nrmse)
  }
)
#about 10 seconds to calculate NRMSE for 20 combinations
#for 9^8 = 43046721 combinations, it will need 5721 hours = 238 days

parallel::detectCores(logical = TRUE)  # Returns the number of logical cores
parallel::detectCores(logical = FALSE) # Returns the number of physical cores



# Set up parallel backend
num_cores <- 18  
cl <- makeCluster(num_cores)
registerDoParallel(cl)
# Parallel execution
system.time(
  imp_comb_results <- foreach(i = 1:500, .combine = rbind) %dopar% {
    temp <- data.frame()
    for (j in 1:(ncol(imp_comb)-2)) {
      temp_name <- paste(colnames(imp_comb)[j], "trs", imp_comb[i,j], sep = "_")
      if (temp_name %in% ls()) {
        temp1 <- get(temp_name)
        if (nrow(temp) == 0) {
          temp <- temp1
        } else {
          temp <- rbind(temp, temp1)
        }
      }
    }
    
    temp_nrmse <- data.frame(sample = colnames(temp),
                             nrmse = sapply(1:ncol(temp), function(k) {
                               sqrt(mean((temp[,k] - whole_trs_real[,k])^2, na.rm = TRUE)) / 
                                 sd(whole_trs_real[,k], na.rm = TRUE)
                             }))
    
    c(nrmse_mean = mean(temp_nrmse$nrmse),
      nrmse_sd = sd(temp_nrmse$nrmse))
  }
)
# Stop the cluster
stopCluster(cl)
#with multi-threading using 18 cores, 344 combinations took 53.25 seconds
#test all combinations will take 1842 hours



#real datasets
hihm_real <- hihm_na
himm_real <- himm_na
hilm_real <- hilm_na
mihm_real <- mihm_na
mimm_real <- mimm_na
milm_real <- milm_na
lihm_real <- lihm_na
limm_real <- limm_na

#generate random sample points, 20%
na_pos_hihm <- sample(length(!hihm_na), 0.2*length(!hihm_na))
na_pos_himm <- sample(length(!himm_na), 0.2*length(!himm_na))
na_pos_hilm <- sample(length(!hilm_na), 0.2*length(!hilm_na))
na_pos_mihm <- sample(length(!mihm_na), 0.2*length(!mihm_na))
na_pos_mimm <- sample(length(!mimm_na), 0.2*length(!mimm_na))
na_pos_milm <- sample(length(!milm_na), 0.2*length(!milm_na))
na_pos_lihm <- sample(length(!lihm_na), 0.2*length(!lihm_na))
na_pos_limm <- sample(length(!limm_na), 0.2*length(!limm_na))

#replace true values with NA
hihm_na <- as.matrix(hihm_na)
hihm_na[na_pos_hihm] <- NA
hihm_na <- as.data.frame(hihm_na)
himm_na <- as.matrix(himm_na)
himm_na[na_pos_himm] <- NA
himm_na <- as.data.frame(himm_na)
hilm_na <- as.matrix(hilm_na)
hilm_na[na_pos_hilm] <- NA
hilm_na <- as.data.frame(hilm_na)
mihm_na <- as.matrix(mihm_na)
mihm_na[na_pos_mihm] <- NA
mihm_na <- as.data.frame(mihm_na)
mimm_na <- as.matrix(mimm_na)
mimm_na[na_pos_mimm] <- NA
mimm_na <- as.data.frame(mimm_na)
milm_na <- as.matrix(milm_na)
milm_na[na_pos_milm] <- NA
milm_na <- as.data.frame(milm_na)
lihm_na <- as.matrix(lihm_na)
lihm_na[na_pos_lihm] <- NA
lihm_na <- as.data.frame(lihm_na)
limm_na <- as.matrix(limm_na)
limm_na[na_pos_limm] <- NA
limm_na <- as.data.frame(limm_na)

#count na percentage
hihm_na <- hihm_na[-which(rowSums(is.na(hihm_na)) == ncol(hihm_na)),]
#himm_na <- himm_na[-which(rowSums(is.na(himm_na)) == ncol(himm_na)),]
#hilm_na <- hilm_na[-which(rowSums(is.na(hilm_na)) == ncol(hilm_na)),]
mihm_na <- mihm_na[-which(rowSums(is.na(mihm_na)) == ncol(mihm_na)),]
#mimm_na <- mimm_na[-which(rowSums(is.na(mimm_na)) == ncol(mimm_na)),]
#milm_na <- milm_na[-which(rowSums(is.na(milm_na)) == ncol(milm_na)),]
lihm_na <- lihm_na[-which(rowSums(is.na(lihm_na)) == ncol(lihm_na)),]
#limm_na <- limm_na[-which(rowSums(is.na(limm_na)) == ncol(limm_na)),]


#impute by intensity and missing rate
#bpca, svd, knn (SEQKNN, TRKNN, KNNIMPUTE), median, missForrest, RSN
#bpca
hihm_bpca <- pca(hihm_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
himm_bpca <- pca(himm_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
hilm_bpca <- pca(hilm_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
mihm_bpca <- pca(mihm_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
mimm_bpca <- pca(mimm_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
milm_bpca <- pca(milm_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
lihm_bpca <- pca(lihm_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
limm_bpca <- pca(limm_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
hihm_bpca <- as.data.frame(hihm_bpca@completeObs)
himm_bpca <- as.data.frame(himm_bpca@completeObs)
hilm_bpca <- as.data.frame(hilm_bpca@completeObs)
mihm_bpca <- as.data.frame(mihm_bpca@completeObs)
mimm_bpca <- as.data.frame(mimm_bpca@completeObs)
milm_bpca <- as.data.frame(milm_bpca@completeObs)
lihm_bpca <- as.data.frame(lihm_bpca@completeObs)
limm_bpca <- as.data.frame(limm_bpca@completeObs)
#svd
hihm_svd <- as.data.frame(impute.wrapper.SVD(hihm_na, K = 10))
himm_svd <- as.data.frame(impute.wrapper.SVD(himm_na, K = 10))
hilm_svd <- as.data.frame(impute.wrapper.SVD(hilm_na, K = 10))
mihm_svd <- as.data.frame(impute.wrapper.SVD(mihm_na, K = 10))
mimm_svd <- as.data.frame(impute.wrapper.SVD(mimm_na, K = 10))
milm_svd <- as.data.frame(impute.wrapper.SVD(milm_na, K = 10))
lihm_svd <- as.data.frame(impute.wrapper.SVD(lihm_na, K = 10))
limm_svd <- as.data.frame(impute.wrapper.SVD(limm_na, K = 10))
#knn
hihm_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(hihm_na), K = 5))
himm_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(himm_na), K = 5))
hilm_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(hilm_na), K = 5))
mihm_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(mihm_na), K = 5))
mimm_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(mimm_na), K = 5))
milm_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(milm_na), K = 5))
lihm_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(lihm_na), K = 5))
limm_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(limm_na), K = 5))
#mle
hihm_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(hihm_na)))
himm_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(himm_na)))
hilm_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(hilm_na)))
mihm_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(mihm_na)))
mimm_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(mimm_na)))
milm_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(milm_na)))
lihm_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(lihm_na)))
limm_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(limm_na)))
#lls
hihm_lls <- llsImpute(hihm_na, k = 5, center = TRUE, completeObs = TRUE)
himm_lls <- llsImpute(himm_na, k = 5, center = TRUE, completeObs = TRUE)
hilm_lls <- llsImpute(hilm_na, k = 5, center = TRUE, completeObs = TRUE)
mihm_lls <- llsImpute(mihm_na, k = 5, center = TRUE, completeObs = TRUE)
mimm_lls <- llsImpute(mimm_na, k = 5, center = TRUE, completeObs = TRUE)
milm_lls <- llsImpute(milm_na, k = 5, center = TRUE, completeObs = TRUE)
lihm_lls <- llsImpute(lihm_na, k = 5, center = TRUE, completeObs = TRUE)
limm_lls <- llsImpute(limm_na, k = 5, center = TRUE, completeObs = TRUE)
hihm_lls <- as.data.frame(hihm_lls@completeObs)
himm_lls <- as.data.frame(himm_lls@completeObs)
hilm_lls <- as.data.frame(hilm_lls@completeObs)
mihm_lls <- as.data.frame(mihm_lls@completeObs)
mimm_lls <- as.data.frame(mimm_lls@completeObs)
milm_lls <- as.data.frame(milm_lls@completeObs)
lihm_lls <- as.data.frame(lihm_lls@completeObs)
limm_lls <- as.data.frame(limm_lls@completeObs)
#rf
#register number of cores for parallel
registerDoParallel(cores = 18)
hihm_rf <- missForest(hihm_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
himm_rf <- missForest(himm_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
hilm_rf <- missForest(hilm_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
mihm_rf <- missForest(mihm_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
mimm_rf <- missForest(mimm_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
milm_rf <- missForest(milm_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
lihm_rf <- missForest(lihm_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
limm_rf <- missForest(limm_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
hihm_rf <- hihm_rf$ximp
himm_rf <- himm_rf$ximp
hilm_rf <- hilm_rf$ximp
mihm_rf <- mihm_rf$ximp
mimm_rf <- mimm_rf$ximp
milm_rf <- milm_rf$ximp
lihm_rf <- lihm_rf$ximp
limm_rf <- limm_rf$ximp

#transpose data frame for pimms
write.csv(t(hihm_na), file = "hihm_na_80.csv")
write.csv(t(himm_na), file = "himm_na_80.csv")
write.csv(t(hilm_na), file = "hilm_na_80.csv")
write.csv(t(mihm_na), file = "mihm_na_80.csv")
write.csv(t(mimm_na), file = "mimm_na_80.csv")
write.csv(t(milm_na), file = "milm_na_80.csv")
write.csv(t(lihm_na), file = "lihm_na_80.csv")
write.csv(t(limm_na), file = "limm_na_80.csv")

#import pimms result
#cf
himm_cf <- read.csv("himm_80_cf.csv")
hilm_cf <- read.csv("hilm_80_cf.csv")
mihm_cf <- read.csv("mihm_80_cf.csv")
mimm_cf <- read.csv("mimm_80_cf.csv")
milm_cf <- read.csv("milm_80_cf.csv")
lihm_cf <- read.csv("lihm_80_cf.csv")
limm_cf <- read.csv("limm_80_cf.csv")
#vae
hihm_vae <- read.csv("hihm_80_vae.csv")
himm_vae <- read.csv("himm_80_vae.csv")
hilm_vae <- read.csv("hilm_80_vae.csv")
mihm_vae <- read.csv("mihm_80_vae.csv")
mimm_vae <- read.csv("mimm_80_vae.csv")
milm_vae <- read.csv("milm_80_vae.csv")
lihm_vae <- read.csv("lihm_80_vae.csv")
limm_vae <- read.csv("limm_80_vae.csv")
#dae
hihm_dae <- read.csv("hihm_80_dae.csv")
himm_dae <- read.csv("himm_80_dae.csv")
hilm_dae <- read.csv("hilm_80_dae.csv")
mihm_dae <- read.csv("mihm_80_dae.csv")
mimm_dae <- read.csv("mimm_80_dae.csv")
milm_dae <- read.csv("milm_80_dae.csv")
lihm_dae <- read.csv("lihm_80_dae.csv")
limm_dae <- read.csv("limm_80_dae.csv")

#set colname
colnames(hilm_dae) <- hilm_dae[1,]
colnames(mimm_dae) <- mimm_dae[1,]
colnames(milm_dae) <- milm_dae[1,]
colnames(lihm_dae) <- lihm_dae[1,]
colnames(hilm_vae) <- hilm_vae[1,]
colnames(mimm_vae) <- mimm_vae[1,]
colnames(milm_vae) <- milm_vae[1,]
colnames(lihm_vae) <- lihm_vae[1,]

#remove excessive rows
hilm_dae <- hilm_dae[-c(1,2),]
mimm_dae <- mimm_dae[-c(1,2),]
milm_dae <- milm_dae[-c(1,2),]
lihm_dae <- lihm_dae[-c(1,2),]
hilm_vae <- hilm_vae[-c(1,2),]
mimm_vae <- mimm_vae[-c(1,2),]
milm_vae <- milm_vae[-c(1,2),]
lihm_vae <- lihm_vae[-c(1,2),]

#turn chr into num
hilm_dae[,-1] <- sapply(hilm_dae[,-1], as.numeric)
mimm_dae[,-1] <- sapply(mimm_dae[,-1], as.numeric)
milm_dae[,-1] <- sapply(milm_dae[,-1], as.numeric)
lihm_dae[,-1] <- sapply(lihm_dae[,-1], as.numeric)
hilm_vae[,-1] <- sapply(hilm_vae[,-1], as.numeric)
mimm_vae[,-1] <- sapply(mimm_vae[,-1], as.numeric)
milm_vae[,-1] <- sapply(milm_vae[,-1], as.numeric)
lihm_vae[,-1] <- sapply(lihm_vae[,-1], as.numeric)

#set sample id as row name
#cf
rownames(himm_cf) <- himm_cf[,1]
rownames(hilm_cf) <- hilm_cf[,1]
rownames(mihm_cf) <- mihm_cf[,1]
rownames(mimm_cf) <- mimm_cf[,1]
rownames(milm_cf) <- milm_cf[,1]
rownames(lihm_cf) <- lihm_cf[,1]
rownames(limm_cf) <- limm_cf[,1]
#vae
rownames(hihm_vae) <- hihm_vae[,1]
rownames(himm_vae) <- himm_vae[,1]
rownames(hilm_vae) <- hilm_vae[,1]
rownames(mihm_vae) <- mihm_vae[,1]
rownames(mimm_vae) <- mimm_vae[,1]
rownames(milm_vae) <- milm_vae[,1]
rownames(lihm_vae) <- lihm_vae[,1]
rownames(limm_vae) <- limm_vae[,1]
#dae
rownames(hihm_dae) <- hihm_dae[,1]
rownames(himm_dae) <- himm_dae[,1]
rownames(hilm_dae) <- hilm_dae[,1]
rownames(mihm_dae) <- mihm_dae[,1]
rownames(mimm_dae) <- mimm_dae[,1]
rownames(milm_dae) <- milm_dae[,1]
rownames(lihm_dae) <- lihm_dae[,1]
rownames(limm_dae) <- limm_dae[,1]

#transpose
#cf
himm_cf <- as.data.frame(t(himm_cf[,-1]))
hilm_cf <- as.data.frame(t(hilm_cf[,-1]))
mihm_cf <- as.data.frame(t(mihm_cf[,-1]))
mimm_cf <- as.data.frame(t(mimm_cf[,-1]))
milm_cf <- as.data.frame(t(milm_cf[,-1]))
lihm_cf <- as.data.frame(t(lihm_cf[,-1]))
limm_cf <- as.data.frame(t(limm_cf[,-1]))
#vae
hihm_vae <- as.data.frame(t(hihm_vae[,-1]))
himm_vae <- as.data.frame(t(himm_vae[,-1]))
hilm_vae <- as.data.frame(t(hilm_vae[,-1]))
mihm_vae <- as.data.frame(t(mihm_vae[,-1]))
mimm_vae <- as.data.frame(t(mimm_vae[,-1]))
milm_vae <- as.data.frame(t(milm_vae[,-1]))
lihm_vae <- as.data.frame(t(lihm_vae[,-1]))
limm_vae <- as.data.frame(t(limm_vae[,-1]))
#dae
hihm_dae <- as.data.frame(t(hihm_dae[,-1]))
himm_dae <- as.data.frame(t(himm_dae[,-1]))
hilm_dae <- as.data.frame(t(hilm_dae[,-1]))
mihm_dae <- as.data.frame(t(mihm_dae[,-1]))
mimm_dae <- as.data.frame(t(mimm_dae[,-1]))
milm_dae <- as.data.frame(t(milm_dae[,-1]))
lihm_dae <- as.data.frame(t(lihm_dae[,-1]))
limm_dae <- as.data.frame(t(limm_dae[,-1]))


#evaluation by region
id <- ""
regions <- c("hihm", "himm", "hilm", "mihm", "mimm", "milm", "lihm", "limm", "lilm")
methods <- c("bpca", "knn", "lls", "mle", "svd", "rf", "cf", "dae", "vae")

for (i in 1:length(regions)) {
  #get na data
  comp_na <- get(paste(regions[i], "na", sep = "_"))
  #get real data
  comp_real <- get(paste(regions[i], "real", sep = "_"))
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
      geom_text(aes(label = round(mean, digits = 2)), vjust = -2, size = 5) +
      scale_fill_discrete(labels = c("RF", "BPCA", "kNN", "LLS", "MLE", "SVD", "CF", "DAE", "VAE")) +
      scale_y_break(c(1.6, 200), scales = 0.5) +
      labs(x = regions[i],
           y = "NRMSE") +
      theme_bw() +
      theme(axis.title = element_text(face = "bold", size = 12),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = 10),
            panel.grid.major.x = element_blank(),
            panel.grid.major = element_line(linewidth = 1),
            legend.position = "none")
  } else {
    ggplot(data = nrmse_stat, aes(x = method, y = mean, fill = method)) +
      geom_bar(stat = "identity") +
      geom_errorbar(aes(ymin = mean - 1.96*se, ymax = mean + 1.96*se, width = 0.2)) + #95% confidence interval
      geom_text(aes(label = round(mean, digits = 2)), vjust = -2, size = 5) +
      scale_fill_discrete(labels = c("RF", "BPCA", "kNN", "LLS", "MLE", "SVD", "CF", "DAE", "VAE")) +
      ylim(0,1.6) +
      labs(x = regions[i],
           y = "NRMSE") +
      theme_bw() +
      theme(axis.title = element_text(face = "bold", size = 12),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = 10),
            panel.grid.major.x = element_blank(),
            panel.grid.major = element_line(linewidth = 1),
            legend.position = "none")
  }
  ggsave(paste("figures/pimms/segregation/PXD000279/NRMSE_barplot_", regions[i], ".png", sep = ""), device = "png", width = 1500, height = 1500, units = "px")
  print(regions[i])
  print(colMeans(nrmse_sum[,2:length(nrmse_sum)], na.rm = TRUE))
  quants <- c(0, 0.25, 0.50, 0.75, 1)
  print(apply(nrmse_sum[,2:length(nrmse_sum)], 2 , quantile , probs = quants , na.rm = TRUE))
  print(paste(gsub("nrmse_", "", nrmse_stat$method[which(nrmse_stat$mean==min(nrmse_stat$mean, na.rm = TRUE))]),
              "is the imputation method with lowest NRMSE in",
              regions[i],
              sep = " "))
}











#hihm
comp_real <- hihm_real[which(rownames(hihm_real) %in% rownames(hihm_na)),] %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "real")
comp_na <- hihm_na %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "na")
comp_knn <- hihm_knn %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "knn")
comp_mle <- hihm_mle %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "mle")
comp_svd <- hihm_svd %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "svd")
comp_bpca <- hihm_bpca %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "bpca")
comp_lls <- hihm_lls %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "lls")
comp_rf <- hihm_rf %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "rf")
comp_vae <- hihm_vae %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "vae")
comp_dae <- hihm_dae %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "dae")

comp <- data.frame(sample = comp_real$sample,
                   real = comp_real$real,
                   na = comp_na$na,
                   knn = comp_knn$knn,
                   mle = comp_mle$mle,
                   svd = comp_svd$svd,
                   bpca = comp_bpca$bpca,
                   lls = comp_lls$lls,
                   rf = comp_rf$rf,
                   vae = comp_vae$vae,
                   dae = comp_dae$dae)
comp_long <- comp %>% pivot_longer(., 2:length(.), names_to = "method", values_to = "value")

#NRMSE
nrmse_sum <- data.frame(sample = unique(comp_real$sample),
                        nrmse_knn = NA,
                        nrmse_mle = NA,
                        nrmse_svd = NA, 
                        nrmse_bpca = NA, 
                        nrmse_lls = NA,
                        nrmse_rf = NA,
                        nrmse_cf = NA,
                        nrmse_vae = NA,
                        nrmse_dae = NA)

for (i in 1:length(nrmse_sum$sample)) {
  z <- hihm_real[which(rownames(hihm_real) %in% rownames(hihm_na) & rownames(hihm_real) %in% rownames(hihm_vae)),grep(nrmse_sum$sample[i], colnames(hihm_real))]
  #normalized by standard deviation of zero imputed matrix
  nrmse_sum$nrmse_knn[i] <- sqrt(mean((hihm_knn[which(rownames(hihm_knn) %in% rownames(hihm_vae)),grep(nrmse_sum$sample[i], colnames(hihm_knn))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_mle[i] <- sqrt(mean((hihm_mle[rownames(hihm_mle) %in% rownames(hihm_vae),grep(nrmse_sum$sample[i], colnames(hihm_mle))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_svd[i] <- sqrt(mean((hihm_svd[rownames(hihm_svd) %in% rownames(hihm_vae),grep(nrmse_sum$sample[i], colnames(hihm_svd))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_bpca[i] <- sqrt(mean((hihm_bpca[rownames(hihm_bpca) %in% rownames(hihm_vae),grep(nrmse_sum$sample[i], colnames(hihm_bpca))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_lls[i] <- sqrt(mean((hihm_lls[rownames(hihm_lls) %in% rownames(hihm_vae),grep(nrmse_sum$sample[i], colnames(hihm_lls))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_rf[i] <- sqrt(mean((hihm_rf[rownames(hihm_rf) %in% rownames(hihm_vae),grep(nrmse_sum$sample[i], colnames(hihm_rf))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_vae[i] <- sqrt(mean((hihm_vae[rownames(hihm_vae) %in% rownames(hihm_vae),grep(nrmse_sum$sample[i], colnames(hihm_vae))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_dae[i] <- sqrt(mean((hihm_dae[rownames(hihm_dae) %in% rownames(hihm_vae),grep(nrmse_sum$sample[i], colnames(hihm_dae))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
}

nrmse_long <- nrmse_sum %>% pivot_longer(cols = 2:length(.), values_to = "NRMSE", names_to = "method")
nrmse_long$method <- factor(nrmse_long$method, levels = c("nrmse_bpca", "nrmse_knn", "nrmse_lls", "nrmse_mle", 
                                                          "nrmse_svd", "nrmse_rf", "nrmse_cf", "nrmse_dae", "nrmse_vae"))
#boxplot
ggplot(data = nrmse_long, aes(x = method, y = NRMSE)) +
  geom_boxplot() +
  #geom_jitter(shape = 16, position = position_jitter(0.2)) +
  labs(title = "NRMSE of different missing value imputation methods \nHIHM",
       x = "Imputation Methods") +
  scale_x_discrete(labels = c("BPCA", "kNN", "LLS", "MLE", "SVD", "RF", "CF", "DAE", "VAE")) +
  theme_bw() +
  theme(axis.text.x = element_text(face = "bold", size = 10), 
        axis.title = element_text(size = 10),
        title = element_text(size = 12)) +
  ylim(0,2)

#barplot
#calculate mean and standard error
#drop mle
nrmse_stat <- data.frame(method = colnames(nrmse_sum)[-c(1,3)],
                         mean = colMeans(nrmse_sum[,-c(1,3)]),
                         se = colSds(as.matrix(nrmse_sum[,-c(1,3)]))/sqrt(nrow(nrmse_sum)))
nrmse_stat$method <- factor(nrmse_stat$method, levels = c("nrmse_rf", "nrmse_bpca", "nrmse_knn", "nrmse_lls",  
                                                          "nrmse_svd", "nrmse_cf", "nrmse_dae", "nrmse_vae"))
ggplot(data = nrmse_stat, aes(x = method, y = mean, fill = method)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = mean - 1.96*se, ymax = mean + 1.96*se, width = 0.2)) + #95% confidence interval
  geom_text(aes(label = round(mean, digits = 2)), vjust = -5, size = 5) +
  scale_fill_discrete(labels = c("RF", "BPCA", "kNN", "LLS", "SVD", "CF", "DAE", "VAE")) +
  ylim(0, 1.6) +
  labs(x = "HIHM",
       y = "NRMSE") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        panel.grid.major.x = element_blank(),
        panel.grid.major = element_line(linewidth = 1),
        legend.position = "none")
colMeans(nrmse_sum[,2:length(nrmse_sum)])
quants <- c(0, 0.25, 0.50, 0.75, 1)
apply(nrmse_sum[,2:length(nrmse_sum)], 2 , quantile , probs = quants , na.rm = TRUE )


#get legend
bap_limm <- ggplot(data = nrmse_stat, aes(x = method, y = mean, fill = method)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = mean - 1.96*se, ymax = mean + 1.96*se, width = 0.2)) + #95% confidence interval
  geom_text(aes(label = round(mean, digits = 2)), vjust = -1, size = 5) +
  scale_fill_discrete(labels = c("RF", "BPCA", "kNN", "LLS", "SVD", "CF", "DAE", "VAE")) +
  ylim(0, 1.6) +
  labs(x = "LIMM",
       y = "NRMSE") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        panel.grid.major.x = element_blank(),
        panel.grid.major = element_line(linewidth = 1),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12, face = "bold"))

legend <- get_legend(bap_limm)
grid.newpage()
grid.draw(legend)


#impute as whole
whole_na <- rbind(hihm_na, himm_na, hilm_na,
                  mihm_na, mimm_na, milm_na,
                  lihm_na, limm_na)
#real
whole_real <- rbind(hihm_real[rownames(hihm_real) %in% rownames(hihm_na),],
                    himm_real[rownames(himm_real) %in% rownames(himm_na),],
                    hilm_real[rownames(hilm_real) %in% rownames(hilm_na),],
                    mihm_real[rownames(mihm_real) %in% rownames(mihm_na),],
                    mimm_real[rownames(mimm_real) %in% rownames(mimm_na),],
                    milm_real[rownames(milm_real) %in% rownames(milm_na),],
                    lihm_real[rownames(lihm_real) %in% rownames(lihm_na),],
                    limm_real[rownames(limm_real) %in% rownames(limm_na),])
#bpca
whole_bpca <- pca(whole_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
whole_bpca <- as.data.frame(whole_bpca@completeObs)
#svd
whole_svd <- as.data.frame(impute.wrapper.SVD(whole_na, K = 10))
#knn
whole_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(whole_na), K = 5))
#mle
whole_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(whole_na)))
#lls
whole_lls <- llsImpute(whole_na, k = 5, center = TRUE, completeObs = TRUE)
whole_lls <- as.data.frame(whole_lls@completeObs)
#rf
whole_rf <- missForest(whole_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
whole_rf <- whole_rf$ximp

#for pimms
write.csv(t(whole_na), file = "whole_na_80.csv")
#import pimms result
#fread for fast read
whole_cf <- fread("whole_80_cf.csv", header = TRUE)
whole_dae <- fread("whole_80_dae.csv", header = TRUE)
whole_vae <- fread("whole_80_vae.csv", header = TRUE)

whole_cf <- as.data.frame(whole_cf)
whole_dae <- as.data.frame(whole_dae)
whole_vae <- as.data.frame(whole_vae)

#set colname
colnames(whole_dae) <- whole_dae[1,]
colnames(whole_vae) <- whole_vae[1,]

#remove excessive rows
whole_dae <- whole_dae[-c(1,2),]
whole_vae <- whole_vae[-c(1,2),]

#turn chr into num
whole_dae[,-1] <- sapply(whole_dae[,-1], as.numeric)
whole_vae[,-1] <- sapply(whole_vae[,-1], as.numeric)

#set sample id as row name
rownames(whole_cf) <- whole_cf[,1]
rownames(whole_dae) <- whole_dae[,1]
rownames(whole_vae) <- whole_vae[,1]

#transpose
whole_cf <- as.data.frame(t(whole_cf[,-1]))
whole_dae <- as.data.frame(t(whole_dae[,-1]))
whole_vae <- as.data.frame(t(whole_vae[,-1]))

#assemble the optimal method matrix
whole_mix <- rbind(hihm_rf, himm_rf, hilm_rf,
                   mihm_bpca, mimm_rf, milm_rf,
                   lihm_vae, limm_vae)


#get na data
comp_na <- whole_na
#get real data
comp_real <- whole_real[which(rownames(whole_real) %in% rownames(comp_na)),]

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

methods <- c("bpca", "knn", "lls", "mle", "svd", "rf", "cf", "dae", "vae", "mix")
#get data for different methods
for (j in 1:length(methods)) {
  if (paste("whole", methods[j], sep = "_") %in% ls() == TRUE) {
    temp <- get(paste("whole", methods[j], sep = "_"))
    for (k in 1:length(nrmse_sum$sample)) {
      nrmse_sum[k, grep(methods[j], colnames(nrmse_sum))] <- sqrt(mean((temp[,k] - comp_real[,k])^2, na.rm = TRUE))/sd(comp_real[,k], na.rm = TRUE)
    }
  } else {
    next
  }
  #if highest value > 4 fold of median, remove values higher than 5*(50% values) in a method
  if (max(nrmse_sum[,grep(methods[j], colnames(nrmse_sum))]) > 4*quantile(nrmse_sum[,grep(methods[j], colnames(nrmse_sum))], probs = c(0,0.5,0.9,0.95,1))[2]) {
    print(methods[j])
    print(quantile(nrmse_sum[,grep(methods[j], colnames(nrmse_sum))], probs = c(0,0.5,0.9,0.95,1)))
    nrmse_sum[which(nrmse_sum[,grep(methods[j], colnames(nrmse_sum))] > 3*quantile(nrmse_sum[,grep(methods[j], colnames(nrmse_sum))], probs = c(0,0.5,0.9,0.95,1))[2]),
              grep(methods[j], colnames(nrmse_sum))] <- NA
  }
}


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
        axis.text = element_text(size = 15, color = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.major = element_line(linewidth = 1),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15))

ggsave("figure2B_mix_compare_PXD000279.tiff",
       plot = last_plot(),
       path = "figures/manuscript/",
       device = "tiff",
       width = 2244,
       height = 950,
       units = "px",
       dpi = 300)
ggsave("figure2B_mix_compare_PXD000279.jpeg",
       plot = last_plot(),
       path = "figures/manuscript/",
       device = "jpeg",
       width = 2244,
       height = 950,
       units = "px",
       dpi = 300)



##validation data
#PXD029525
PXD029525 <- read_delim("reference/Pietz. 2024/data/datasets/PXD029525_prostate_cancer/peptides.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
PXD029525_meta <- read_delim("reference/Pietz. 2024/data/datasets/PXD029525_prostate_cancer/experimentalDesignTemplate.proteome.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

#keep intensity
PXD029525_peptide <- as.data.frame(PXD029525[,c(1,119:136)])
#set row name
rownames(PXD029525_peptide) <- PXD029525_peptide[,1]
PXD029525_peptide[,-1] <- log2(PXD029525_peptide[,-1])
PXD029525_peptide <- PXD029525_peptide[,-1]
PXD029525_peptide[PXD029525_peptide == -Inf] <- NA

#calculate NA rate, mean intensity
PXD029525_peptide <- add_column(PXD029525_peptide, mean_intensity = NA, .before = "LFQ intensity C42 CTRL rep01")
PXD029525_peptide <- add_column(PXD029525_peptide, NA_rate = NA, .after = "mean_intensity")

PXD029525_peptide$mean_intensity <- rowMeans(PXD029525_peptide[, -c(1,2)], na.rm = TRUE)
PXD029525_peptide$NA_rate <- rowSums(is.na(PXD029525_peptide[,-c(1,2)]))/(ncol(PXD029525_peptide) - 2)

#NA rate vs intensity
ggplot(data = PXD029525_peptide, aes(x = mean_intensity, y = NA_rate)) + 
  geom_point(alpha = 0.1, shape = 20) +
  labs(title = paste("PXD029525 Prostate cancer", "Relationship between peptide mean intensity and missing rate", sep = "\n"),
       x = bquote(paste("mean log"["2"]*"(intensity)")),
       y = "Missing rate") +
  stat_smooth(method = "lm", geom = "smooth", formula = y~x, se = TRUE, color = "black", size = 2) +
  lims(x = c(15,30),
       y = c(0,1)) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        title = element_text(size = 20)) 
last_plot() +
  stat_density2d(aes(fill = ..density..), geom = "tile", contour = FALSE, n=200, alpha = 0.8) +
  scale_fill_continuous(low = "white", high = "dodgerblue4") +
  geom_line(aes(x = PXD029525_int_high), color = "red", size = 1) +
  geom_line(aes(x = PXD029525_int_low), color = "red", size = 1) +
  geom_line(aes(y = 0.75), color = "red", size = 1) +
  geom_line(aes(y = 0.25), color = "red", size = 1) 




PXD029525_int_low <- quantile(PXD029525_peptide$mean_intensity, na.rm = TRUE)[2]
PXD029525_int_high <- quantile(PXD029525_peptide$mean_intensity, na.rm = TRUE)[4]

#high intensity, high missing peptides
PXD029525_hihm <- PXD029525_peptide[which(PXD029525_peptide$mean_intensity >= PXD029525_int_high & 
                                      PXD029525_peptide$NA_rate >= 0.75 &
                                      PXD029525_peptide$NA_rate < 1),]
#high intensity, medium missing peptides
PXD029525_himm <- PXD029525_peptide[which(PXD029525_peptide$mean_intensity >= PXD029525_int_high & 
                                      PXD029525_peptide$NA_rate >= 0.25 &
                                      PXD029525_peptide$NA_rate < 0.75),]
#high intensity, low missing peptides
PXD029525_hilm <- PXD029525_peptide[which(PXD029525_peptide$mean_intensity >= PXD029525_int_high & 
                                      PXD029525_peptide$NA_rate < 0.25),]
#medium intensity, high missing peptides
PXD029525_mihm <- PXD029525_peptide[which(PXD029525_peptide$mean_intensity < PXD029525_int_high & 
                                      PXD029525_peptide$mean_intensity >= PXD029525_int_low &
                                      PXD029525_peptide$NA_rate >= 0.75 &
                                      PXD029525_peptide$NA_rate < 1),]
#medium intensity, medium missing peptides
PXD029525_mimm <- PXD029525_peptide[which(PXD029525_peptide$mean_intensity < PXD029525_int_high & 
                                      PXD029525_peptide$mean_intensity >= PXD029525_int_low &
                                      PXD029525_peptide$NA_rate >= 0.25 &
                                      PXD029525_peptide$NA_rate < 0.75),]
#medium intensity, low missing peptides
PXD029525_milm <- PXD029525_peptide[which(PXD029525_peptide$mean_intensity < PXD029525_int_high & 
                                      PXD029525_peptide$mean_intensity >= PXD029525_int_low &
                                      PXD029525_peptide$NA_rate < 0.25),]
#low intensity, high missing peptides
PXD029525_lihm <- PXD029525_peptide[which(PXD029525_peptide$mean_intensity < PXD029525_int_low & 
                                      PXD029525_peptide$NA_rate >= 0.75 &
                                      PXD029525_peptide$NA_rate < 1),]
#low intensity, medium missing peptides
PXD029525_limm <- PXD029525_peptide[which(PXD029525_peptide$mean_intensity < PXD029525_int_low & 
                                      PXD029525_peptide$NA_rate >= 0.25 &
                                      PXD029525_peptide$NA_rate < 0.75),]
#low intensity, low missing peptides
PXD029525_lilm <- PXD029525_peptide[which(PXD029525_peptide$mean_intensity < PXD029525_int_low & 
                                      PXD029525_peptide$NA_rate < 0.25),]

#keep intensity
PXD029525_hihm <- PXD029525_hihm[,-c(1,2)]
PXD029525_himm <- PXD029525_himm[,-c(1,2)]
PXD029525_hilm <- PXD029525_hilm[,-c(1,2)]
PXD029525_mihm <- PXD029525_mihm[,-c(1,2)]
PXD029525_mimm <- PXD029525_mimm[,-c(1,2)]
PXD029525_milm <- PXD029525_milm[,-c(1,2)]
PXD029525_lihm <- PXD029525_lihm[,-c(1,2)]
PXD029525_limm <- PXD029525_limm[,-c(1,2)]
PXD029525_lilm <- PXD029525_lilm[,-c(1,2)]

#prepare datasets for imputation, remove sample and pepetide that are all NA
PXD029525_hihm_na <- PXD029525_hihm[which(rowSums(is.na(PXD029525_hihm)) < ncol(PXD029525_hihm)),
                                    which(colSums(is.na(PXD029525_hihm)) < nrow(PXD029525_hihm))]
PXD029525_himm_na <- PXD029525_himm[which(rowSums(is.na(PXD029525_himm)) < ncol(PXD029525_himm)),
                                    which(colSums(is.na(PXD029525_himm)) < nrow(PXD029525_himm))]
PXD029525_hilm_na <- PXD029525_hilm[which(rowSums(is.na(PXD029525_hilm)) < ncol(PXD029525_hilm)),
                                    which(colSums(is.na(PXD029525_hilm)) < nrow(PXD029525_hilm))]
PXD029525_mihm_na <- PXD029525_mihm[which(rowSums(is.na(PXD029525_mihm)) < ncol(PXD029525_mihm)),
                                    which(colSums(is.na(PXD029525_mihm)) < nrow(PXD029525_mihm))]
PXD029525_mimm_na <- PXD029525_mimm[which(rowSums(is.na(PXD029525_mimm)) < ncol(PXD029525_mimm)),
                                    which(colSums(is.na(PXD029525_mimm)) < nrow(PXD029525_mimm))]
PXD029525_milm_na <- PXD029525_milm[which(rowSums(is.na(PXD029525_milm)) < ncol(PXD029525_milm)),
                                    which(colSums(is.na(PXD029525_milm)) < nrow(PXD029525_milm))]
PXD029525_lihm_na <- PXD029525_lihm[which(rowSums(is.na(PXD029525_lihm)) < ncol(PXD029525_lihm)),
                                    which(colSums(is.na(PXD029525_lihm)) < nrow(PXD029525_lihm))]
PXD029525_limm_na <- PXD029525_limm[which(rowSums(is.na(PXD029525_limm)) < ncol(PXD029525_limm)),
                                    which(colSums(is.na(PXD029525_limm)) < nrow(PXD029525_limm))]
PXD029525_lilm_na <- PXD029525_lilm[which(rowSums(is.na(PXD029525_lilm)) < ncol(PXD029525_lilm)),
                                    which(colSums(is.na(PXD029525_lilm)) < nrow(PXD029525_lilm))]

#real datasets
PXD029525_hihm_real <- PXD029525_hihm_na
PXD029525_himm_real <- PXD029525_himm_na
PXD029525_hilm_real <- PXD029525_hilm_na
PXD029525_mihm_real <- PXD029525_mihm_na
PXD029525_mimm_real <- PXD029525_mimm_na
PXD029525_milm_real <- PXD029525_milm_na
PXD029525_lihm_real <- PXD029525_lihm_na
PXD029525_limm_real <- PXD029525_limm_na
PXD029525_lilm_real <- PXD029525_lilm_na

#generate random sample points, 20%
which(is.na(PXD029525_hihm_na) == FALSE)

na_pos_PXD029525_hihm <- sample(which(is.na(PXD029525_hihm_na) == FALSE), 
                                0.2*length(which(is.na(PXD029525_hihm_na) == FALSE)))
na_pos_PXD029525_himm <- sample(which(is.na(PXD029525_himm_na) == FALSE), 
                                0.2*length(which(is.na(PXD029525_himm_na) == FALSE)))
na_pos_PXD029525_hilm <- sample(which(is.na(PXD029525_hilm_na) == FALSE), 
                                0.2*length(which(is.na(PXD029525_hilm_na) == FALSE)))
na_pos_PXD029525_mihm <- sample(which(is.na(PXD029525_mihm_na) == FALSE), 
                                0.2*length(which(is.na(PXD029525_mihm_na) == FALSE)))
na_pos_PXD029525_mimm <- sample(which(is.na(PXD029525_mimm_na) == FALSE), 
                                0.2*length(which(is.na(PXD029525_mimm_na) == FALSE)))
na_pos_PXD029525_milm <- sample(which(is.na(PXD029525_milm_na) == FALSE), 
                                0.2*length(which(is.na(PXD029525_milm_na) == FALSE)))
na_pos_PXD029525_lihm <- sample(which(is.na(PXD029525_lihm_na) == FALSE), 
                                0.2*length(which(is.na(PXD029525_lihm_na) == FALSE)))
na_pos_PXD029525_limm <- sample(which(is.na(PXD029525_limm_na) == FALSE), 
                                0.2*length(which(is.na(PXD029525_limm_na) == FALSE)))
na_pos_PXD029525_lilm <- sample(which(is.na(PXD029525_lilm_na) == FALSE), 
                                0.2*length(which(is.na(PXD029525_lilm_na) == FALSE)))

#replace true values with NA
PXD029525_hihm_na <- as.matrix(PXD029525_hihm_na)
PXD029525_hihm_na[na_pos_PXD029525_hihm] <- NA
PXD029525_hihm_na <- as.data.frame(PXD029525_hihm_na)
PXD029525_himm_na <- as.matrix(PXD029525_himm_na)
PXD029525_himm_na[na_pos_PXD029525_himm] <- NA
PXD029525_himm_na <- as.data.frame(PXD029525_himm_na)
PXD029525_hilm_na <- as.matrix(PXD029525_hilm_na)
PXD029525_hilm_na[na_pos_PXD029525_hilm] <- NA
PXD029525_hilm_na <- as.data.frame(PXD029525_hilm_na)
PXD029525_mihm_na <- as.matrix(PXD029525_mihm_na)
PXD029525_mihm_na[na_pos_PXD029525_mihm] <- NA
PXD029525_mihm_na <- as.data.frame(PXD029525_mihm_na)
PXD029525_mimm_na <- as.matrix(PXD029525_mimm_na)
PXD029525_mimm_na[na_pos_PXD029525_mimm] <- NA
PXD029525_mimm_na <- as.data.frame(PXD029525_mimm_na)
PXD029525_milm_na <- as.matrix(PXD029525_milm_na)
PXD029525_milm_na[na_pos_PXD029525_milm] <- NA
PXD029525_milm_na <- as.data.frame(PXD029525_milm_na)
PXD029525_lihm_na <- as.matrix(PXD029525_lihm_na)
PXD029525_lihm_na[na_pos_PXD029525_lihm] <- NA
PXD029525_lihm_na <- as.data.frame(PXD029525_lihm_na)
PXD029525_limm_na <- as.matrix(PXD029525_limm_na)
PXD029525_limm_na[na_pos_PXD029525_limm] <- NA
PXD029525_limm_na <- as.data.frame(PXD029525_limm_na)
PXD029525_lilm_na <- as.matrix(PXD029525_lilm_na)
PXD029525_lilm_na[na_pos_PXD029525_lilm] <- NA
PXD029525_lilm_na <- as.data.frame(PXD029525_lilm_na)


####
#change na pos: hihm, mihm, mimm, lihm, limm
#count na percentage
PXD029525_hihm_na <- PXD029525_hihm_na[-which(rowSums(is.na(PXD029525_hihm_na)) == ncol(PXD029525_hihm_na)),]
PXD029525_mihm_na <- PXD029525_mihm_na[-which(rowSums(is.na(PXD029525_mihm_na)) == ncol(PXD029525_mihm_na)),]
PXD029525_mimm_na <- PXD029525_mimm_na[-which(rowSums(is.na(PXD029525_mimm_na)) == ncol(PXD029525_mimm_na)),]
PXD029525_lihm_na <- PXD029525_lihm_na[-which(rowSums(is.na(PXD029525_lihm_na)) == ncol(PXD029525_lihm_na)),]
PXD029525_limm_na <- PXD029525_limm_na[-which(rowSums(is.na(PXD029525_limm_na)) == ncol(PXD029525_limm_na)),]

#impute by intensity and missing rate
#bpca, svd, knn (SEQKNN, TRKNN, KNNIMPUTE), median, missForrest, RSN
#bpca
PXD029525_hihm_bpca <- pca(PXD029525_hihm_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
PXD029525_himm_bpca <- pca(PXD029525_himm_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
PXD029525_hilm_bpca <- pca(PXD029525_hilm_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
PXD029525_mihm_bpca <- pca(PXD029525_mihm_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
PXD029525_mimm_bpca <- pca(PXD029525_mimm_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
PXD029525_milm_bpca <- pca(PXD029525_milm_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
PXD029525_lihm_bpca <- pca(PXD029525_lihm_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
PXD029525_limm_bpca <- pca(PXD029525_limm_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
PXD029525_lilm_bpca <- pca(PXD029525_lilm_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
PXD029525_hihm_bpca <- as.data.frame(PXD029525_hihm_bpca@completeObs)
PXD029525_himm_bpca <- as.data.frame(PXD029525_himm_bpca@completeObs)
PXD029525_hilm_bpca <- as.data.frame(PXD029525_hilm_bpca@completeObs)
PXD029525_mihm_bpca <- as.data.frame(PXD029525_mihm_bpca@completeObs)
PXD029525_mimm_bpca <- as.data.frame(PXD029525_mimm_bpca@completeObs)
PXD029525_milm_bpca <- as.data.frame(PXD029525_milm_bpca@completeObs)
PXD029525_lihm_bpca <- as.data.frame(PXD029525_lihm_bpca@completeObs)
PXD029525_limm_bpca <- as.data.frame(PXD029525_limm_bpca@completeObs)
PXD029525_lilm_bpca <- as.data.frame(PXD029525_lilm_bpca@completeObs)

#svd
PXD029525_hihm_svd <- as.data.frame(impute.wrapper.SVD(PXD029525_hihm_na, K = 10))
PXD029525_himm_svd <- as.data.frame(impute.wrapper.SVD(PXD029525_himm_na, K = 10))
PXD029525_hilm_svd <- as.data.frame(impute.wrapper.SVD(PXD029525_hilm_na, K = 10))
PXD029525_mihm_svd <- as.data.frame(impute.wrapper.SVD(PXD029525_mihm_na, K = 10))
PXD029525_mimm_svd <- as.data.frame(impute.wrapper.SVD(PXD029525_mimm_na, K = 10))
PXD029525_milm_svd <- as.data.frame(impute.wrapper.SVD(PXD029525_milm_na, K = 10))
PXD029525_lihm_svd <- as.data.frame(impute.wrapper.SVD(PXD029525_lihm_na, K = 10))
PXD029525_limm_svd <- as.data.frame(impute.wrapper.SVD(PXD029525_limm_na, K = 10))
PXD029525_lilm_svd <- as.data.frame(impute.wrapper.SVD(PXD029525_lilm_na, K = 10))

#knn
PXD029525_hihm_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(PXD029525_hihm_na), K = 5))
PXD029525_himm_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(PXD029525_himm_na), K = 5))
PXD029525_hilm_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(PXD029525_hilm_na), K = 5))
PXD029525_mihm_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(PXD029525_mihm_na), K = 5))
PXD029525_mimm_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(PXD029525_mimm_na), K = 5))
PXD029525_milm_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(PXD029525_milm_na), K = 5))
PXD029525_lihm_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(PXD029525_lihm_na), K = 5))
PXD029525_limm_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(PXD029525_limm_na), K = 5))
PXD029525_lilm_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(PXD029525_lilm_na), K = 5))

#mle
PXD029525_hihm_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(PXD029525_hihm_na)))
PXD029525_himm_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(PXD029525_himm_na)))
PXD029525_hilm_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(PXD029525_hilm_na)))
PXD029525_mihm_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(PXD029525_mihm_na)))
PXD029525_mimm_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(PXD029525_mimm_na)))
PXD029525_milm_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(PXD029525_milm_na)))
PXD029525_lihm_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(PXD029525_lihm_na)))
PXD029525_limm_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(PXD029525_limm_na)))
PXD029525_lilm_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(PXD029525_lilm_na)))

#lls
PXD029525_hihm_lls <- llsImpute(PXD029525_hihm_na, k = 5, center = TRUE, completeObs = TRUE)
PXD029525_himm_lls <- llsImpute(PXD029525_himm_na, k = 5, center = TRUE, completeObs = TRUE)
PXD029525_hilm_lls <- llsImpute(PXD029525_hilm_na, k = 5, center = TRUE, completeObs = TRUE)
PXD029525_mihm_lls <- llsImpute(PXD029525_mihm_na, k = 5, center = TRUE, completeObs = TRUE)
PXD029525_mimm_lls <- llsImpute(PXD029525_mimm_na, k = 5, center = TRUE, completeObs = TRUE)
PXD029525_milm_lls <- llsImpute(PXD029525_milm_na, k = 5, center = TRUE, completeObs = TRUE)
PXD029525_lihm_lls <- llsImpute(PXD029525_lihm_na, k = 5, center = TRUE, completeObs = TRUE)
PXD029525_limm_lls <- llsImpute(PXD029525_limm_na, k = 5, center = TRUE, completeObs = TRUE)
PXD029525_lilm_lls <- llsImpute(PXD029525_lilm_na, k = 5, center = TRUE, completeObs = TRUE)
PXD029525_hihm_lls <- as.data.frame(PXD029525_hihm_lls@completeObs)
PXD029525_himm_lls <- as.data.frame(PXD029525_himm_lls@completeObs)
PXD029525_hilm_lls <- as.data.frame(PXD029525_hilm_lls@completeObs)
PXD029525_mihm_lls <- as.data.frame(PXD029525_mihm_lls@completeObs)
PXD029525_mimm_lls <- as.data.frame(PXD029525_mimm_lls@completeObs)
PXD029525_milm_lls <- as.data.frame(PXD029525_milm_lls@completeObs)
PXD029525_lihm_lls <- as.data.frame(PXD029525_lihm_lls@completeObs)
PXD029525_limm_lls <- as.data.frame(PXD029525_limm_lls@completeObs)
PXD029525_lilm_lls <- as.data.frame(PXD029525_lilm_lls@completeObs)

#rf
#register number of cores for parallel
registerDoParallel(cores = 18)
PXD029525_hihm_rf <- missForest(PXD029525_hihm_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
PXD029525_himm_rf <- missForest(PXD029525_himm_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
PXD029525_hilm_rf <- missForest(PXD029525_hilm_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
PXD029525_mihm_rf <- missForest(PXD029525_mihm_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
PXD029525_mimm_rf <- missForest(PXD029525_mimm_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
PXD029525_milm_rf <- missForest(PXD029525_milm_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
PXD029525_lihm_rf <- missForest(PXD029525_lihm_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
PXD029525_limm_rf <- missForest(PXD029525_limm_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
PXD029525_lilm_rf <- missForest(PXD029525_lilm_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
PXD029525_hihm_rf <- PXD029525_hihm_rf$ximp
PXD029525_himm_rf <- PXD029525_himm_rf$ximp
PXD029525_hilm_rf <- PXD029525_hilm_rf$ximp
PXD029525_mihm_rf <- PXD029525_mihm_rf$ximp
PXD029525_mimm_rf <- PXD029525_mimm_rf$ximp
PXD029525_milm_rf <- PXD029525_milm_rf$ximp
PXD029525_lihm_rf <- PXD029525_lihm_rf$ximp
PXD029525_limm_rf <- PXD029525_limm_rf$ximp
PXD029525_lilm_rf <- PXD029525_lilm_rf$ximp

#transpose data frame for pimms
write.csv(t(PXD029525_hihm_na), file = "PXD029525_hihm_na.csv")
write.csv(t(PXD029525_himm_na), file = "PXD029525_himm_na.csv")
write.csv(t(PXD029525_hilm_na), file = "PXD029525_hilm_na.csv")
write.csv(t(PXD029525_mihm_na), file = "PXD029525_mihm_na.csv")
write.csv(t(PXD029525_mimm_na), file = "PXD029525_mimm_na.csv")
write.csv(t(PXD029525_milm_na), file = "PXD029525_milm_na.csv")
write.csv(t(PXD029525_lihm_na), file = "PXD029525_lihm_na.csv")
write.csv(t(PXD029525_limm_na), file = "PXD029525_limm_na.csv")
write.csv(t(PXD029525_lilm_na), file = "PXD029525_lilm_na.csv")

#import pimms result
#cf
PXD029525_himm_cf <- read.csv("PXD029525_himm_cf.csv")
PXD029525_hilm_cf <- read.csv("PXD029525_hilm_cf.csv")
PXD029525_mihm_cf <- read.csv("PXD029525_mihm_cf.csv")
PXD029525_mimm_cf <- read.csv("PXD029525_mimm_cf.csv")
PXD029525_milm_cf <- read.csv("PXD029525_milm_cf.csv")
PXD029525_lihm_cf <- read.csv("PXD029525_lihm_cf.csv")
PXD029525_limm_cf <- read.csv("PXD029525_limm_cf.csv")
PXD029525_lilm_cf <- read.csv("PXD029525_lilm_cf.csv")
#vae
PXD029525_hihm_vae <- read.csv("PXD029525_hihm_vae.csv")
PXD029525_himm_vae <- read.csv("PXD029525_himm_vae.csv")
PXD029525_hilm_vae <- read.csv("PXD029525_hilm_vae.csv")
PXD029525_mihm_vae <- read.csv("PXD029525_mihm_vae.csv")
PXD029525_mimm_vae <- read.csv("PXD029525_mimm_vae.csv")
PXD029525_milm_vae <- read.csv("PXD029525_milm_vae.csv")
PXD029525_lihm_vae <- read.csv("PXD029525_lihm_vae.csv")
PXD029525_limm_vae <- read.csv("PXD029525_limm_vae.csv")
PXD029525_lilm_vae <- read.csv("PXD029525_lilm_vae.csv")
#dae
PXD029525_hihm_dae <- read.csv("PXD029525_hihm_dae.csv")
PXD029525_himm_dae <- read.csv("PXD029525_himm_dae.csv")
PXD029525_hilm_dae <- read.csv("PXD029525_hilm_dae.csv")
PXD029525_mihm_dae <- read.csv("PXD029525_mihm_dae.csv")
PXD029525_mimm_dae <- read.csv("PXD029525_mimm_dae.csv")
PXD029525_milm_dae <- read.csv("PXD029525_milm_dae.csv")
PXD029525_lihm_dae <- read.csv("PXD029525_lihm_dae.csv")
PXD029525_limm_dae <- read.csv("PXD029525_limm_dae.csv")
PXD029525_lilm_dae <- read.csv("PXD029525_lilm_dae.csv")

#set colname
colnames(PXD029525_hihm_dae) <- PXD029525_hihm_dae[1,]
colnames(PXD029525_himm_dae) <- PXD029525_himm_dae[1,]
colnames(PXD029525_hilm_dae) <- PXD029525_hilm_dae[1,]
colnames(PXD029525_mihm_dae) <- PXD029525_mihm_dae[1,]
colnames(PXD029525_mimm_dae) <- PXD029525_mimm_dae[1,]
colnames(PXD029525_milm_dae) <- PXD029525_milm_dae[1,]
colnames(PXD029525_lihm_dae) <- PXD029525_lihm_dae[1,]
colnames(PXD029525_limm_dae) <- PXD029525_limm_dae[1,]
colnames(PXD029525_lilm_dae) <- PXD029525_lilm_dae[1,]

colnames(PXD029525_hihm_vae) <- PXD029525_hihm_vae[1,]
colnames(PXD029525_himm_vae) <- PXD029525_himm_vae[1,]
colnames(PXD029525_hilm_vae) <- PXD029525_hilm_vae[1,]
colnames(PXD029525_mihm_vae) <- PXD029525_mihm_vae[1,]
colnames(PXD029525_mimm_vae) <- PXD029525_mimm_vae[1,]
colnames(PXD029525_milm_vae) <- PXD029525_milm_vae[1,]
colnames(PXD029525_lihm_vae) <- PXD029525_lihm_vae[1,]
colnames(PXD029525_limm_vae) <- PXD029525_limm_vae[1,]
colnames(PXD029525_lilm_vae) <- PXD029525_lilm_vae[1,]

#remove excessive rows
PXD029525_hihm_dae <- PXD029525_hihm_dae[-c(1,2),]
PXD029525_himm_dae <- PXD029525_himm_dae[-c(1,2),]
PXD029525_hilm_dae <- PXD029525_hilm_dae[-c(1,2),]
PXD029525_mihm_dae <- PXD029525_mihm_dae[-c(1,2),]
PXD029525_mimm_dae <- PXD029525_mimm_dae[-c(1,2),]
PXD029525_milm_dae <- PXD029525_milm_dae[-c(1,2),]
PXD029525_lihm_dae <- PXD029525_lihm_dae[-c(1,2),]
PXD029525_limm_dae <- PXD029525_limm_dae[-c(1,2),]
PXD029525_lilm_dae <- PXD029525_lilm_dae[-c(1,2),]

PXD029525_hihm_vae <- PXD029525_hihm_vae[-c(1,2),]
PXD029525_himm_vae <- PXD029525_himm_vae[-c(1,2),]
PXD029525_hilm_vae <- PXD029525_hilm_vae[-c(1,2),]
PXD029525_mihm_vae <- PXD029525_mihm_vae[-c(1,2),]
PXD029525_mimm_vae <- PXD029525_mimm_vae[-c(1,2),]
PXD029525_milm_vae <- PXD029525_milm_vae[-c(1,2),]
PXD029525_lihm_vae <- PXD029525_lihm_vae[-c(1,2),]
PXD029525_limm_vae <- PXD029525_limm_vae[-c(1,2),]
PXD029525_lilm_vae <- PXD029525_lilm_vae[-c(1,2),]

#turn chr into num
PXD029525_hihm_dae[,-1] <- sapply(PXD029525_hihm_dae[,-1], as.numeric)
PXD029525_himm_dae[,-1] <- sapply(PXD029525_himm_dae[,-1], as.numeric)
PXD029525_hilm_dae[,-1] <- sapply(PXD029525_hilm_dae[,-1], as.numeric)
PXD029525_mihm_dae[,-1] <- sapply(PXD029525_mihm_dae[,-1], as.numeric)
PXD029525_mimm_dae[,-1] <- sapply(PXD029525_mimm_dae[,-1], as.numeric)
PXD029525_milm_dae[,-1] <- sapply(PXD029525_milm_dae[,-1], as.numeric)
PXD029525_lihm_dae[,-1] <- sapply(PXD029525_lihm_dae[,-1], as.numeric)
PXD029525_limm_dae[,-1] <- sapply(PXD029525_limm_dae[,-1], as.numeric)
PXD029525_lilm_dae[,-1] <- sapply(PXD029525_lilm_dae[,-1], as.numeric)

PXD029525_hihm_vae[,-1] <- sapply(PXD029525_hihm_vae[,-1], as.numeric)
PXD029525_himm_vae[,-1] <- sapply(PXD029525_himm_vae[,-1], as.numeric)
PXD029525_hilm_vae[,-1] <- sapply(PXD029525_hilm_vae[,-1], as.numeric)
PXD029525_mihm_vae[,-1] <- sapply(PXD029525_mihm_vae[,-1], as.numeric)
PXD029525_mimm_vae[,-1] <- sapply(PXD029525_mimm_vae[,-1], as.numeric)
PXD029525_milm_vae[,-1] <- sapply(PXD029525_milm_vae[,-1], as.numeric)
PXD029525_lihm_vae[,-1] <- sapply(PXD029525_lihm_vae[,-1], as.numeric)
PXD029525_limm_vae[,-1] <- sapply(PXD029525_limm_vae[,-1], as.numeric)
PXD029525_lilm_vae[,-1] <- sapply(PXD029525_lilm_vae[,-1], as.numeric)

#set sample id as row name
#cf
rownames(PXD029525_himm_cf) <- PXD029525_himm_cf[,1]
rownames(PXD029525_hilm_cf) <- PXD029525_hilm_cf[,1]
rownames(PXD029525_mihm_cf) <- PXD029525_mihm_cf[,1]
rownames(PXD029525_mimm_cf) <- PXD029525_mimm_cf[,1]
rownames(PXD029525_milm_cf) <- PXD029525_milm_cf[,1]
rownames(PXD029525_lihm_cf) <- PXD029525_lihm_cf[,1]
rownames(PXD029525_limm_cf) <- PXD029525_limm_cf[,1]
rownames(PXD029525_lilm_cf) <- PXD029525_lilm_cf[,1]
#vae
rownames(PXD029525_hihm_vae) <- PXD029525_hihm_vae[,1]
rownames(PXD029525_himm_vae) <- PXD029525_himm_vae[,1]
rownames(PXD029525_hilm_vae) <- PXD029525_hilm_vae[,1]
rownames(PXD029525_mihm_vae) <- PXD029525_mihm_vae[,1]
rownames(PXD029525_mimm_vae) <- PXD029525_mimm_vae[,1]
rownames(PXD029525_milm_vae) <- PXD029525_milm_vae[,1]
rownames(PXD029525_lihm_vae) <- PXD029525_lihm_vae[,1]
rownames(PXD029525_limm_vae) <- PXD029525_limm_vae[,1]
rownames(PXD029525_lilm_vae) <- PXD029525_lilm_vae[,1]
#dae
rownames(PXD029525_hihm_dae) <- PXD029525_hihm_dae[,1]
rownames(PXD029525_himm_dae) <- PXD029525_himm_dae[,1]
rownames(PXD029525_hilm_dae) <- PXD029525_hilm_dae[,1]
rownames(PXD029525_mihm_dae) <- PXD029525_mihm_dae[,1]
rownames(PXD029525_mimm_dae) <- PXD029525_mimm_dae[,1]
rownames(PXD029525_milm_dae) <- PXD029525_milm_dae[,1]
rownames(PXD029525_lihm_dae) <- PXD029525_lihm_dae[,1]
rownames(PXD029525_limm_dae) <- PXD029525_limm_dae[,1]
rownames(PXD029525_lilm_dae) <- PXD029525_lilm_dae[,1]

#transpose
#cf
PXD029525_himm_cf <- as.data.frame(t(PXD029525_himm_cf[,-1]))
PXD029525_hilm_cf <- as.data.frame(t(PXD029525_hilm_cf[,-1]))
PXD029525_mihm_cf <- as.data.frame(t(PXD029525_mihm_cf[,-1]))
PXD029525_mimm_cf <- as.data.frame(t(PXD029525_mimm_cf[,-1]))
PXD029525_milm_cf <- as.data.frame(t(PXD029525_milm_cf[,-1]))
PXD029525_lihm_cf <- as.data.frame(t(PXD029525_lihm_cf[,-1]))
PXD029525_limm_cf <- as.data.frame(t(PXD029525_limm_cf[,-1]))
PXD029525_lilm_cf <- as.data.frame(t(PXD029525_lilm_cf[,-1]))
#vae
PXD029525_hihm_vae <- as.data.frame(t(PXD029525_hihm_vae[,-1]))
PXD029525_himm_vae <- as.data.frame(t(PXD029525_himm_vae[,-1]))
PXD029525_hilm_vae <- as.data.frame(t(PXD029525_hilm_vae[,-1]))
PXD029525_mihm_vae <- as.data.frame(t(PXD029525_mihm_vae[,-1]))
PXD029525_mimm_vae <- as.data.frame(t(PXD029525_mimm_vae[,-1]))
PXD029525_milm_vae <- as.data.frame(t(PXD029525_milm_vae[,-1]))
PXD029525_lihm_vae <- as.data.frame(t(PXD029525_lihm_vae[,-1]))
PXD029525_limm_vae <- as.data.frame(t(PXD029525_limm_vae[,-1]))
PXD029525_lilm_vae <- as.data.frame(t(PXD029525_lilm_vae[,-1]))
#dae
PXD029525_hihm_dae <- as.data.frame(t(PXD029525_hihm_dae[,-1]))
PXD029525_himm_dae <- as.data.frame(t(PXD029525_himm_dae[,-1]))
PXD029525_hilm_dae <- as.data.frame(t(PXD029525_hilm_dae[,-1]))
PXD029525_mihm_dae <- as.data.frame(t(PXD029525_mihm_dae[,-1]))
PXD029525_mimm_dae <- as.data.frame(t(PXD029525_mimm_dae[,-1]))
PXD029525_milm_dae <- as.data.frame(t(PXD029525_milm_dae[,-1]))
PXD029525_lihm_dae <- as.data.frame(t(PXD029525_lihm_dae[,-1]))
PXD029525_limm_dae <- as.data.frame(t(PXD029525_limm_dae[,-1]))
PXD029525_lilm_dae <- as.data.frame(t(PXD029525_lilm_dae[,-1]))

#evaluation by region
#hihm
comp_real <- PXD029525_hihm_real[which(rownames(PXD029525_hihm_real) %in% rownames(PXD029525_hihm_na)),] %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "real")
comp_na <- PXD029525_hihm_na %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "na")
comp_knn <- PXD029525_hihm_knn %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "knn")
comp_mle <- PXD029525_hihm_mle %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "mle")
comp_svd <- PXD029525_hihm_svd %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "svd")
comp_bpca <- PXD029525_hihm_bpca %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "bpca")
comp_lls <- PXD029525_hihm_lls %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "lls")
comp_rf <- PXD029525_hihm_rf %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "rf")
comp_vae <- PXD029525_hihm_vae %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "vae")
comp_dae <- PXD029525_hihm_dae %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "dae")

comp <- data.frame(sample = comp_real$sample,
                   real = comp_real$real,
                   na = comp_na$na,
                   knn = comp_knn$knn,
                   mle = comp_mle$mle,
                   svd = comp_svd$svd,
                   bpca = comp_bpca$bpca,
                   lls = comp_lls$lls,
                   rf = comp_rf$rf,
                   vae = comp_vae$vae,
                   dae = comp_dae$dae)
comp_long <- comp %>% pivot_longer(., 2:length(.), names_to = "method", values_to = "value")

#NRMSE
nrmse_sum <- data.frame(sample = unique(comp_real$sample),
                        nrmse_knn = NA,
                        nrmse_mle = NA,
                        nrmse_svd = NA, 
                        nrmse_bpca = NA, 
                        nrmse_lls = NA,
                        nrmse_rf = NA,
                        nrmse_cf = NA,
                        nrmse_vae = NA,
                        nrmse_dae = NA)

for (i in 1:length(nrmse_sum$sample)) {
  z <- PXD029525_hihm_real[which(rownames(PXD029525_hihm_real) %in% rownames(PXD029525_hihm_na) & rownames(PXD029525_hihm_real) %in% rownames(PXD029525_hihm_vae)),grep(nrmse_sum$sample[i], colnames(PXD029525_hihm_real))]
  #normalized by standard deviation of zero imputed matrix
  nrmse_sum$nrmse_knn[i] <- sqrt(mean((PXD029525_hihm_knn[which(rownames(PXD029525_hihm_knn) %in% rownames(PXD029525_hihm_vae)),grep(nrmse_sum$sample[i], colnames(PXD029525_hihm_knn))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_mle[i] <- sqrt(mean((PXD029525_hihm_mle[rownames(PXD029525_hihm_mle) %in% rownames(PXD029525_hihm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_hihm_mle))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_svd[i] <- sqrt(mean((PXD029525_hihm_svd[rownames(PXD029525_hihm_svd) %in% rownames(PXD029525_hihm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_hihm_svd))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_bpca[i] <- sqrt(mean((PXD029525_hihm_bpca[rownames(PXD029525_hihm_bpca) %in% rownames(PXD029525_hihm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_hihm_bpca))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_lls[i] <- sqrt(mean((PXD029525_hihm_lls[rownames(PXD029525_hihm_lls) %in% rownames(PXD029525_hihm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_hihm_lls))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_rf[i] <- sqrt(mean((PXD029525_hihm_rf[rownames(PXD029525_hihm_rf) %in% rownames(PXD029525_hihm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_hihm_rf))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_vae[i] <- sqrt(mean((PXD029525_hihm_vae[rownames(PXD029525_hihm_vae) %in% rownames(PXD029525_hihm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_hihm_vae))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_dae[i] <- sqrt(mean((PXD029525_hihm_dae[rownames(PXD029525_hihm_dae) %in% rownames(PXD029525_hihm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_hihm_dae))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
}

nrmse_long <- nrmse_sum %>% pivot_longer(cols = 2:length(.), values_to = "NRMSE", names_to = "method")
nrmse_long$method <- factor(nrmse_long$method, levels = c("nrmse_bpca", "nrmse_knn", "nrmse_lls", "nrmse_mle", 
                                                          "nrmse_svd", "nrmse_rf", "nrmse_cf", "nrmse_dae", "nrmse_vae"))

#barplot
#calculate mean and standard error
#drop mle
nrmse_stat <- data.frame(method = colnames(nrmse_sum)[-c(1)],
                         mean = colMeans(nrmse_sum[,-c(1)]),
                         se = colSds(as.matrix(nrmse_sum[,-c(1)]))/sqrt(nrow(nrmse_sum)))
nrmse_stat$method <- factor(nrmse_stat$method, levels = c("nrmse_rf", "nrmse_bpca", "nrmse_knn", "nrmse_lls",  "nrmse_mle",
                                                          "nrmse_svd", "nrmse_cf", "nrmse_dae", "nrmse_vae"))
ggplot(data = nrmse_stat, aes(x = method, y = mean, fill = method)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = mean - 1.96*se, ymax = mean + 1.96*se, width = 0.2)) + #95% confidence interval
  geom_text(aes(label = round(mean, digits = 2)), vjust = -5, size = 5) +
  scale_fill_discrete(labels = c("RF", "BPCA", "kNN", "LLS", "MLE", "SVD", "CF", "DAE", "VAE")) +
  ylim(0, 0.8) +
  labs(x = "HIHM",
       y = "NRMSE") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        panel.grid.major.x = element_blank(),
        panel.grid.major = element_line(linewidth = 1),
        legend.position = "none")
colMeans(nrmse_sum[,2:length(nrmse_sum)])
quants <- c(0, 0.25, 0.50, 0.75, 1)
apply(nrmse_sum[,2:length(nrmse_sum)], 2 , quantile , probs = quants , na.rm = TRUE )

#himm
comp_real <- PXD029525_himm_real[which(rownames(PXD029525_himm_real) %in% rownames(PXD029525_himm_na)),] %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "real")
comp_na <- PXD029525_himm_na %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "na")
comp_knn <- PXD029525_himm_knn %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "knn")
comp_mle <- PXD029525_himm_mle %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "mle")
comp_svd <- PXD029525_himm_svd %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "svd")
comp_bpca <- PXD029525_himm_bpca %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "bpca")
comp_lls <- PXD029525_himm_lls %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "lls")
comp_rf <- PXD029525_himm_rf %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "rf")
comp_cf <- PXD029525_himm_cf %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "cf")
comp_vae <- PXD029525_himm_vae %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "vae")
comp_dae <- PXD029525_himm_dae %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "dae")

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
                   dae = comp_dae$dae)
comp_long <- comp %>% pivot_longer(., 2:length(.), names_to = "method", values_to = "value")

#NRMSE
nrmse_sum <- data.frame(sample = unique(comp_real$sample),
                        nrmse_knn = NA,
                        nrmse_mle = NA,
                        nrmse_svd = NA, 
                        nrmse_bpca = NA, 
                        nrmse_lls = NA,
                        nrmse_rf = NA,
                        nrmse_cf = NA,
                        nrmse_vae = NA,
                        nrmse_dae = NA)

for (i in 1:length(nrmse_sum$sample)) {
  z <- PXD029525_himm_real[which(rownames(PXD029525_himm_real) %in% rownames(PXD029525_himm_na) & rownames(PXD029525_himm_real) %in% rownames(PXD029525_himm_vae)),grep(nrmse_sum$sample[i], colnames(PXD029525_himm_real))]
  #normalized by standard deviation of zero imputed matrix
  nrmse_sum$nrmse_knn[i] <- sqrt(mean((PXD029525_himm_knn[which(rownames(PXD029525_himm_knn) %in% rownames(PXD029525_himm_vae)),grep(nrmse_sum$sample[i], colnames(PXD029525_himm_knn))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_mle[i] <- sqrt(mean((PXD029525_himm_mle[rownames(PXD029525_himm_mle) %in% rownames(PXD029525_himm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_himm_mle))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_svd[i] <- sqrt(mean((PXD029525_himm_svd[rownames(PXD029525_himm_svd) %in% rownames(PXD029525_himm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_himm_svd))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_bpca[i] <- sqrt(mean((PXD029525_himm_bpca[rownames(PXD029525_himm_bpca) %in% rownames(PXD029525_himm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_himm_bpca))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_lls[i] <- sqrt(mean((PXD029525_himm_lls[rownames(PXD029525_himm_lls) %in% rownames(PXD029525_himm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_himm_lls))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_rf[i] <- sqrt(mean((PXD029525_himm_rf[rownames(PXD029525_himm_rf) %in% rownames(PXD029525_himm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_himm_rf))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_cf[i] <- sqrt(mean((PXD029525_himm_cf[rownames(PXD029525_himm_cf) %in% rownames(PXD029525_himm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_himm_cf))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_vae[i] <- sqrt(mean((PXD029525_himm_vae[rownames(PXD029525_himm_vae) %in% rownames(PXD029525_himm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_himm_vae))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_dae[i] <- sqrt(mean((PXD029525_himm_dae[rownames(PXD029525_himm_dae) %in% rownames(PXD029525_himm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_himm_dae))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
}

nrmse_long <- nrmse_sum %>% pivot_longer(cols = 2:length(.), values_to = "NRMSE", names_to = "method")
nrmse_long$method <- factor(nrmse_long$method, levels = c("nrmse_bpca", "nrmse_knn", "nrmse_lls", "nrmse_mle", 
                                                          "nrmse_svd", "nrmse_rf", "nrmse_cf", "nrmse_dae", "nrmse_vae"))

#barplot
#calculate mean and standard error
nrmse_stat <- data.frame(method = colnames(nrmse_sum)[-c(1)],
                         mean = colMeans(nrmse_sum[,-c(1)]),
                         se = colSds(as.matrix(nrmse_sum[,-c(1)]))/sqrt(nrow(nrmse_sum)))
nrmse_stat$method <- factor(nrmse_stat$method, levels = c("nrmse_rf", "nrmse_bpca", "nrmse_knn", "nrmse_lls",  "nrmse_mle",
                                                          "nrmse_svd", "nrmse_cf", "nrmse_dae", "nrmse_vae"))
ggplot(data = nrmse_stat, aes(x = method, y = mean, fill = method)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = mean - 1.96*se, ymax = mean + 1.96*se, width = 0.2)) + #95% confidence interval
  geom_text(aes(label = round(mean, digits = 2)), vjust = -5, size = 5) +
  scale_fill_discrete(labels = c("RF", "BPCA", "kNN", "LLS", "MLE", "SVD", "CF", "DAE", "VAE")) +
  ylim(0, 0.8) +
  labs(x = "HIMM",
       y = "NRMSE") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        panel.grid.major.x = element_blank(),
        panel.grid.major = element_line(linewidth = 1),
        legend.position = "none")
colMeans(nrmse_sum[,2:length(nrmse_sum)])
quants <- c(0, 0.25, 0.50, 0.75, 1)
apply(nrmse_sum[,2:length(nrmse_sum)], 2 , quantile , probs = quants , na.rm = TRUE )

#hilm
comp_real <- PXD029525_hilm_real[which(rownames(PXD029525_hilm_real) %in% rownames(PXD029525_hilm_na)),] %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "real")
comp_na <- PXD029525_hilm_na %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "na")
comp_knn <- PXD029525_hilm_knn %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "knn")
comp_mle <- PXD029525_hilm_mle %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "mle")
comp_svd <- PXD029525_hilm_svd %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "svd")
comp_bpca <- PXD029525_hilm_bpca %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "bpca")
comp_lls <- PXD029525_hilm_lls %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "lls")
comp_rf <- PXD029525_hilm_rf %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "rf")
comp_cf <- PXD029525_hilm_cf %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "cf")
comp_vae <- PXD029525_hilm_vae %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "vae")
comp_dae <- PXD029525_hilm_dae %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "dae")

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
                   dae = comp_dae$dae)
comp_long <- comp %>% pivot_longer(., 2:length(.), names_to = "method", values_to = "value")

#NRMSE
nrmse_sum <- data.frame(sample = unique(comp_real$sample),
                        nrmse_knn = NA,
                        nrmse_mle = NA,
                        nrmse_svd = NA, 
                        nrmse_bpca = NA, 
                        nrmse_lls = NA,
                        nrmse_rf = NA,
                        nrmse_cf = NA,
                        nrmse_vae = NA,
                        nrmse_dae = NA)

for (i in 1:length(nrmse_sum$sample)) {
  z <- PXD029525_hilm_real[which(rownames(PXD029525_hilm_real) %in% rownames(PXD029525_hilm_na) & rownames(PXD029525_hilm_real) %in% rownames(PXD029525_hilm_vae)),grep(nrmse_sum$sample[i], colnames(PXD029525_hilm_real))]
  #normalized by standard deviation of zero imputed matrix
  nrmse_sum$nrmse_knn[i] <- sqrt(mean((PXD029525_hilm_knn[which(rownames(PXD029525_hilm_knn) %in% rownames(PXD029525_hilm_vae)),grep(nrmse_sum$sample[i], colnames(PXD029525_hilm_knn))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_mle[i] <- sqrt(mean((PXD029525_hilm_mle[rownames(PXD029525_hilm_mle) %in% rownames(PXD029525_hilm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_hilm_mle))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_svd[i] <- sqrt(mean((PXD029525_hilm_svd[rownames(PXD029525_hilm_svd) %in% rownames(PXD029525_hilm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_hilm_svd))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_bpca[i] <- sqrt(mean((PXD029525_hilm_bpca[rownames(PXD029525_hilm_bpca) %in% rownames(PXD029525_hilm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_hilm_bpca))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_lls[i] <- sqrt(mean((PXD029525_hilm_lls[rownames(PXD029525_hilm_lls) %in% rownames(PXD029525_hilm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_hilm_lls))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_rf[i] <- sqrt(mean((PXD029525_hilm_rf[rownames(PXD029525_hilm_rf) %in% rownames(PXD029525_hilm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_hilm_rf))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_cf[i] <- sqrt(mean((PXD029525_hilm_cf[rownames(PXD029525_hilm_cf) %in% rownames(PXD029525_hilm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_hilm_cf))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_vae[i] <- sqrt(mean((PXD029525_hilm_vae[rownames(PXD029525_hilm_vae) %in% rownames(PXD029525_hilm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_hilm_vae))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_dae[i] <- sqrt(mean((PXD029525_hilm_dae[rownames(PXD029525_hilm_dae) %in% rownames(PXD029525_hilm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_hilm_dae))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
}

nrmse_long <- nrmse_sum %>% pivot_longer(cols = 2:length(.), values_to = "NRMSE", names_to = "method")
nrmse_long$method <- factor(nrmse_long$method, levels = c("nrmse_bpca", "nrmse_knn", "nrmse_lls", "nrmse_mle", 
                                                          "nrmse_svd", "nrmse_rf", "nrmse_cf", "nrmse_dae", "nrmse_vae"))

#barplot
#calculate mean and standard error
nrmse_stat <- data.frame(method = colnames(nrmse_sum)[-c(1)],
                         mean = colMeans(nrmse_sum[,-c(1)]),
                         se = colSds(as.matrix(nrmse_sum[,-c(1)]))/sqrt(nrow(nrmse_sum)))
nrmse_stat$method <- factor(nrmse_stat$method, levels = c("nrmse_rf", "nrmse_bpca", "nrmse_knn", "nrmse_lls",  "nrmse_mle",
                                                          "nrmse_svd", "nrmse_cf", "nrmse_dae", "nrmse_vae"))
ggplot(data = nrmse_stat, aes(x = method, y = mean, fill = method)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = mean - 1.96*se, ymax = mean + 1.96*se, width = 0.2)) + #95% confidence interval
  geom_text(aes(label = round(mean, digits = 2)), vjust = -5, size = 5) +
  scale_fill_discrete(labels = c("RF", "BPCA", "kNN", "LLS", "MLE", "SVD", "CF", "DAE", "VAE")) +
  ylim(0, 0.8) +
  labs(x = "HILM",
       y = "NRMSE") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        panel.grid.major.x = element_blank(),
        panel.grid.major = element_line(linewidth = 1),
        legend.position = "none")
colMeans(nrmse_sum[,2:length(nrmse_sum)])
quants <- c(0, 0.25, 0.50, 0.75, 1)
apply(nrmse_sum[,2:length(nrmse_sum)], 2 , quantile , probs = quants , na.rm = TRUE )

#mihm
comp_real <- PXD029525_mihm_real[which(rownames(PXD029525_mihm_real) %in% rownames(PXD029525_mihm_na)),] %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "real")
comp_na <- PXD029525_mihm_na %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "na")
comp_knn <- PXD029525_mihm_knn %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "knn")
comp_mle <- PXD029525_mihm_mle %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "mle")
comp_svd <- PXD029525_mihm_svd %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "svd")
comp_bpca <- PXD029525_mihm_bpca %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "bpca")
comp_lls <- PXD029525_mihm_lls %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "lls")
comp_rf <- PXD029525_mihm_rf %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "rf")
comp_cf <- PXD029525_mihm_cf %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "cf")
comp_vae <- PXD029525_mihm_vae %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "vae")
comp_dae <- PXD029525_mihm_dae %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "dae")

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
                   dae = comp_dae$dae)
comp_long <- comp %>% pivot_longer(., 2:length(.), names_to = "method", values_to = "value")

#NRMSE
nrmse_sum <- data.frame(sample = unique(comp_real$sample),
                        nrmse_knn = NA,
                        nrmse_mle = NA,
                        nrmse_svd = NA, 
                        nrmse_bpca = NA, 
                        nrmse_lls = NA,
                        nrmse_rf = NA,
                        nrmse_cf = NA,
                        nrmse_vae = NA,
                        nrmse_dae = NA)

for (i in 1:length(nrmse_sum$sample)) {
  z <- PXD029525_mihm_real[which(rownames(PXD029525_mihm_real) %in% rownames(PXD029525_mihm_na) & rownames(PXD029525_mihm_real) %in% rownames(PXD029525_mihm_vae)),grep(nrmse_sum$sample[i], colnames(PXD029525_mihm_real))]
  #normalized by standard deviation of zero imputed matrix
  nrmse_sum$nrmse_knn[i] <- sqrt(mean((PXD029525_mihm_knn[which(rownames(PXD029525_mihm_knn) %in% rownames(PXD029525_mihm_vae)),grep(nrmse_sum$sample[i], colnames(PXD029525_mihm_knn))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_mle[i] <- sqrt(mean((PXD029525_mihm_mle[rownames(PXD029525_mihm_mle) %in% rownames(PXD029525_mihm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_mihm_mle))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_svd[i] <- sqrt(mean((PXD029525_mihm_svd[rownames(PXD029525_mihm_svd) %in% rownames(PXD029525_mihm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_mihm_svd))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_bpca[i] <- sqrt(mean((PXD029525_mihm_bpca[rownames(PXD029525_mihm_bpca) %in% rownames(PXD029525_mihm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_mihm_bpca))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_lls[i] <- sqrt(mean((PXD029525_mihm_lls[rownames(PXD029525_mihm_lls) %in% rownames(PXD029525_mihm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_mihm_lls))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_rf[i] <- sqrt(mean((PXD029525_mihm_rf[rownames(PXD029525_mihm_rf) %in% rownames(PXD029525_mihm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_mihm_rf))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_cf[i] <- sqrt(mean((PXD029525_mihm_cf[rownames(PXD029525_mihm_cf) %in% rownames(PXD029525_mihm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_mihm_cf))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_vae[i] <- sqrt(mean((PXD029525_mihm_vae[rownames(PXD029525_mihm_vae) %in% rownames(PXD029525_mihm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_mihm_vae))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_dae[i] <- sqrt(mean((PXD029525_mihm_dae[rownames(PXD029525_mihm_dae) %in% rownames(PXD029525_mihm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_mihm_dae))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
}

nrmse_long <- nrmse_sum %>% pivot_longer(cols = 2:length(.), values_to = "NRMSE", names_to = "method")
nrmse_long$method <- factor(nrmse_long$method, levels = c("nrmse_bpca", "nrmse_knn", "nrmse_lls", "nrmse_mle", 
                                                          "nrmse_svd", "nrmse_rf", "nrmse_cf", "nrmse_dae", "nrmse_vae"))

#barplot
#calculate mean and standard error
nrmse_stat <- data.frame(method = colnames(nrmse_sum)[-c(1)],
                         mean = colMeans(nrmse_sum[,-c(1)]),
                         se = colSds(as.matrix(nrmse_sum[,-c(1)]))/sqrt(nrow(nrmse_sum)))
nrmse_stat$method <- factor(nrmse_stat$method, levels = c("nrmse_rf", "nrmse_bpca", "nrmse_knn", "nrmse_lls",  "nrmse_mle",
                                                          "nrmse_svd", "nrmse_cf", "nrmse_dae", "nrmse_vae"))
ggplot(data = nrmse_stat, aes(x = method, y = mean, fill = method)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = mean - 1.96*se, ymax = mean + 1.96*se, width = 0.2)) + #95% confidence interval
  geom_text(aes(label = round(mean, digits = 2)), vjust = -5, size = 5) +
  scale_fill_discrete(labels = c("RF", "BPCA", "kNN", "LLS", "MLE", "SVD", "CF", "DAE", "VAE")) +
  ylim(0, 1.5) +
  labs(x = "MIHM",
       y = "NRMSE") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        panel.grid.major.x = element_blank(),
        panel.grid.major = element_line(linewidth = 1),
        legend.position = "none")
colMeans(nrmse_sum[,2:length(nrmse_sum)])
quants <- c(0, 0.25, 0.50, 0.75, 1)
apply(nrmse_sum[,2:length(nrmse_sum)], 2 , quantile , probs = quants , na.rm = TRUE )

#mimm
comp_real <- PXD029525_mimm_real[which(rownames(PXD029525_mimm_real) %in% rownames(PXD029525_mimm_na)),] %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "real")
comp_na <- PXD029525_mimm_na %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "na")
comp_knn <- PXD029525_mimm_knn %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "knn")
comp_mle <- PXD029525_mimm_mle %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "mle")
comp_svd <- PXD029525_mimm_svd %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "svd")
comp_bpca <- PXD029525_mimm_bpca %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "bpca")
comp_lls <- PXD029525_mimm_lls %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "lls")
comp_rf <- PXD029525_mimm_rf %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "rf")
comp_cf <- PXD029525_mimm_cf %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "cf")
comp_vae <- PXD029525_mimm_vae %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "vae")
comp_dae <- PXD029525_mimm_dae %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "dae")

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
                   dae = comp_dae$dae)
comp_long <- comp %>% pivot_longer(., 2:length(.), names_to = "method", values_to = "value")

#NRMSE
nrmse_sum <- data.frame(sample = unique(comp_real$sample),
                        nrmse_knn = NA,
                        nrmse_mle = NA,
                        nrmse_svd = NA, 
                        nrmse_bpca = NA, 
                        nrmse_lls = NA,
                        nrmse_rf = NA,
                        nrmse_cf = NA,
                        nrmse_vae = NA,
                        nrmse_dae = NA)

for (i in 1:length(nrmse_sum$sample)) {
  z <- PXD029525_mimm_real[which(rownames(PXD029525_mimm_real) %in% rownames(PXD029525_mimm_na) & rownames(PXD029525_mimm_real) %in% rownames(PXD029525_mimm_vae)),grep(nrmse_sum$sample[i], colnames(PXD029525_mimm_real))]
  #normalized by standard deviation of zero imputed matrix
  nrmse_sum$nrmse_knn[i] <- sqrt(mean((PXD029525_mimm_knn[which(rownames(PXD029525_mimm_knn) %in% rownames(PXD029525_mimm_vae)),grep(nrmse_sum$sample[i], colnames(PXD029525_mimm_knn))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_mle[i] <- sqrt(mean((PXD029525_mimm_mle[rownames(PXD029525_mimm_mle) %in% rownames(PXD029525_mimm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_mimm_mle))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_svd[i] <- sqrt(mean((PXD029525_mimm_svd[rownames(PXD029525_mimm_svd) %in% rownames(PXD029525_mimm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_mimm_svd))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_bpca[i] <- sqrt(mean((PXD029525_mimm_bpca[rownames(PXD029525_mimm_bpca) %in% rownames(PXD029525_mimm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_mimm_bpca))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_lls[i] <- sqrt(mean((PXD029525_mimm_lls[rownames(PXD029525_mimm_lls) %in% rownames(PXD029525_mimm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_mimm_lls))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_rf[i] <- sqrt(mean((PXD029525_mimm_rf[rownames(PXD029525_mimm_rf) %in% rownames(PXD029525_mimm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_mimm_rf))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_cf[i] <- sqrt(mean((PXD029525_mimm_cf[rownames(PXD029525_mimm_cf) %in% rownames(PXD029525_mimm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_mimm_cf))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_vae[i] <- sqrt(mean((PXD029525_mimm_vae[rownames(PXD029525_mimm_vae) %in% rownames(PXD029525_mimm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_mimm_vae))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_dae[i] <- sqrt(mean((PXD029525_mimm_dae[rownames(PXD029525_mimm_dae) %in% rownames(PXD029525_mimm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_mimm_dae))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
}

nrmse_long <- nrmse_sum %>% pivot_longer(cols = 2:length(.), values_to = "NRMSE", names_to = "method")
nrmse_long$method <- factor(nrmse_long$method, levels = c("nrmse_bpca", "nrmse_knn", "nrmse_lls", "nrmse_mle", 
                                                          "nrmse_svd", "nrmse_rf", "nrmse_cf", "nrmse_dae", "nrmse_vae"))

#barplot
#calculate mean and standard error
nrmse_stat <- data.frame(method = colnames(nrmse_sum)[-c(1)],
                         mean = colMeans(nrmse_sum[,-c(1)]),
                         se = colSds(as.matrix(nrmse_sum[,-c(1)]))/sqrt(nrow(nrmse_sum)))
nrmse_stat$method <- factor(nrmse_stat$method, levels = c("nrmse_rf", "nrmse_bpca", "nrmse_knn", "nrmse_lls",  "nrmse_mle",
                                                          "nrmse_svd", "nrmse_cf", "nrmse_dae", "nrmse_vae"))
ggplot(data = nrmse_stat, aes(x = method, y = mean, fill = method)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = mean - 1.96*se, ymax = mean + 1.96*se, width = 0.2)) + #95% confidence interval
  geom_text(aes(label = round(mean, digits = 2)), vjust = -5, size = 5) +
  scale_fill_discrete(labels = c("RF", "BPCA", "kNN", "LLS", "MLE", "SVD", "CF", "DAE", "VAE")) +
  ylim(0, 1) +
  labs(x = "MIMM",
       y = "NRMSE") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        panel.grid.major.x = element_blank(),
        panel.grid.major = element_line(linewidth = 1),
        legend.position = "none")
colMeans(nrmse_sum[,2:length(nrmse_sum)])
quants <- c(0, 0.25, 0.50, 0.75, 1)
apply(nrmse_sum[,2:length(nrmse_sum)], 2 , quantile , probs = quants , na.rm = TRUE )

#hilm
comp_real <- PXD029525_hilm_real[which(rownames(PXD029525_hilm_real) %in% rownames(PXD029525_hilm_na)),] %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "real")
comp_na <- PXD029525_hilm_na %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "na")
comp_knn <- PXD029525_hilm_knn %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "knn")
comp_mle <- PXD029525_hilm_mle %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "mle")
comp_svd <- PXD029525_hilm_svd %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "svd")
comp_bpca <- PXD029525_hilm_bpca %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "bpca")
comp_lls <- PXD029525_hilm_lls %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "lls")
comp_rf <- PXD029525_hilm_rf %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "rf")
comp_cf <- PXD029525_hilm_cf %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "cf")
comp_vae <- PXD029525_hilm_vae %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "vae")
comp_dae <- PXD029525_hilm_dae %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "dae")

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
                   dae = comp_dae$dae)
comp_long <- comp %>% pivot_longer(., 2:length(.), names_to = "method", values_to = "value")

#NRMSE
nrmse_sum <- data.frame(sample = unique(comp_real$sample),
                        nrmse_knn = NA,
                        nrmse_mle = NA,
                        nrmse_svd = NA, 
                        nrmse_bpca = NA, 
                        nrmse_lls = NA,
                        nrmse_rf = NA,
                        nrmse_cf = NA,
                        nrmse_vae = NA,
                        nrmse_dae = NA)

for (i in 1:length(nrmse_sum$sample)) {
  z <- PXD029525_hilm_real[which(rownames(PXD029525_hilm_real) %in% rownames(PXD029525_hilm_na) & rownames(PXD029525_hilm_real) %in% rownames(PXD029525_hilm_vae)),grep(nrmse_sum$sample[i], colnames(PXD029525_hilm_real))]
  #normalized by standard deviation of zero imputed matrix
  nrmse_sum$nrmse_knn[i] <- sqrt(mean((PXD029525_hilm_knn[which(rownames(PXD029525_hilm_knn) %in% rownames(PXD029525_hilm_vae)),grep(nrmse_sum$sample[i], colnames(PXD029525_hilm_knn))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_mle[i] <- sqrt(mean((PXD029525_hilm_mle[rownames(PXD029525_hilm_mle) %in% rownames(PXD029525_hilm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_hilm_mle))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_svd[i] <- sqrt(mean((PXD029525_hilm_svd[rownames(PXD029525_hilm_svd) %in% rownames(PXD029525_hilm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_hilm_svd))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_bpca[i] <- sqrt(mean((PXD029525_hilm_bpca[rownames(PXD029525_hilm_bpca) %in% rownames(PXD029525_hilm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_hilm_bpca))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_lls[i] <- sqrt(mean((PXD029525_hilm_lls[rownames(PXD029525_hilm_lls) %in% rownames(PXD029525_hilm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_hilm_lls))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_rf[i] <- sqrt(mean((PXD029525_hilm_rf[rownames(PXD029525_hilm_rf) %in% rownames(PXD029525_hilm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_hilm_rf))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_cf[i] <- sqrt(mean((PXD029525_hilm_cf[rownames(PXD029525_hilm_cf) %in% rownames(PXD029525_hilm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_hilm_cf))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_vae[i] <- sqrt(mean((PXD029525_hilm_vae[rownames(PXD029525_hilm_vae) %in% rownames(PXD029525_hilm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_hilm_vae))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_dae[i] <- sqrt(mean((PXD029525_hilm_dae[rownames(PXD029525_hilm_dae) %in% rownames(PXD029525_hilm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_hilm_dae))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
}

nrmse_long <- nrmse_sum %>% pivot_longer(cols = 2:length(.), values_to = "NRMSE", names_to = "method")
nrmse_long$method <- factor(nrmse_long$method, levels = c("nrmse_bpca", "nrmse_knn", "nrmse_lls", "nrmse_mle", 
                                                          "nrmse_svd", "nrmse_rf", "nrmse_cf", "nrmse_dae", "nrmse_vae"))

#barplot
#calculate mean and standard error
nrmse_stat <- data.frame(method = colnames(nrmse_sum)[-c(1)],
                         mean = colMeans(nrmse_sum[,-c(1)]),
                         se = colSds(as.matrix(nrmse_sum[,-c(1)]))/sqrt(nrow(nrmse_sum)))
nrmse_stat$method <- factor(nrmse_stat$method, levels = c("nrmse_rf", "nrmse_bpca", "nrmse_knn", "nrmse_lls",  "nrmse_mle",
                                                          "nrmse_svd", "nrmse_cf", "nrmse_dae", "nrmse_vae"))
ggplot(data = nrmse_stat, aes(x = method, y = mean, fill = method)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = mean - 1.96*se, ymax = mean + 1.96*se, width = 0.2)) + #95% confidence interval
  geom_text(aes(label = round(mean, digits = 2)), vjust = -5, size = 5) +
  scale_fill_discrete(labels = c("RF", "BPCA", "kNN", "LLS", "MLE", "SVD", "CF", "DAE", "VAE")) +
  ylim(0, 0.8) +
  labs(x = "MILM",
       y = "NRMSE") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        panel.grid.major.x = element_blank(),
        panel.grid.major = element_line(linewidth = 1),
        legend.position = "none")
colMeans(nrmse_sum[,2:length(nrmse_sum)])
quants <- c(0, 0.25, 0.50, 0.75, 1)
apply(nrmse_sum[,2:length(nrmse_sum)], 2 , quantile , probs = quants , na.rm = TRUE )

#lihm
comp_real <- PXD029525_lihm_real[which(rownames(PXD029525_lihm_real) %in% rownames(PXD029525_lihm_na)),] %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "real")
comp_na <- PXD029525_lihm_na %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "na")
comp_mle <- PXD029525_lihm_mle %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "mle")
comp_svd <- PXD029525_lihm_svd %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "svd")
comp_bpca <- PXD029525_lihm_bpca %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "bpca")
comp_lls <- PXD029525_lihm_lls %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "lls")
comp_rf <- PXD029525_lihm_rf %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "rf")
comp_cf <- PXD029525_lihm_cf %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "cf")
comp_vae <- PXD029525_lihm_vae %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "vae")
comp_dae <- PXD029525_lihm_dae %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "dae")

comp <- data.frame(sample = comp_real$sample,
                   real = comp_real$real,
                   na = comp_na$na,
                   mle = comp_mle$mle,
                   svd = comp_svd$svd,
                   bpca = comp_bpca$bpca,
                   lls = comp_lls$lls,
                   rf = comp_rf$rf,
                   cf = comp_cf$cf,
                   vae = comp_vae$vae,
                   dae = comp_dae$dae)
comp_long <- comp %>% pivot_longer(., 2:length(.), names_to = "method", values_to = "value")

#NRMSE
nrmse_sum <- data.frame(sample = unique(comp_real$sample),
                        nrmse_knn = NA,
                        nrmse_mle = NA,
                        nrmse_svd = NA, 
                        nrmse_bpca = NA, 
                        nrmse_lls = NA,
                        nrmse_rf = NA,
                        nrmse_cf = NA,
                        nrmse_vae = NA,
                        nrmse_dae = NA)

for (i in 1:length(nrmse_sum$sample)) {
  z <- PXD029525_lihm_real[which(rownames(PXD029525_lihm_real) %in% rownames(PXD029525_lihm_na) & rownames(PXD029525_lihm_real) %in% rownames(PXD029525_lihm_vae)),grep(nrmse_sum$sample[i], colnames(PXD029525_lihm_real))]
  #normalized by standard deviation of zero imputed matrix
  nrmse_sum$nrmse_mle[i] <- sqrt(mean((PXD029525_lihm_mle[rownames(PXD029525_lihm_mle) %in% rownames(PXD029525_lihm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_lihm_mle))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_svd[i] <- sqrt(mean((PXD029525_lihm_svd[rownames(PXD029525_lihm_svd) %in% rownames(PXD029525_lihm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_lihm_svd))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_bpca[i] <- sqrt(mean((PXD029525_lihm_bpca[rownames(PXD029525_lihm_bpca) %in% rownames(PXD029525_lihm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_lihm_bpca))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_lls[i] <- sqrt(mean((PXD029525_lihm_lls[rownames(PXD029525_lihm_lls) %in% rownames(PXD029525_lihm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_lihm_lls))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_rf[i] <- sqrt(mean((PXD029525_lihm_rf[rownames(PXD029525_lihm_rf) %in% rownames(PXD029525_lihm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_lihm_rf))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_cf[i] <- sqrt(mean((PXD029525_lihm_cf[rownames(PXD029525_lihm_cf) %in% rownames(PXD029525_lihm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_lihm_cf))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_vae[i] <- sqrt(mean((PXD029525_lihm_vae[rownames(PXD029525_lihm_vae) %in% rownames(PXD029525_lihm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_lihm_vae))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_dae[i] <- sqrt(mean((PXD029525_lihm_dae[rownames(PXD029525_lihm_dae) %in% rownames(PXD029525_lihm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_lihm_dae))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
}

nrmse_long <- nrmse_sum %>% pivot_longer(cols = 2:length(.), values_to = "NRMSE", names_to = "method")
nrmse_long$method <- factor(nrmse_long$method, levels = c("nrmse_bpca", "nrmse_knn", "nrmse_lls", "nrmse_mle", 
                                                          "nrmse_svd", "nrmse_rf", "nrmse_cf", "nrmse_dae", "nrmse_vae"))

#barplot
#calculate mean and standard error
nrmse_stat <- data.frame(method = colnames(nrmse_sum)[-c(1)],
                         mean = colMeans(nrmse_sum[,-c(1)]),
                         se = colSds(as.matrix(nrmse_sum[,-c(1)]))/sqrt(nrow(nrmse_sum)))
nrmse_stat$method <- factor(nrmse_stat$method, levels = c("nrmse_rf", "nrmse_bpca", "nrmse_knn", "nrmse_lls",  "nrmse_mle",
                                                          "nrmse_svd", "nrmse_cf", "nrmse_dae", "nrmse_vae"))
ggplot(data = nrmse_stat, aes(x = method, y = mean, fill = method)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = mean - 1.96*se, ymax = mean + 1.96*se, width = 0.2)) + #95% confidence interval
  geom_text(aes(label = round(mean, digits = 2)), vjust = -5, size = 5) +
  scale_fill_discrete(labels = c("RF", "BPCA", "kNN", "LLS", "MLE", "SVD", "CF", "DAE", "VAE")) +
  ylim(0, 0.8) +
  labs(x = "LIHM",
       y = "NRMSE") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        panel.grid.major.x = element_blank(),
        panel.grid.major = element_line(linewidth = 1),
        legend.position = "none")
colMeans(nrmse_sum[,2:length(nrmse_sum)])
quants <- c(0, 0.25, 0.50, 0.75, 1)
apply(nrmse_sum[,2:length(nrmse_sum)], 2 , quantile , probs = quants , na.rm = TRUE )

#limm
comp_real <- PXD029525_limm_real[which(rownames(PXD029525_limm_real) %in% rownames(PXD029525_limm_na)),] %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "real")
comp_na <- PXD029525_limm_na %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "na")
comp_mle <- PXD029525_limm_mle %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "mle")
comp_svd <- PXD029525_limm_svd %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "svd")
comp_bpca <- PXD029525_limm_bpca %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "bpca")
comp_lls <- PXD029525_limm_lls %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "lls")
comp_rf <- PXD029525_limm_rf %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "rf")
comp_cf <- PXD029525_limm_cf %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "cf")
comp_vae <- PXD029525_limm_vae %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "vae")
comp_dae <- PXD029525_limm_dae %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "dae")

comp <- data.frame(sample = comp_real$sample,
                   real = comp_real$real,
                   na = comp_na$na,
                   mle = comp_mle$mle,
                   svd = comp_svd$svd,
                   bpca = comp_bpca$bpca,
                   lls = comp_lls$lls,
                   rf = comp_rf$rf,
                   cf = comp_cf$cf,
                   vae = comp_vae$vae,
                   dae = comp_dae$dae)
comp_long <- comp %>% pivot_longer(., 2:length(.), names_to = "method", values_to = "value")

#NRMSE
nrmse_sum <- data.frame(sample = unique(comp_real$sample),
                        nrmse_knn = NA,
                        nrmse_mle = NA,
                        nrmse_svd = NA, 
                        nrmse_bpca = NA, 
                        nrmse_lls = NA,
                        nrmse_rf = NA,
                        nrmse_cf = NA,
                        nrmse_vae = NA,
                        nrmse_dae = NA)

for (i in 1:length(nrmse_sum$sample)) {
  z <- PXD029525_limm_real[which(rownames(PXD029525_limm_real) %in% rownames(PXD029525_limm_na) & rownames(PXD029525_limm_real) %in% rownames(PXD029525_limm_vae)),grep(nrmse_sum$sample[i], colnames(PXD029525_limm_real))]
  #normalized by standard deviation of zero imputed matrix
  nrmse_sum$nrmse_mle[i] <- sqrt(mean((PXD029525_limm_mle[rownames(PXD029525_limm_mle) %in% rownames(PXD029525_limm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_limm_mle))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_svd[i] <- sqrt(mean((PXD029525_limm_svd[rownames(PXD029525_limm_svd) %in% rownames(PXD029525_limm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_limm_svd))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_bpca[i] <- sqrt(mean((PXD029525_limm_bpca[rownames(PXD029525_limm_bpca) %in% rownames(PXD029525_limm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_limm_bpca))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_lls[i] <- sqrt(mean((PXD029525_limm_lls[rownames(PXD029525_limm_lls) %in% rownames(PXD029525_limm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_limm_lls))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_rf[i] <- sqrt(mean((PXD029525_limm_rf[rownames(PXD029525_limm_rf) %in% rownames(PXD029525_limm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_limm_rf))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_cf[i] <- sqrt(mean((PXD029525_limm_cf[rownames(PXD029525_limm_cf) %in% rownames(PXD029525_limm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_limm_cf))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_vae[i] <- sqrt(mean((PXD029525_limm_vae[rownames(PXD029525_limm_vae) %in% rownames(PXD029525_limm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_limm_vae))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_dae[i] <- sqrt(mean((PXD029525_limm_dae[rownames(PXD029525_limm_dae) %in% rownames(PXD029525_limm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_limm_dae))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
}

nrmse_long <- nrmse_sum %>% pivot_longer(cols = 2:length(.), values_to = "NRMSE", names_to = "method")
nrmse_long$method <- factor(nrmse_long$method, levels = c("nrmse_bpca", "nrmse_knn", "nrmse_lls", "nrmse_mle", 
                                                          "nrmse_svd", "nrmse_rf", "nrmse_cf", "nrmse_dae", "nrmse_vae"))

#barplot
#calculate mean and standard error
nrmse_stat <- data.frame(method = colnames(nrmse_sum)[-c(1)],
                         mean = colMeans(nrmse_sum[,-c(1)]),
                         se = colSds(as.matrix(nrmse_sum[,-c(1)]))/sqrt(nrow(nrmse_sum)))
nrmse_stat$method <- factor(nrmse_stat$method, levels = c("nrmse_rf", "nrmse_bpca", "nrmse_knn", "nrmse_lls",  "nrmse_mle",
                                                          "nrmse_svd", "nrmse_cf", "nrmse_dae", "nrmse_vae"))
ggplot(data = nrmse_stat, aes(x = method, y = mean, fill = method)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = mean - 1.96*se, ymax = mean + 1.96*se, width = 0.2)) + #95% confidence interval
  geom_text(aes(label = round(mean, digits = 2)), vjust = -5, size = 5) +
  scale_fill_discrete(labels = c("RF", "BPCA", "kNN", "LLS", "MLE", "SVD", "CF", "DAE", "VAE")) +
  ylim(0, 0.8) +
  labs(x = "LIMM",
       y = "NRMSE") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        panel.grid.major.x = element_blank(),
        panel.grid.major = element_line(linewidth = 1),
        legend.position = "none")
colMeans(nrmse_sum[,2:length(nrmse_sum)])
quants <- c(0, 0.25, 0.50, 0.75, 1)
apply(nrmse_sum[,2:length(nrmse_sum)], 2 , quantile , probs = quants , na.rm = TRUE )

#lilm
comp_real <- PXD029525_lilm_real[which(rownames(PXD029525_lilm_real) %in% rownames(PXD029525_lilm_na)),] %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "real")
comp_na <- PXD029525_lilm_na %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "na")
comp_knn <- PXD029525_lilm_knn %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "knn")
comp_mle <- PXD029525_lilm_mle %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "mle")
comp_svd <- PXD029525_lilm_svd %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "svd")
comp_bpca <- PXD029525_lilm_bpca %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "bpca")
comp_lls <- PXD029525_lilm_lls %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "lls")
comp_rf <- PXD029525_lilm_rf %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "rf")
comp_cf <- PXD029525_lilm_cf %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "cf")
comp_vae <- PXD029525_lilm_vae %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "vae")
comp_dae <- PXD029525_lilm_dae %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "dae")

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
                   dae = comp_dae$dae)
comp_long <- comp %>% pivot_longer(., 2:length(.), names_to = "method", values_to = "value")

#NRMSE
nrmse_sum <- data.frame(sample = unique(comp_real$sample),
                        nrmse_knn = NA,
                        nrmse_mle = NA,
                        nrmse_svd = NA, 
                        nrmse_bpca = NA, 
                        nrmse_lls = NA,
                        nrmse_rf = NA,
                        nrmse_cf = NA,
                        nrmse_vae = NA,
                        nrmse_dae = NA)

for (i in 1:length(nrmse_sum$sample)) {
  z <- PXD029525_lilm_real[which(rownames(PXD029525_lilm_real) %in% rownames(PXD029525_lilm_na) & rownames(PXD029525_lilm_real) %in% rownames(PXD029525_lilm_vae)),grep(nrmse_sum$sample[i], colnames(PXD029525_lilm_real))]
  #normalized by standard deviation of zero imputed matrix
  nrmse_sum$nrmse_knn[i] <- sqrt(mean((PXD029525_lilm_knn[which(rownames(PXD029525_lilm_knn) %in% rownames(PXD029525_lilm_vae)),grep(nrmse_sum$sample[i], colnames(PXD029525_lilm_knn))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_mle[i] <- sqrt(mean((PXD029525_lilm_mle[rownames(PXD029525_lilm_mle) %in% rownames(PXD029525_lilm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_lilm_mle))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_svd[i] <- sqrt(mean((PXD029525_lilm_svd[rownames(PXD029525_lilm_svd) %in% rownames(PXD029525_lilm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_lilm_svd))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_bpca[i] <- sqrt(mean((PXD029525_lilm_bpca[rownames(PXD029525_lilm_bpca) %in% rownames(PXD029525_lilm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_lilm_bpca))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_lls[i] <- sqrt(mean((PXD029525_lilm_lls[rownames(PXD029525_lilm_lls) %in% rownames(PXD029525_lilm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_lilm_lls))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_rf[i] <- sqrt(mean((PXD029525_lilm_rf[rownames(PXD029525_lilm_rf) %in% rownames(PXD029525_lilm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_lilm_rf))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_cf[i] <- sqrt(mean((PXD029525_lilm_cf[rownames(PXD029525_lilm_cf) %in% rownames(PXD029525_lilm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_lilm_cf))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_vae[i] <- sqrt(mean((PXD029525_lilm_vae[rownames(PXD029525_lilm_vae) %in% rownames(PXD029525_lilm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_lilm_vae))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_dae[i] <- sqrt(mean((PXD029525_lilm_dae[rownames(PXD029525_lilm_dae) %in% rownames(PXD029525_lilm_vae),grep(nrmse_sum$sample[i], colnames(PXD029525_lilm_dae))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
}

nrmse_long <- nrmse_sum %>% pivot_longer(cols = 2:length(.), values_to = "NRMSE", names_to = "method")
nrmse_long$method <- factor(nrmse_long$method, levels = c("nrmse_bpca", "nrmse_knn", "nrmse_lls", "nrmse_mle", 
                                                          "nrmse_svd", "nrmse_rf", "nrmse_cf", "nrmse_dae", "nrmse_vae"))

#barplot
#calculate mean and standard error
nrmse_stat <- data.frame(method = colnames(nrmse_sum)[-c(1)],
                         mean = colMeans(nrmse_sum[,-c(1)]),
                         se = colSds(as.matrix(nrmse_sum[,-c(1)]))/sqrt(nrow(nrmse_sum)))
nrmse_stat$method <- factor(nrmse_stat$method, levels = c("nrmse_rf", "nrmse_bpca", "nrmse_knn", "nrmse_lls",  "nrmse_mle",
                                                          "nrmse_svd", "nrmse_cf", "nrmse_dae", "nrmse_vae"))
ggplot(data = nrmse_stat, aes(x = method, y = mean, fill = method)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = mean - 1.96*se, ymax = mean + 1.96*se, width = 0.2)) + #95% confidence interval
  geom_text(aes(label = round(mean, digits = 2)), vjust = -5, size = 5) +
  scale_fill_discrete(labels = c("RF", "BPCA", "kNN", "LLS", "MLE", "SVD", "CF", "DAE", "VAE")) +
  ylim(0, 1.5) +
  labs(x = "LILM",
       y = "NRMSE") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        panel.grid.major.x = element_blank(),
        panel.grid.major = element_line(linewidth = 1),
        legend.position = "none")
colMeans(nrmse_sum[,2:length(nrmse_sum)])
quants <- c(0, 0.25, 0.50, 0.75, 1)
apply(nrmse_sum[,2:length(nrmse_sum)], 2 , quantile , probs = quants , na.rm = TRUE )

bap_lilm <- ggplot(data = nrmse_stat, aes(x = method, y = mean, fill = method)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = mean - 1.96*se, ymax = mean + 1.96*se, width = 0.2)) + #95% confidence interval
  geom_text(aes(label = round(mean, digits = 2)), vjust = -5, size = 5) +
  scale_fill_discrete(labels = c("RF", "BPCA", "kNN", "LLS", "MLE", "SVD", "CF", "DAE", "VAE")) +
  ylim(0, 1.5) +
  labs(x = "LILM",
       y = "NRMSE") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        panel.grid.major.x = element_blank(),
        panel.grid.major = element_line(linewidth = 1), 
        legend.title = element_text(face = "bold"))
legend <- get_legend(bap_lilm)
grid.newpage()
grid.draw(legend)

#impute as whole
PXD029525_whole_na <- rbind(PXD029525_hihm_na, PXD029525_himm_na, PXD029525_hilm_na,
                            PXD029525_mihm_na, PXD029525_mimm_na, PXD029525_milm_na,
                            PXD029525_lihm_na, PXD029525_limm_na, PXD029525_lilm_na)
#real
PXD029525_whole_real <- rbind(PXD029525_hihm_real[rownames(PXD029525_hihm_real) %in% rownames(PXD029525_hihm_na),],
                              PXD029525_himm_real[rownames(PXD029525_himm_real) %in% rownames(PXD029525_himm_na),],
                              PXD029525_hilm_real[rownames(PXD029525_hilm_real) %in% rownames(PXD029525_hilm_na),],
                              PXD029525_mihm_real[rownames(PXD029525_mihm_real) %in% rownames(PXD029525_mihm_na),],
                              PXD029525_mimm_real[rownames(PXD029525_mimm_real) %in% rownames(PXD029525_mimm_na),],
                              PXD029525_milm_real[rownames(PXD029525_milm_real) %in% rownames(PXD029525_milm_na),],
                              PXD029525_lihm_real[rownames(PXD029525_lihm_real) %in% rownames(PXD029525_lihm_na),],
                              PXD029525_limm_real[rownames(PXD029525_limm_real) %in% rownames(PXD029525_limm_na),],
                              PXD029525_lilm_real[rownames(PXD029525_lilm_real) %in% rownames(PXD029525_lilm_na),])
#bpca
PXD029525_whole_bpca <- pca(PXD029525_whole_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
PXD029525_whole_bpca <- as.data.frame(PXD029525_whole_bpca@completeObs)
#svd
PXD029525_whole_svd <- as.data.frame(impute.wrapper.SVD(PXD029525_whole_na, K = 10))
#knn
PXD029525_whole_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(PXD029525_whole_na), K = 5))
#mle
PXD029525_whole_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(PXD029525_whole_na)))
#lls
PXD029525_whole_lls <- llsImpute(PXD029525_whole_na, k = 5, center = TRUE, completeObs = TRUE)
PXD029525_whole_lls <- as.data.frame(PXD029525_whole_lls@completeObs)
#rf
PXD029525_whole_rf <- missForest(PXD029525_whole_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
PXD029525_whole_rf <- PXD029525_whole_rf$ximp

#for pimms
write.csv(t(PXD029525_whole_na), file = "PXD029525_whole_na.csv")
#import pimms result
#fread for fast read
PXD029525_whole_cf <- fread("PXD029525_whole_cf.csv", header = TRUE)
PXD029525_whole_dae <- fread("PXD029525_whole_dae.csv", header = TRUE)
PXD029525_whole_vae <- fread("PXD029525_whole_vae.csv", header = TRUE)

PXD029525_whole_cf <- as.data.frame(PXD029525_whole_cf)
PXD029525_whole_dae <- as.data.frame(PXD029525_whole_dae)
PXD029525_whole_vae <- as.data.frame(PXD029525_whole_vae)

#set colname
colnames(PXD029525_whole_dae) <- PXD029525_whole_dae[1,]
colnames(PXD029525_whole_vae) <- PXD029525_whole_vae[1,]

#remove excessive rows
PXD029525_whole_dae <- PXD029525_whole_dae[-c(1,2),]
PXD029525_whole_vae <- PXD029525_whole_vae[-c(1,2),]

#turn chr into num
PXD029525_whole_dae[,-1] <- sapply(PXD029525_whole_dae[,-1], as.numeric)
PXD029525_whole_vae[,-1] <- sapply(PXD029525_whole_vae[,-1], as.numeric)

#set sample id as row name
rownames(PXD029525_whole_cf) <- PXD029525_whole_cf[,1]
rownames(PXD029525_whole_dae) <- PXD029525_whole_dae[,1]
rownames(PXD029525_whole_vae) <- PXD029525_whole_vae[,1]

#transpose
PXD029525_whole_cf <- as.data.frame(t(PXD029525_whole_cf[,-1]))
PXD029525_whole_dae <- as.data.frame(t(PXD029525_whole_dae[,-1]))
PXD029525_whole_vae <- as.data.frame(t(PXD029525_whole_vae[,-1]))

#assemble the optimal method matrix
PXD029525_whole_mix <- rbind(PXD029525_hihm_rf, PXD029525_himm_rf, PXD029525_hilm_rf,
                             PXD029525_mihm_bpca, PXD029525_mimm_rf, PXD029525_milm_rf,
                             PXD029525_lihm_bpca, PXD029525_limm_bpca, PXD029525_lilm_rf)

comp_real <- PXD029525_whole_real %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "real")
comp_na <- PXD029525_whole_na %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "na")
comp_knn <- PXD029525_whole_knn %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "knn")
comp_mle <- PXD029525_whole_mle %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "mle")
comp_svd <- PXD029525_whole_svd %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "svd")
comp_bpca <- PXD029525_whole_bpca %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "bpca")
comp_lls <- PXD029525_whole_lls %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "lls")
comp_rf <- PXD029525_whole_rf %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "rf")
comp_cf <- PXD029525_whole_cf %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "cf")
comp_vae <- PXD029525_whole_vae %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "vae")
comp_dae <- PXD029525_whole_dae %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "dae")
comp_mix <- PXD029525_whole_mix %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "mix")

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
comp_long <- comp %>% pivot_longer(., 2:length(.), names_to = "method", values_to = "value")

#NRMSE
nrmse_sum <- data.frame(sample = unique(comp_real$sample),
                        nrmse_knn = NA,
                        nrmse_mle = NA,
                        nrmse_svd = NA, 
                        nrmse_bpca = NA, 
                        nrmse_lls = NA,
                        nrmse_rf = NA,
                        nrmse_cf = NA,
                        nrmse_vae = NA,
                        nrmse_dae = NA,
                        nrmse_mix = NA)

for (i in 1:length(nrmse_sum$sample)) {
  z <- PXD029525_whole_real[,grep(nrmse_sum$sample[i], colnames(PXD029525_whole_real))]
  #normalized by standard deviation of zero imputed matrix
  nrmse_sum$nrmse_knn[i] <- sqrt(mean((PXD029525_whole_knn[,grep(nrmse_sum$sample[i], colnames(PXD029525_whole_knn))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_mle[i] <- sqrt(mean((PXD029525_whole_mle[,grep(nrmse_sum$sample[i], colnames(PXD029525_whole_mle))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_svd[i] <- sqrt(mean((PXD029525_whole_svd[,grep(nrmse_sum$sample[i], colnames(PXD029525_whole_svd))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_bpca[i] <- sqrt(mean((PXD029525_whole_bpca[,grep(nrmse_sum$sample[i], colnames(PXD029525_whole_bpca))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_lls[i] <- sqrt(mean((PXD029525_whole_lls[,grep(nrmse_sum$sample[i], colnames(PXD029525_whole_lls))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_rf[i] <- sqrt(mean((PXD029525_whole_rf[,grep(nrmse_sum$sample[i], colnames(PXD029525_whole_rf))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_cf[i] <- sqrt(mean((PXD029525_whole_cf[,grep(nrmse_sum$sample[i], colnames(PXD029525_whole_cf))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_vae[i] <- sqrt(mean((PXD029525_whole_vae[,grep(nrmse_sum$sample[i], colnames(PXD029525_whole_vae))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_dae[i] <- sqrt(mean((PXD029525_whole_dae[,grep(nrmse_sum$sample[i], colnames(PXD029525_whole_dae))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
  nrmse_sum$nrmse_mix[i] <- sqrt(mean((PXD029525_whole_mix[,grep(nrmse_sum$sample[i], colnames(PXD029525_whole_dae))] - z)^2, na.rm = TRUE))/sd(z, na.rm = TRUE)
}

nrmse_long <- nrmse_sum %>% pivot_longer(cols = 2:length(.), values_to = "NRMSE", names_to = "method")
nrmse_long$method <- factor(nrmse_long$method, levels = c("nrmse_bpca", "nrmse_knn", "nrmse_lls", "nrmse_mle", "nrmse_svd", "nrmse_rf", 
                                                          "nrmse_cf", "nrmse_dae", "nrmse_vae", "nrmse_mix"))
ggplot(data = nrmse_long, aes(x = method, y = NRMSE)) +
  geom_boxplot() +
  #geom_jitter(shape = 16, position = position_jitter(0.2)) +
  labs(title = "NRMSE of different missing value imputation methods \nWhole",
       x = "Imputation Methods") +
  scale_x_discrete(labels = c("BPCA", "kNN", "LLS", "MLE", "SVD", "RF", "CF", "DAE", "VAE", "Mix")) +
  theme_bw() +
  theme(axis.text.x = element_text(face = "bold", size = 10), 
        axis.title = element_text(size = 10),
        title = element_text(size = 12)) +
  ylim(0,0.8)
colMeans(nrmse_sum[,2:length(nrmse_sum)])
quants <- c(0, 0.25, 0.50, 0.75, 1)
apply(nrmse_sum[,2:length(nrmse_sum)], 2 , quantile , probs = quants , na.rm = TRUE )

#barplot
#calculate mean and standard error
#drop mle
nrmse_stat <- data.frame(method = colnames(nrmse_sum)[-c(1,3)],
                         mean = colMeans(nrmse_sum[,-c(1,3)]),
                         se = colSds(as.matrix(nrmse_sum[,-c(1,3)]))/sqrt(nrow(nrmse_sum)))
nrmse_stat$method <- factor(nrmse_stat$method, levels = c("nrmse_rf", "nrmse_bpca", "nrmse_knn", "nrmse_lls",  
                                                          "nrmse_svd", "nrmse_cf", "nrmse_dae", "nrmse_vae", "nrmse_mix"))

ggplot(data = nrmse_stat, aes(x = method, y = mean, fill = method)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = mean - 1.96*se, ymax = mean + 1.96*se, width = 0.2)) + #95% confidence interval
  geom_text(aes(label = round(mean, digits = 2)), vjust = -4, size = 5) +
  scale_fill_discrete(labels = c("RF", "BPCA", "kNN", "LLS", "SVD", "CF", "DAE", "VAE", "Mix")) +
  scale_x_discrete(labels = c("RF", "BPCA", "kNN", "LLS", "SVD", "CF", "DAE", "VAE", "Mix")) +
  ylim(0, 0.5) +
  labs(x = "Method",
       y = "NRMSE") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12),
        #axis.ticks.x = element_blank(),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        panel.grid.major.x = element_blank(),
        panel.grid.major = element_line(linewidth = 1),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12, face = "bold"))


##differential expression analysis
PXD029525_meta <- data.frame("sample" = colnames(PXD029525_whole_real),
                             "cell" = str_split_i(colnames(PXD029525_whole_real), " ", 3),
                             "treatment" = str_split_i(colnames(PXD029525_whole_real), " ", 4))

DE_result <- data.frame("df" = ls()[grep("PXD029525_whole_", ls())],
                        "method" = NA,
                        "sig_DEP" = NA,
                        "sig_DF" = NA)
DE_result$method <- str_split_i(DE_result$df, "_", 3)

#limma
design <- model.matrix(~ 0 + cell + treatment, data = PXD029525_meta)

for (i in 1:length(DE_result$df)) {
  data <- get(DE_result$df[i])
  rid <- rownames(data)
  cid <- colnames(data)
  data <- normalize.quantiles(as.matrix(data))
  rownames(data) <- rid
  colnames(data) <- cid
  fit <- lmFit(data, design)
  con <- makeContrasts(treatmentTRTM, levels = design)
  fit.con <- contrasts.fit(fit, con)
  bayes <- eBayes(fit.con)
  result <- topTable(bayes, sort.by="none", n=Inf, adjust.method = "BH")
  sig <- result[result$adj.P.Val < 0.05 & abs(result$logFC) > 1,]
  sig <- sig[!is.na(sig$logFC),]
  dfname <- paste(DE_result$df[i], "sigDEP", sep = "_")
  assign(dfname, sig)
  DE_result$sig_DEP[i] <- length(rownames(sig))
  DE_result$sig_DF[i] <- dfname
  print(paste(DE_result$df[i], length(sig$logFC), "sigDEP", sep = " "))
}

DE_result <- DE_result[DE_result$method != "na",]
DE_result$method <- factor(DE_result$method, levels = c("real", "rf", "bpca", "knn", "lls", "mle", 
                                                        "svd", "cf", "dae", "vae", "mix"))

ggplot(data = DE_result, aes(x = method, y = sig_DEP, fill = method)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sig_DEP, vjust = -0.5)) +
  scale_fill_discrete(labels = c("Real", "RF", "BPCA", "kNN", "LLS", "MLE", "SVD", "CF", "DAE", "VAE", "Mix")) +
  scale_x_discrete(labels = c("Real", "RF", "BPCA", "kNN", "LLS", "MLE", "SVD", "CF", "DAE", "VAE", "Mix")) +
  labs(x = "Method",
       y = "Number of significant differentially expressed peptides") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12, face = "bold"))


#PXD002882
PXD002882 <- read_delim("reference/Pietz.2024/data/datasets/PXD002882_crohns_disease/peptides.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
PXD002882_meta <- read_delim("reference/Pietz.2024/data/datasets/PXD002882_crohns_disease/experimentalDesign_Dec92014.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

#keep intensity
PXD002882_peptide <- as.data.frame(PXD002882[,c(1,230:322)])
#set row name
rownames(PXD002882_peptide) <- PXD002882_peptide[,1]
PXD002882_peptide[,-1] <- log2(PXD002882_peptide[,-1])
PXD002882_peptide <- PXD002882_peptide[,-1]
PXD002882_peptide[PXD002882_peptide == -Inf] <- NA

#calculate NA rate, mean intensity
PXD002882_peptide <- add_column(PXD002882_peptide, mean_intensity = NA, .before = "Intensity CD_CoA_230_Severe")
PXD002882_peptide <- add_column(PXD002882_peptide, NA_rate = NA, .after = "mean_intensity")

PXD002882_peptide$mean_intensity <- rowMeans(PXD002882_peptide[, -c(1,2)], na.rm = TRUE)
PXD002882_peptide$NA_rate <- rowSums(is.na(PXD002882_peptide[,-c(1,2)]))/(ncol(PXD002882_peptide) - 2)

PXD002882_int_low <- quantile(PXD002882_peptide$mean_intensity, na.rm = TRUE)[2]
PXD002882_int_high <- quantile(PXD002882_peptide$mean_intensity, na.rm = TRUE)[4]

#NA rate vs intensity
ggplot(data = PXD002882_peptide, aes(x = mean_intensity, y = NA_rate)) + 
  geom_point(alpha = 0.1, shape = 20) +
  labs(title = paste("C. PXD002882 Crohn's disease", "Peptide mean intensity and missing rate", sep = "\n"),
       x = bquote(paste("mean log"["2"]*"(intensity)")),
       y = "Missing rate") +
  lims(x = c(15,30),
       y = c(0,1)) +
  theme_classic() +
  theme(axis.text = element_text(size = 15, color = "black"),
        axis.title = element_text(size = 15),
        title = element_text(size = 15)) 
last_plot() +
  stat_density2d(aes(fill = ..density..), geom = "tile", contour = FALSE, n=200, alpha = 0.8) +
  scale_fill_continuous(low = "white", high = "dodgerblue4") +
  stat_smooth(method = "lm", geom = "smooth", formula = y~x, se = TRUE, color = "black", size = 0.5) +
  labs(fill = "Density") +
  geom_line(aes(x = PXD002882_int_high), color = "red", size = 0.5) +
  geom_line(aes(x = PXD002882_int_low), color = "red", size = 0.5) +
  geom_line(aes(y = 0.75), color = "red", size = 0.5) +
  geom_line(aes(y = 0.25), color = "red", size = 0.5) +
  theme(legend.title = element_text(size = 13),
        legend.text = element_text(size = 12))
ggsave("figure6_int_mr_PXD002282.tiff",
       plot = last_plot(),
       path = "figures/manuscript/",
       device = "tiff",
       width = 2244,
       height = 1496,
       units = "px",
       dpi = 300)
ggsave("figure6_int_mr_PXD002282.jpeg",
       plot = last_plot(),
       path = "figures/manuscript/",
       device = "jpeg",
       width = 2244,
       height = 1496,
       units = "px",
       dpi = 300)
summary(lm(NA_rate~mean_intensity, data = PXD002882_peptide))






#high intensity, high missing peptides
PXD002882_hihm <- PXD002882_peptide[which(PXD002882_peptide$mean_intensity >= PXD002882_int_high & 
                                            PXD002882_peptide$NA_rate >= 0.75 &
                                            PXD002882_peptide$NA_rate < 1),]
#high intensity, medium missing peptides
PXD002882_himm <- PXD002882_peptide[which(PXD002882_peptide$mean_intensity >= PXD002882_int_high & 
                                            PXD002882_peptide$NA_rate >= 0.25 &
                                            PXD002882_peptide$NA_rate < 0.75),]
#high intensity, low missing peptides
PXD002882_hilm <- PXD002882_peptide[which(PXD002882_peptide$mean_intensity >= PXD002882_int_high & 
                                            PXD002882_peptide$NA_rate < 0.25),]
#medium intensity, high missing peptides
PXD002882_mihm <- PXD002882_peptide[which(PXD002882_peptide$mean_intensity < PXD002882_int_high & 
                                            PXD002882_peptide$mean_intensity >= PXD002882_int_low &
                                            PXD002882_peptide$NA_rate >= 0.75 &
                                            PXD002882_peptide$NA_rate < 1),]
#medium intensity, medium missing peptides
PXD002882_mimm <- PXD002882_peptide[which(PXD002882_peptide$mean_intensity < PXD002882_int_high & 
                                            PXD002882_peptide$mean_intensity >= PXD002882_int_low &
                                            PXD002882_peptide$NA_rate >= 0.25 &
                                            PXD002882_peptide$NA_rate < 0.75),]
#medium intensity, low missing peptides
PXD002882_milm <- PXD002882_peptide[which(PXD002882_peptide$mean_intensity < PXD002882_int_high & 
                                            PXD002882_peptide$mean_intensity >= PXD002882_int_low &
                                            PXD002882_peptide$NA_rate < 0.25),]
#low intensity, high missing peptides
PXD002882_lihm <- PXD002882_peptide[which(PXD002882_peptide$mean_intensity < PXD002882_int_low & 
                                            PXD002882_peptide$NA_rate >= 0.75 &
                                            PXD002882_peptide$NA_rate < 1),]
#low intensity, medium missing peptides
PXD002882_limm <- PXD002882_peptide[which(PXD002882_peptide$mean_intensity < PXD002882_int_low & 
                                            PXD002882_peptide$NA_rate >= 0.25 &
                                            PXD002882_peptide$NA_rate < 0.75),]
#low intensity, low missing peptides
PXD002882_lilm <- PXD002882_peptide[which(PXD002882_peptide$mean_intensity < PXD002882_int_low & 
                                            PXD002882_peptide$NA_rate < 0.25),]

#keep intensity
PXD002882_hihm <- PXD002882_hihm[,-c(1,2)]
PXD002882_himm <- PXD002882_himm[,-c(1,2)]
PXD002882_hilm <- PXD002882_hilm[,-c(1,2)]
PXD002882_mihm <- PXD002882_mihm[,-c(1,2)]
PXD002882_mimm <- PXD002882_mimm[,-c(1,2)]
PXD002882_milm <- PXD002882_milm[,-c(1,2)]
PXD002882_lihm <- PXD002882_lihm[,-c(1,2)]
PXD002882_limm <- PXD002882_limm[,-c(1,2)]
PXD002882_lilm <- PXD002882_lilm[,-c(1,2)]

#prepare datasets for imputation, remove sample and pepetide that are all NA
PXD002882_hihm_na <- PXD002882_hihm[which(rowSums(is.na(PXD002882_hihm)) < ncol(PXD002882_hihm)),
                                    which(colSums(is.na(PXD002882_hihm)) < nrow(PXD002882_hihm))]
PXD002882_himm_na <- PXD002882_himm[which(rowSums(is.na(PXD002882_himm)) < ncol(PXD002882_himm)),
                                    which(colSums(is.na(PXD002882_himm)) < nrow(PXD002882_himm))]
PXD002882_hilm_na <- PXD002882_hilm[which(rowSums(is.na(PXD002882_hilm)) < ncol(PXD002882_hilm)),
                                    which(colSums(is.na(PXD002882_hilm)) < nrow(PXD002882_hilm))]
PXD002882_mihm_na <- PXD002882_mihm[which(rowSums(is.na(PXD002882_mihm)) < ncol(PXD002882_mihm)),
                                    which(colSums(is.na(PXD002882_mihm)) < nrow(PXD002882_mihm))]
PXD002882_mimm_na <- PXD002882_mimm[which(rowSums(is.na(PXD002882_mimm)) < ncol(PXD002882_mimm)),
                                    which(colSums(is.na(PXD002882_mimm)) < nrow(PXD002882_mimm))]
PXD002882_milm_na <- PXD002882_milm[which(rowSums(is.na(PXD002882_milm)) < ncol(PXD002882_milm)),
                                    which(colSums(is.na(PXD002882_milm)) < nrow(PXD002882_milm))]
PXD002882_lihm_na <- PXD002882_lihm[which(rowSums(is.na(PXD002882_lihm)) < ncol(PXD002882_lihm)),
                                    which(colSums(is.na(PXD002882_lihm)) < nrow(PXD002882_lihm))]
PXD002882_limm_na <- PXD002882_limm[which(rowSums(is.na(PXD002882_limm)) < ncol(PXD002882_limm)),
                                    which(colSums(is.na(PXD002882_limm)) < nrow(PXD002882_limm))]
PXD002882_lilm_na <- PXD002882_lilm[which(rowSums(is.na(PXD002882_lilm)) < ncol(PXD002882_lilm)),
                                    which(colSums(is.na(PXD002882_lilm)) < nrow(PXD002882_lilm))]

#real datasets
PXD002882_hihm_real <- PXD002882_hihm_na
PXD002882_himm_real <- PXD002882_himm_na
PXD002882_hilm_real <- PXD002882_hilm_na
PXD002882_mihm_real <- PXD002882_mihm_na
PXD002882_mimm_real <- PXD002882_mimm_na
PXD002882_milm_real <- PXD002882_milm_na
PXD002882_lihm_real <- PXD002882_lihm_na
PXD002882_limm_real <- PXD002882_limm_na
PXD002882_lilm_real <- PXD002882_lilm_na

#generate random sample points, 20%
which(is.na(PXD002882_hihm_na) == FALSE)

na_pos_PXD002882_hihm <- sample(which(is.na(PXD002882_hihm_na) == FALSE), 
                                0.2*length(which(is.na(PXD002882_hihm_na) == FALSE)))
na_pos_PXD002882_himm <- sample(which(is.na(PXD002882_himm_na) == FALSE), 
                                0.2*length(which(is.na(PXD002882_himm_na) == FALSE)))
na_pos_PXD002882_hilm <- sample(which(is.na(PXD002882_hilm_na) == FALSE), 
                                0.2*length(which(is.na(PXD002882_hilm_na) == FALSE)))
na_pos_PXD002882_mihm <- sample(which(is.na(PXD002882_mihm_na) == FALSE), 
                                0.2*length(which(is.na(PXD002882_mihm_na) == FALSE)))
na_pos_PXD002882_mimm <- sample(which(is.na(PXD002882_mimm_na) == FALSE), 
                                0.2*length(which(is.na(PXD002882_mimm_na) == FALSE)))
na_pos_PXD002882_milm <- sample(which(is.na(PXD002882_milm_na) == FALSE), 
                                0.2*length(which(is.na(PXD002882_milm_na) == FALSE)))
na_pos_PXD002882_lihm <- sample(which(is.na(PXD002882_lihm_na) == FALSE), 
                                0.2*length(which(is.na(PXD002882_lihm_na) == FALSE)))
na_pos_PXD002882_limm <- sample(which(is.na(PXD002882_limm_na) == FALSE), 
                                0.2*length(which(is.na(PXD002882_limm_na) == FALSE)))
na_pos_PXD002882_lilm <- sample(which(is.na(PXD002882_lilm_na) == FALSE), 
                                0.2*length(which(is.na(PXD002882_lilm_na) == FALSE)))

#replace true values with NA
PXD002882_hihm_na <- as.matrix(PXD002882_hihm_na)
PXD002882_hihm_na[na_pos_PXD002882_hihm] <- NA
PXD002882_hihm_na <- as.data.frame(PXD002882_hihm_na)
PXD002882_himm_na <- as.matrix(PXD002882_himm_na)
PXD002882_himm_na[na_pos_PXD002882_himm] <- NA
PXD002882_himm_na <- as.data.frame(PXD002882_himm_na)
PXD002882_hilm_na <- as.matrix(PXD002882_hilm_na)
PXD002882_hilm_na[na_pos_PXD002882_hilm] <- NA
PXD002882_hilm_na <- as.data.frame(PXD002882_hilm_na)
PXD002882_mihm_na <- as.matrix(PXD002882_mihm_na)
PXD002882_mihm_na[na_pos_PXD002882_mihm] <- NA
PXD002882_mihm_na <- as.data.frame(PXD002882_mihm_na)
PXD002882_mimm_na <- as.matrix(PXD002882_mimm_na)
PXD002882_mimm_na[na_pos_PXD002882_mimm] <- NA
PXD002882_mimm_na <- as.data.frame(PXD002882_mimm_na)
PXD002882_milm_na <- as.matrix(PXD002882_milm_na)
PXD002882_milm_na[na_pos_PXD002882_milm] <- NA
PXD002882_milm_na <- as.data.frame(PXD002882_milm_na)
PXD002882_lihm_na <- as.matrix(PXD002882_lihm_na)
PXD002882_lihm_na[na_pos_PXD002882_lihm] <- NA
PXD002882_lihm_na <- as.data.frame(PXD002882_lihm_na)
PXD002882_limm_na <- as.matrix(PXD002882_limm_na)
PXD002882_limm_na[na_pos_PXD002882_limm] <- NA
PXD002882_limm_na <- as.data.frame(PXD002882_limm_na)
PXD002882_lilm_na <- as.matrix(PXD002882_lilm_na)
PXD002882_lilm_na[na_pos_PXD002882_lilm] <- NA
PXD002882_lilm_na <- as.data.frame(PXD002882_lilm_na)

#change na pos: hihm, mihm, lihm, lilm
#count na percentage
PXD002882_hihm_na <- PXD002882_hihm_na[-which(rowSums(is.na(PXD002882_hihm_na)) == ncol(PXD002882_hihm_na)),]
PXD002882_mihm_na <- PXD002882_mihm_na[-which(rowSums(is.na(PXD002882_mihm_na)) == ncol(PXD002882_mihm_na)),]
PXD002882_lihm_na <- PXD002882_lihm_na[-which(rowSums(is.na(PXD002882_lihm_na)) == ncol(PXD002882_lihm_na)),]
#ignore lilm

#impute by intensity and missing rate
#bpca, svd, knn (SEQKNN, TRKNN, KNNIMPUTE), median, missForrest, RSN
#bpca
PXD002882_hihm_bpca <- pca(PXD002882_hihm_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
PXD002882_himm_bpca <- pca(PXD002882_himm_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
PXD002882_hilm_bpca <- pca(PXD002882_hilm_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
PXD002882_mihm_bpca <- pca(PXD002882_mihm_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
PXD002882_mimm_bpca <- pca(PXD002882_mimm_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
PXD002882_milm_bpca <- pca(PXD002882_milm_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
PXD002882_lihm_bpca <- pca(PXD002882_lihm_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
PXD002882_limm_bpca <- pca(PXD002882_limm_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
PXD002882_hihm_bpca <- as.data.frame(PXD002882_hihm_bpca@completeObs)
PXD002882_himm_bpca <- as.data.frame(PXD002882_himm_bpca@completeObs)
PXD002882_hilm_bpca <- as.data.frame(PXD002882_hilm_bpca@completeObs)
PXD002882_mihm_bpca <- as.data.frame(PXD002882_mihm_bpca@completeObs)
PXD002882_mimm_bpca <- as.data.frame(PXD002882_mimm_bpca@completeObs)
PXD002882_milm_bpca <- as.data.frame(PXD002882_milm_bpca@completeObs)
PXD002882_lihm_bpca <- as.data.frame(PXD002882_lihm_bpca@completeObs)
PXD002882_limm_bpca <- as.data.frame(PXD002882_limm_bpca@completeObs)

#svd
PXD002882_hihm_svd <- as.data.frame(impute.wrapper.SVD(PXD002882_hihm_na, K = 10))
PXD002882_himm_svd <- as.data.frame(impute.wrapper.SVD(PXD002882_himm_na, K = 10))
PXD002882_hilm_svd <- as.data.frame(impute.wrapper.SVD(PXD002882_hilm_na, K = 10))
PXD002882_mihm_svd <- as.data.frame(impute.wrapper.SVD(PXD002882_mihm_na, K = 10))
PXD002882_mimm_svd <- as.data.frame(impute.wrapper.SVD(PXD002882_mimm_na, K = 10))
PXD002882_milm_svd <- as.data.frame(impute.wrapper.SVD(PXD002882_milm_na, K = 10))
PXD002882_lihm_svd <- as.data.frame(impute.wrapper.SVD(PXD002882_lihm_na, K = 10))
PXD002882_limm_svd <- as.data.frame(impute.wrapper.SVD(PXD002882_limm_na, K = 10))

#knn
PXD002882_hihm_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(PXD002882_hihm_na), K = 5))
PXD002882_himm_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(PXD002882_himm_na), K = 5))
PXD002882_hilm_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(PXD002882_hilm_na), K = 5))
PXD002882_mihm_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(PXD002882_mihm_na), K = 5))
PXD002882_mimm_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(PXD002882_mimm_na), K = 5))
PXD002882_milm_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(PXD002882_milm_na), K = 5))
PXD002882_lihm_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(PXD002882_lihm_na), K = 5))
PXD002882_limm_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(PXD002882_limm_na), K = 5))

#mle
PXD002882_hihm_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(PXD002882_hihm_na)))
PXD002882_himm_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(PXD002882_himm_na)))
PXD002882_hilm_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(PXD002882_hilm_na)))
PXD002882_mihm_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(PXD002882_mihm_na)))
PXD002882_mimm_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(PXD002882_mimm_na)))
PXD002882_milm_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(PXD002882_milm_na)))
PXD002882_lihm_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(PXD002882_lihm_na)))
PXD002882_limm_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(PXD002882_limm_na)))

#lls
PXD002882_hihm_lls <- llsImpute(PXD002882_hihm_na, k = 5, center = TRUE, completeObs = TRUE)
PXD002882_himm_lls <- llsImpute(PXD002882_himm_na, k = 5, center = TRUE, completeObs = TRUE)
PXD002882_hilm_lls <- llsImpute(PXD002882_hilm_na, k = 5, center = TRUE, completeObs = TRUE)
PXD002882_mihm_lls <- llsImpute(PXD002882_mihm_na, k = 5, center = TRUE, completeObs = TRUE)
PXD002882_mimm_lls <- llsImpute(PXD002882_mimm_na, k = 5, center = TRUE, completeObs = TRUE)
PXD002882_milm_lls <- llsImpute(PXD002882_milm_na, k = 5, center = TRUE, completeObs = TRUE)
PXD002882_lihm_lls <- llsImpute(PXD002882_lihm_na, k = 5, center = TRUE, completeObs = TRUE)
PXD002882_limm_lls <- llsImpute(PXD002882_limm_na, k = 5, center = TRUE, completeObs = TRUE)
PXD002882_hihm_lls <- as.data.frame(PXD002882_hihm_lls@completeObs)
PXD002882_himm_lls <- as.data.frame(PXD002882_himm_lls@completeObs)
PXD002882_hilm_lls <- as.data.frame(PXD002882_hilm_lls@completeObs)
PXD002882_mihm_lls <- as.data.frame(PXD002882_mihm_lls@completeObs)
PXD002882_mimm_lls <- as.data.frame(PXD002882_mimm_lls@completeObs)
PXD002882_milm_lls <- as.data.frame(PXD002882_milm_lls@completeObs)
PXD002882_lihm_lls <- as.data.frame(PXD002882_lihm_lls@completeObs)
PXD002882_limm_lls <- as.data.frame(PXD002882_limm_lls@completeObs)

#rf
#register number of cores for parallel
registerDoParallel(cores = 20)
PXD002882_hihm_rf <- missForest(PXD002882_hihm_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
PXD002882_himm_rf <- missForest(PXD002882_himm_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
PXD002882_hilm_rf <- missForest(PXD002882_hilm_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
PXD002882_mihm_rf <- missForest(PXD002882_mihm_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
PXD002882_mimm_rf <- missForest(PXD002882_mimm_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
PXD002882_milm_rf <- missForest(PXD002882_milm_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
PXD002882_lihm_rf <- missForest(PXD002882_lihm_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
PXD002882_limm_rf <- missForest(PXD002882_limm_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
PXD002882_hihm_rf <- PXD002882_hihm_rf$ximp
PXD002882_himm_rf <- PXD002882_himm_rf$ximp
PXD002882_hilm_rf <- PXD002882_hilm_rf$ximp
PXD002882_mihm_rf <- PXD002882_mihm_rf$ximp
PXD002882_mimm_rf <- PXD002882_mimm_rf$ximp
PXD002882_milm_rf <- PXD002882_milm_rf$ximp
PXD002882_lihm_rf <- PXD002882_lihm_rf$ximp
PXD002882_limm_rf <- PXD002882_limm_rf$ximp

#transpose data frame for pimms
write.csv(t(PXD002882_hihm_na), file = "PXD002882_hihm_na.csv")
write.csv(t(PXD002882_himm_na), file = "PXD002882_himm_na.csv")
write.csv(t(PXD002882_hilm_na), file = "PXD002882_hilm_na.csv")
write.csv(t(PXD002882_mihm_na), file = "PXD002882_mihm_na.csv")
write.csv(t(PXD002882_mimm_na), file = "PXD002882_mimm_na.csv")
write.csv(t(PXD002882_milm_na), file = "PXD002882_milm_na.csv")
write.csv(t(PXD002882_lihm_na), file = "PXD002882_lihm_na.csv")
write.csv(t(PXD002882_limm_na), file = "PXD002882_limm_na.csv")

#import pimms result
#cf
PXD002882_hihm_cf <- read.csv("PXD002882_hihm_cf.csv")
PXD002882_himm_cf <- read.csv("PXD002882_himm_cf.csv")
PXD002882_hilm_cf <- read.csv("PXD002882_hilm_cf.csv")
PXD002882_mihm_cf <- read.csv("PXD002882_mihm_cf.csv")
PXD002882_mimm_cf <- read.csv("PXD002882_mimm_cf.csv")
PXD002882_milm_cf <- read.csv("PXD002882_milm_cf.csv")
PXD002882_lihm_cf <- read.csv("PXD002882_lihm_cf.csv")
PXD002882_limm_cf <- read.csv("PXD002882_limm_cf.csv")
#vae
PXD002882_hihm_vae <- read.csv("PXD002882_hihm_vae.csv")
PXD002882_himm_vae <- read.csv("PXD002882_himm_vae.csv")
PXD002882_hilm_vae <- read.csv("PXD002882_hilm_vae.csv")
PXD002882_mihm_vae <- read.csv("PXD002882_mihm_vae.csv")
PXD002882_mimm_vae <- read.csv("PXD002882_mimm_vae.csv")
PXD002882_milm_vae <- read.csv("PXD002882_milm_vae.csv")
PXD002882_lihm_vae <- read.csv("PXD002882_lihm_vae.csv")
PXD002882_limm_vae <- read.csv("PXD002882_limm_vae.csv")
#dae
PXD002882_hihm_dae <- read.csv("PXD002882_hihm_dae.csv")
PXD002882_himm_dae <- read.csv("PXD002882_himm_dae.csv")
PXD002882_hilm_dae <- read.csv("PXD002882_hilm_dae.csv")
PXD002882_mihm_dae <- read.csv("PXD002882_mihm_dae.csv")
PXD002882_mimm_dae <- read.csv("PXD002882_mimm_dae.csv")
PXD002882_milm_dae <- read.csv("PXD002882_milm_dae.csv")
PXD002882_lihm_dae <- read.csv("PXD002882_lihm_dae.csv")
PXD002882_limm_dae <- read.csv("PXD002882_limm_dae.csv")

#set colname
colnames(PXD002882_hihm_dae) <- PXD002882_hihm_dae[1,]
colnames(PXD002882_himm_dae) <- PXD002882_himm_dae[1,]
colnames(PXD002882_hilm_dae) <- PXD002882_hilm_dae[1,]
colnames(PXD002882_mihm_dae) <- PXD002882_mihm_dae[1,]
colnames(PXD002882_mimm_dae) <- PXD002882_mimm_dae[1,]
colnames(PXD002882_milm_dae) <- PXD002882_milm_dae[1,]
colnames(PXD002882_lihm_dae) <- PXD002882_lihm_dae[1,]
colnames(PXD002882_limm_dae) <- PXD002882_limm_dae[1,]

colnames(PXD002882_hihm_vae) <- PXD002882_hihm_vae[1,]
colnames(PXD002882_himm_vae) <- PXD002882_himm_vae[1,]
colnames(PXD002882_hilm_vae) <- PXD002882_hilm_vae[1,]
colnames(PXD002882_mihm_vae) <- PXD002882_mihm_vae[1,]
colnames(PXD002882_mimm_vae) <- PXD002882_mimm_vae[1,]
colnames(PXD002882_milm_vae) <- PXD002882_milm_vae[1,]
colnames(PXD002882_lihm_vae) <- PXD002882_lihm_vae[1,]
colnames(PXD002882_limm_vae) <- PXD002882_limm_vae[1,]

#remove excessive rows
PXD002882_hihm_dae <- PXD002882_hihm_dae[-c(1,2),]
PXD002882_himm_dae <- PXD002882_himm_dae[-c(1,2),]
PXD002882_hilm_dae <- PXD002882_hilm_dae[-c(1,2),]
PXD002882_mihm_dae <- PXD002882_mihm_dae[-c(1,2),]
PXD002882_mimm_dae <- PXD002882_mimm_dae[-c(1,2),]
PXD002882_milm_dae <- PXD002882_milm_dae[-c(1,2),]
PXD002882_lihm_dae <- PXD002882_lihm_dae[-c(1,2),]
PXD002882_limm_dae <- PXD002882_limm_dae[-c(1,2),]

PXD002882_hihm_vae <- PXD002882_hihm_vae[-c(1,2),]
PXD002882_himm_vae <- PXD002882_himm_vae[-c(1,2),]
PXD002882_hilm_vae <- PXD002882_hilm_vae[-c(1,2),]
PXD002882_mihm_vae <- PXD002882_mihm_vae[-c(1,2),]
PXD002882_mimm_vae <- PXD002882_mimm_vae[-c(1,2),]
PXD002882_milm_vae <- PXD002882_milm_vae[-c(1,2),]
PXD002882_lihm_vae <- PXD002882_lihm_vae[-c(1,2),]
PXD002882_limm_vae <- PXD002882_limm_vae[-c(1,2),]

#turn chr into num
PXD002882_hihm_dae[,-1] <- sapply(PXD002882_hihm_dae[,-1], as.numeric)
PXD002882_himm_dae[,-1] <- sapply(PXD002882_himm_dae[,-1], as.numeric)
PXD002882_hilm_dae[,-1] <- sapply(PXD002882_hilm_dae[,-1], as.numeric)
PXD002882_mihm_dae[,-1] <- sapply(PXD002882_mihm_dae[,-1], as.numeric)
PXD002882_mimm_dae[,-1] <- sapply(PXD002882_mimm_dae[,-1], as.numeric)
PXD002882_milm_dae[,-1] <- sapply(PXD002882_milm_dae[,-1], as.numeric)
PXD002882_lihm_dae[,-1] <- sapply(PXD002882_lihm_dae[,-1], as.numeric)
PXD002882_limm_dae[,-1] <- sapply(PXD002882_limm_dae[,-1], as.numeric)

PXD002882_hihm_vae[,-1] <- sapply(PXD002882_hihm_vae[,-1], as.numeric)
PXD002882_himm_vae[,-1] <- sapply(PXD002882_himm_vae[,-1], as.numeric)
PXD002882_hilm_vae[,-1] <- sapply(PXD002882_hilm_vae[,-1], as.numeric)
PXD002882_mihm_vae[,-1] <- sapply(PXD002882_mihm_vae[,-1], as.numeric)
PXD002882_mimm_vae[,-1] <- sapply(PXD002882_mimm_vae[,-1], as.numeric)
PXD002882_milm_vae[,-1] <- sapply(PXD002882_milm_vae[,-1], as.numeric)
PXD002882_lihm_vae[,-1] <- sapply(PXD002882_lihm_vae[,-1], as.numeric)
PXD002882_limm_vae[,-1] <- sapply(PXD002882_limm_vae[,-1], as.numeric)

#set sample id as row name
#cf
rownames(PXD002882_hihm_cf) <- PXD002882_hihm_cf[,1]
rownames(PXD002882_himm_cf) <- PXD002882_himm_cf[,1]
rownames(PXD002882_hilm_cf) <- PXD002882_hilm_cf[,1]
rownames(PXD002882_mihm_cf) <- PXD002882_mihm_cf[,1]
rownames(PXD002882_mimm_cf) <- PXD002882_mimm_cf[,1]
rownames(PXD002882_milm_cf) <- PXD002882_milm_cf[,1]
rownames(PXD002882_lihm_cf) <- PXD002882_lihm_cf[,1]
rownames(PXD002882_limm_cf) <- PXD002882_limm_cf[,1]
#vae
rownames(PXD002882_hihm_vae) <- PXD002882_hihm_vae[,1]
rownames(PXD002882_himm_vae) <- PXD002882_himm_vae[,1]
rownames(PXD002882_hilm_vae) <- PXD002882_hilm_vae[,1]
rownames(PXD002882_mihm_vae) <- PXD002882_mihm_vae[,1]
rownames(PXD002882_mimm_vae) <- PXD002882_mimm_vae[,1]
rownames(PXD002882_milm_vae) <- PXD002882_milm_vae[,1]
rownames(PXD002882_lihm_vae) <- PXD002882_lihm_vae[,1]
rownames(PXD002882_limm_vae) <- PXD002882_limm_vae[,1]
#dae
rownames(PXD002882_hihm_dae) <- PXD002882_hihm_dae[,1]
rownames(PXD002882_himm_dae) <- PXD002882_himm_dae[,1]
rownames(PXD002882_hilm_dae) <- PXD002882_hilm_dae[,1]
rownames(PXD002882_mihm_dae) <- PXD002882_mihm_dae[,1]
rownames(PXD002882_mimm_dae) <- PXD002882_mimm_dae[,1]
rownames(PXD002882_milm_dae) <- PXD002882_milm_dae[,1]
rownames(PXD002882_lihm_dae) <- PXD002882_lihm_dae[,1]
rownames(PXD002882_limm_dae) <- PXD002882_limm_dae[,1]

#transpose
#cf
PXD002882_hihm_cf <- as.data.frame(t(PXD002882_hihm_cf[,-1]))
PXD002882_himm_cf <- as.data.frame(t(PXD002882_himm_cf[,-1]))
PXD002882_hilm_cf <- as.data.frame(t(PXD002882_hilm_cf[,-1]))
PXD002882_mihm_cf <- as.data.frame(t(PXD002882_mihm_cf[,-1]))
PXD002882_mimm_cf <- as.data.frame(t(PXD002882_mimm_cf[,-1]))
PXD002882_milm_cf <- as.data.frame(t(PXD002882_milm_cf[,-1]))
PXD002882_lihm_cf <- as.data.frame(t(PXD002882_lihm_cf[,-1]))
PXD002882_limm_cf <- as.data.frame(t(PXD002882_limm_cf[,-1]))
#vae
PXD002882_hihm_vae <- as.data.frame(t(PXD002882_hihm_vae[,-1]))
PXD002882_himm_vae <- as.data.frame(t(PXD002882_himm_vae[,-1]))
PXD002882_hilm_vae <- as.data.frame(t(PXD002882_hilm_vae[,-1]))
PXD002882_mihm_vae <- as.data.frame(t(PXD002882_mihm_vae[,-1]))
PXD002882_mimm_vae <- as.data.frame(t(PXD002882_mimm_vae[,-1]))
PXD002882_milm_vae <- as.data.frame(t(PXD002882_milm_vae[,-1]))
PXD002882_lihm_vae <- as.data.frame(t(PXD002882_lihm_vae[,-1]))
PXD002882_limm_vae <- as.data.frame(t(PXD002882_limm_vae[,-1]))
#dae
PXD002882_hihm_dae <- as.data.frame(t(PXD002882_hihm_dae[,-1]))
PXD002882_himm_dae <- as.data.frame(t(PXD002882_himm_dae[,-1]))
PXD002882_hilm_dae <- as.data.frame(t(PXD002882_hilm_dae[,-1]))
PXD002882_mihm_dae <- as.data.frame(t(PXD002882_mihm_dae[,-1]))
PXD002882_mimm_dae <- as.data.frame(t(PXD002882_mimm_dae[,-1]))
PXD002882_milm_dae <- as.data.frame(t(PXD002882_milm_dae[,-1]))
PXD002882_lihm_dae <- as.data.frame(t(PXD002882_lihm_dae[,-1]))
PXD002882_limm_dae <- as.data.frame(t(PXD002882_limm_dae[,-1]))


#impute as whole
PXD002882_whole_na <- rbind(PXD002882_hihm_na, PXD002882_himm_na, PXD002882_hilm_na,
                            PXD002882_mihm_na, PXD002882_mimm_na, PXD002882_milm_na,
                            PXD002882_lihm_na, PXD002882_limm_na, PXD002882_lilm_na)
write.csv(PXD002882_whole_na, "PXD002882_whole_na.csv")


id <- "PXD002882"
regions <- c("hihm", "himm", "hilm", "mihm", "mimm", "milm", "lihm", "limm")
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
      geom_text(aes(label = round(mean, digits = 2)), vjust = -1, size = 5) +
      scale_fill_brewer(palette = "Set3", labels = c("RF", "BPCA", "kNN", "LLS", "MLE", "SVD", "CF", "DAE", "VAE")) +
      scale_y_break(c(2, 200), scales = 0.5) +
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
      geom_text(aes(label = round(mean, digits = 2)), vjust = -1.8, size = 5) +
      scale_fill_brewer(palette = "Set3", labels = c("RF", "BPCA", "kNN", "LLS", "MLE", "SVD", "CF", "DAE", "VAE")) +
      ylim(0,1.6) +
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
  ggsave(paste("figures/manuscript/PXD002882/NRMSE_barplot_", regions[i], ".jpeg", sep = ""), device = "jpeg", width = 1400, height = 1400, units = "px", dpi = 300)
  print(regions[i])
  print(colMeans(nrmse_sum[,2:length(nrmse_sum)], na.rm = TRUE))
  quants <- c(0, 0.25, 0.50, 0.75, 1)
  print(apply(nrmse_sum[,2:length(nrmse_sum)], 2 , quantile , probs = quants , na.rm = TRUE))
  print(paste(gsub("nrmse_", "", nrmse_stat$method[which(nrmse_stat$mean==min(nrmse_stat$mean, na.rm = TRUE))]),
              "is the imputation method with lowest NRMSE in",
              regions[i],
              sep = " "))
}


#impute as whole
PXD002882_whole_na <- rbind(PXD002882_hihm_na, PXD002882_himm_na, PXD002882_hilm_na,
                            PXD002882_mihm_na, PXD002882_mimm_na, PXD002882_milm_na,
                            PXD002882_lihm_na, PXD002882_limm_na)
#real
PXD002882_whole_real <- rbind(PXD002882_hihm_real[rownames(PXD002882_hihm_real) %in% rownames(PXD002882_hihm_na),],
                              PXD002882_himm_real[rownames(PXD002882_himm_real) %in% rownames(PXD002882_himm_na),],
                              PXD002882_hilm_real[rownames(PXD002882_hilm_real) %in% rownames(PXD002882_hilm_na),],
                              PXD002882_mihm_real[rownames(PXD002882_mihm_real) %in% rownames(PXD002882_mihm_na),],
                              PXD002882_mimm_real[rownames(PXD002882_mimm_real) %in% rownames(PXD002882_mimm_na),],
                              PXD002882_milm_real[rownames(PXD002882_milm_real) %in% rownames(PXD002882_milm_na),],
                              PXD002882_lihm_real[rownames(PXD002882_lihm_real) %in% rownames(PXD002882_lihm_na),],
                              PXD002882_limm_real[rownames(PXD002882_limm_real) %in% rownames(PXD002882_limm_na),],
                              PXD002882_lilm_real[rownames(PXD002882_lilm_real) %in% rownames(PXD002882_lilm_na),])
#bpca
PXD002882_whole_bpca <- pca(PXD002882_whole_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
PXD002882_whole_bpca <- as.data.frame(PXD002882_whole_bpca@completeObs)
#svd
PXD002882_whole_svd <- as.data.frame(impute.wrapper.SVD(PXD002882_whole_na, K = 10))
#knn
PXD002882_whole_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(PXD002882_whole_na), K = 5))
#mle
PXD002882_whole_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(PXD002882_whole_na)))
#lls
PXD002882_whole_lls <- llsImpute(PXD002882_whole_na, k = 5, center = TRUE, completeObs = TRUE)
PXD002882_whole_lls <- as.data.frame(PXD002882_whole_lls@completeObs)
#rf
PXD002882_whole_rf <- missForest(PXD002882_whole_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
PXD002882_whole_rf <- PXD002882_whole_rf$ximp

#for pimms
write.csv(t(PXD002882_whole_na), file = "PXD002882_whole_na.csv")
#import pimms result
#fread for fast read
PXD002882_whole_cf <- fread("PXD002882_whole_cf.csv", header = TRUE)
PXD002882_whole_dae <- fread("PXD002882_whole_dae.csv", header = TRUE)
PXD002882_whole_vae <- fread("PXD002882_whole_vae.csv", header = TRUE)

PXD002882_whole_cf <- as.data.frame(PXD002882_whole_cf)
PXD002882_whole_dae <- as.data.frame(PXD002882_whole_dae)
PXD002882_whole_vae <- as.data.frame(PXD002882_whole_vae)

#set colname
colnames(PXD002882_whole_dae) <- PXD002882_whole_dae[1,]
colnames(PXD002882_whole_vae) <- PXD002882_whole_vae[1,]

#remove excessive rows
PXD002882_whole_dae <- PXD002882_whole_dae[-c(1,2),]
PXD002882_whole_vae <- PXD002882_whole_vae[-c(1,2),]

#turn chr into num
PXD002882_whole_dae[,-1] <- sapply(PXD002882_whole_dae[,-1], as.numeric)
PXD002882_whole_vae[,-1] <- sapply(PXD002882_whole_vae[,-1], as.numeric)

#set sample id as row name
rownames(PXD002882_whole_cf) <- PXD002882_whole_cf[,1]
rownames(PXD002882_whole_dae) <- PXD002882_whole_dae[,1]
rownames(PXD002882_whole_vae) <- PXD002882_whole_vae[,1]

PXD002882_whole_cf <- as.data.frame(PXD002882_whole_cf[,-1])
PXD002882_whole_dae <- as.data.frame(PXD002882_whole_dae[,-1])
PXD002882_whole_vae <- as.data.frame(PXD002882_whole_vae[,-1])

#assemble the optimal method matrix
PXD002882_whole_mix <- rbind(PXD002882_hihm_rf, PXD002882_himm_rf, PXD002882_hilm_rf,
                             PXD002882_mihm_cf, PXD002882_mimm_rf, PXD002882_milm_rf,
                             PXD002882_lihm_rf, PXD002882_limm_rf)

#get na data
comp_na <- PXD002882_whole_na
#get real data
comp_real <- PXD002882_whole_real[which(rownames(PXD002882_whole_real) %in% rownames(comp_na)),]

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

methods <- c("bpca", "knn", "lls", "mle", "svd", "rf", "cf", "dae", "vae", "mix")
id <- "PXD002882"
#get data for different methods
for (j in 1:length(methods)) {
  if (paste(id, "whole", methods[j], sep = "_") %in% ls() == TRUE) {
    temp <- get(paste(id, "whole", methods[j], sep = "_"))
    for (k in 1:length(nrmse_sum$sample)) {
      nrmse_sum[k, grep(methods[j], colnames(nrmse_sum))] <- sqrt(mean((temp[,k] - comp_real[,k])^2, na.rm = TRUE))/sd(comp_real[,i], na.rm = TRUE)
    }
  } else {
    next
  }
  #if highest value > 4 fold of median, remove values higher than 5*(50% values) in a method
  if (max(nrmse_sum[,grep(methods[j], colnames(nrmse_sum))]) > 4*quantile(nrmse_sum[,grep(methods[j], colnames(nrmse_sum))], probs = c(0,0.5,0.9,0.95,1))[2]) {
    print(methods[j])
    print(quantile(nrmse_sum[,grep(methods[j], colnames(nrmse_sum))], probs = c(0,0.5,0.9,0.95,1)))
    nrmse_sum[which(nrmse_sum[,grep(methods[j], colnames(nrmse_sum))] > 3*quantile(nrmse_sum[,grep(methods[j], colnames(nrmse_sum))], probs = c(0,0.5,0.9,0.95,1))[2]),
              grep(methods[j], colnames(nrmse_sum))] <- NA
  }
}


#barplot
#calculate mean and standard error
#drop mle
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
  ylim(0,1.4) +
  scale_fill_brewer(labels = c("RF", "BPCA", "kNN", "LLS", "MLE", "SVD", "CF", "DAE", "VAE", "Mix"), palette = "Set3") +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 15, color = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.major = element_line(linewidth = 1),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15))

ggsave("figure7B_mix_compare_PXD002882.tiff",
       plot = last_plot(),
       path = "figures/manuscript/",
       device = "tiff",
       width = 2244,
       height = 1400,
       units = "px",
       dpi = 300)
ggsave("figure7B_mix_compare_PXD002882.jpeg",
       plot = last_plot(),
       path = "figures/manuscript/",
       device = "jpeg",
       width = 2800,
       height = 1400,
       units = "px",
       dpi = 300)


##differential expression analysis
PXD002882_meta <- data.frame("sample" = colnames(PXD002882_whole_real),
                             "dose" = str_split_i(colnames(PXD002882_whole_real), " ", 2),
                             "treatment" = str_split_i(colnames(PXD002882_whole_real), " ", 3),
                             "status" = NA)

for (i in 1:length(PXD002882_meta$treatment)) {
  if (is.na(PXD002882_meta$treatment[i] == TRUE)) {
    PXD002882_meta$treatment[i] <- PXD002882_meta$dose[i]
  }
}
for (i in 1:length(PXD002882_meta$treatment)) {
  if (PXD002882_meta$treatment[i] == PXD002882_meta$dose[i]) {
    PXD002882_meta$dose[i] <- "N"
  }
}
PXD002882_meta$treatment[c(1:3)] <- "CD_CoA_Severe_230"
PXD002882_meta$status <- str_split_i(PXD002882_meta$treatment, "_", 3)
PXD002882_meta$status[is.na(PXD002882_meta$status)] <- "Control"
PXD002882_meta$cell <- str_split_i(PXD002882_meta$treatment, "_", 4)
PXD002882_meta$cell[PXD002882_meta$status == "Control"] <- str_split_i(PXD002882_meta$treatment[PXD002882_meta$status == "Control"], "_", 2)
PXD002882_meta$treatment[1:63] <- "CD_CoA"
PXD002882_meta$treatment[64:93] <- "Control"


DE_result <- data.frame("df" = ls()[grep("PXD002882_whole_", ls())],
                        "method" = NA,
                        "sig_DEP" = NA,
                        "sig_DF" = NA)
DE_result$method <- str_split_i(DE_result$df, "_", 3)

#limma
design <- model.matrix(~ 0 + treatment + dose + status + cell, data = PXD002882_meta)

for (i in 1:length(DE_result$df)) {
  data <- get(DE_result$df[i])
  rid <- rownames(data)
  cid <- colnames(data)
  data <- normalize.quantiles(as.matrix(data))
  rownames(data) <- rid
  colnames(data) <- cid
  fit <- lmFit(data, design)
  con <- makeContrasts(treatmentCD_CoA - treatmentControl, levels = design)
  fit.con <- contrasts.fit(fit, con)
  bayes <- eBayes(fit.con)
  result <- topTable(bayes, sort.by="none", n=Inf, adjust.method = "BH")
  sig <- result[result$adj.P.Val < 0.05 & abs(result$logFC) > 1,]
  sig <- sig[!is.na(sig$logFC),]
  dfname <- paste(DE_result$df[i], "sigDEP", sep = "_")
  assign(dfname, sig)
  DE_result$sig_DEP[i] <- length(rownames(sig))
  DE_result$sig_DF[i] <- dfname
  print(paste(DE_result$df[i], length(sig$logFC), "sigDEP", sep = " "))
}

DE_result <- DE_result[DE_result$method != "na",]
DE_result$method <- factor(DE_result$method, levels = c("real", "rf", "bpca", "knn", "lls", "mle", 
                                                        "svd", "cf", "dae", "vae", "mix"))

ggplot(data = DE_result, aes(x = method, y = sig_DEP, fill = method)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sig_DEP, vjust = -0.5)) +
  scale_fill_discrete(labels = c("Real", "RF", "BPCA", "kNN", "LLS", "MLE", "SVD", "CF", "DAE", "VAE", "Mix")) +
  scale_x_discrete(labels = c("Real", "RF", "BPCA", "kNN", "LLS", "MLE", "SVD", "CF", "DAE", "VAE", "Mix")) +
  labs(x = "Method",
       y = "Number of significant differentially expressed peptides") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12, face = "bold"))

#reproducibility
auc_summary_plot <- readRDS("C:/Users/ys25.stu/OneDrive - UBC/BeeCSI/Proteomics/imputation/auc_summary_plot.rds")
auc_summary_plot$Method <- as.character(auc_summary_plot$Method)
#KNN to kNN
auc_summary_plot$Method[4] <- "kNN"
#MIX to Mix
auc_summary_plot$Method[5] <- "Mix"

#remove NA
auc_summary_plot <- auc_summary_plot[-7,]
auc_summary_plot$Method <- factor(auc_summary_plot$Method, levels = c("RF", "BPCA", "kNN", "LLS", "MLE", 
                                                               "SVD", "CF", "DAE", "VAE", "Mix"))


# accuracy barplot, marked numbers on the top
ggplot(auc_summary_plot, aes(x = Method, y = Accuracy, fill = Method)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = scales::percent(Accuracy, accuracy = 0.01)), vjust = -0.3, size = 5) +  # 在柱上标注准确率
  ylim(0,0.8) +
  theme_bw() +
  labs(title = "", x = "Method", y = "Reproducibility") +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 15, color = "black"),
        panel.grid = element_blank(),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15)) +
  scale_fill_brewer(palette = "Set3")
ggsave("figure8_reproducibility_PXD002882.jpeg",
       plot = last_plot(),
       path = "figures/manuscript/",
       device = "jpeg",
       width = 2400,
       height = 1600,
       units = "px",
       dpi = 300)




#PXD047528
PXD047528_dda <- read_delim("reference/Pietz.2024/data/datasets/blood_hiv_dda_dia/dda/peptides.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
PXD047528_dia <- read_delim("reference/Pietz.2024/data/datasets/blood_hiv_dda_dia/dia/peptides.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(PXD047528_dia)
colnames(PXD047528_dda)
length(which(PXD047528_dda$Sequence %in% PXD047528_dia$Sequence))
length(which(is.na(PXD047528_dda[,-c(1:3)]) == TRUE))/(nrow(PXD047528_dda)*(ncol(PXD047528_dda)-3))
length(which(is.na(PXD047528_dia[,-c(1,2,18)]) == TRUE))/(nrow(PXD047528_dda)*(ncol(PXD047528_dda)-3))

PXD047528_peptide <- as.data.frame(PXD047528_dda[which(PXD047528_dda$Sequence %in% PXD047528_dia$Sequence),c(1,4:18)])
#set row name
rownames(PXD047528_peptide) <- PXD047528_peptide[,1]
PXD047528_peptide[,-1] <- log2(PXD047528_peptide[,-1])
PXD047528_peptide <- PXD047528_peptide[,-1]

#calculate NA rate, mean intensity
PXD047528_peptide <- add_column(PXD047528_peptide, mean_intensity = NA, .before = "10_Slot1")
PXD047528_peptide <- add_column(PXD047528_peptide, NA_rate = NA, .after = "mean_intensity")

PXD047528_peptide$mean_intensity <- rowMeans(PXD047528_peptide[, -c(1,2)], na.rm = TRUE)
PXD047528_peptide$NA_rate <- rowSums(is.na(PXD047528_peptide[,-c(1,2)]))/(ncol(PXD047528_peptide) - 2)

PXD047528_int_low <- quantile(PXD047528_peptide$mean_intensity, na.rm = TRUE)[2]
PXD047528_int_high <- quantile(PXD047528_peptide$mean_intensity, na.rm = TRUE)[4]

#NA rate vs intensity
ggplot(data = PXD047528_peptide, aes(x = mean_intensity, y = NA_rate)) + 
  geom_point(alpha = 0.1, shape = 20) +
  labs(title = paste("PXD047528 Crohn's disease", "Relationship between peptide mean intensity and missing rate", sep = "\n"),
       x = bquote(paste("mean log"["2"]*"(intensity)")),
       y = "Missing rate") +
  stat_smooth(method = "lm", geom = "smooth", formula = y~x, se = TRUE, color = "black", size = 2) +
  lims(x = c(12, 23),
       y = c(0,1)) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        title = element_text(size = 20)) 
last_plot() +
  stat_density2d(aes(fill = ..density..), geom = "tile", contour = FALSE, n=200, alpha = 0.8) +
  scale_fill_continuous(low = "white", high = "dodgerblue4") +
  geom_line(aes(x = PXD047528_int_high), color = "red", size = 1) +
  geom_line(aes(x = PXD047528_int_low), color = "red", size = 1) +
  geom_line(aes(y = 0.75), color = "red", size = 1) +
  geom_line(aes(y = 0.25), color = "red", size = 1) 

#high intensity, high missing peptides
PXD047528_hihm <- PXD047528_peptide[which(PXD047528_peptide$mean_intensity >= PXD047528_int_high & 
                                            PXD047528_peptide$NA_rate >= 0.75 &
                                            PXD047528_peptide$NA_rate < 1),]
#high intensity, medium missing peptides
PXD047528_himm <- PXD047528_peptide[which(PXD047528_peptide$mean_intensity >= PXD047528_int_high & 
                                            PXD047528_peptide$NA_rate >= 0.25 &
                                            PXD047528_peptide$NA_rate < 0.75),]
#high intensity, low missing peptides
PXD047528_hilm <- PXD047528_peptide[which(PXD047528_peptide$mean_intensity >= PXD047528_int_high & 
                                            PXD047528_peptide$NA_rate < 0.25),]
#medium intensity, high missing peptides
PXD047528_mihm <- PXD047528_peptide[which(PXD047528_peptide$mean_intensity < PXD047528_int_high & 
                                            PXD047528_peptide$mean_intensity >= PXD047528_int_low &
                                            PXD047528_peptide$NA_rate >= 0.75 &
                                            PXD047528_peptide$NA_rate < 1),]
#medium intensity, medium missing peptides
PXD047528_mimm <- PXD047528_peptide[which(PXD047528_peptide$mean_intensity < PXD047528_int_high & 
                                            PXD047528_peptide$mean_intensity >= PXD047528_int_low &
                                            PXD047528_peptide$NA_rate >= 0.25 &
                                            PXD047528_peptide$NA_rate < 0.75),]
#medium intensity, low missing peptides
PXD047528_milm <- PXD047528_peptide[which(PXD047528_peptide$mean_intensity < PXD047528_int_high & 
                                            PXD047528_peptide$mean_intensity >= PXD047528_int_low &
                                            PXD047528_peptide$NA_rate < 0.25),]
#low intensity, high missing peptides
PXD047528_lihm <- PXD047528_peptide[which(PXD047528_peptide$mean_intensity < PXD047528_int_low & 
                                            PXD047528_peptide$NA_rate >= 0.75 &
                                            PXD047528_peptide$NA_rate < 1),]
#low intensity, medium missing peptides
PXD047528_limm <- PXD047528_peptide[which(PXD047528_peptide$mean_intensity < PXD047528_int_low & 
                                            PXD047528_peptide$NA_rate >= 0.25 &
                                            PXD047528_peptide$NA_rate < 0.75),]
#low intensity, low missing peptides
PXD047528_lilm <- PXD047528_peptide[which(PXD047528_peptide$mean_intensity < PXD047528_int_low & 
                                            PXD047528_peptide$NA_rate < 0.25),]

#keep intensity
PXD047528_hihm <- PXD047528_hihm[,-c(1,2)]
PXD047528_himm <- PXD047528_himm[,-c(1,2)]
PXD047528_hilm <- PXD047528_hilm[,-c(1,2)]
PXD047528_mihm <- PXD047528_mihm[,-c(1,2)]
PXD047528_mimm <- PXD047528_mimm[,-c(1,2)]
PXD047528_milm <- PXD047528_milm[,-c(1,2)]
PXD047528_lihm <- PXD047528_lihm[,-c(1,2)]
PXD047528_limm <- PXD047528_limm[,-c(1,2)]
PXD047528_lilm <- PXD047528_lilm[,-c(1,2)]

#prepare datasets for imputation, remove sample and pepetide that are all NA
PXD047528_hihm_na <- PXD047528_hihm[which(rowSums(is.na(PXD047528_hihm)) < ncol(PXD047528_hihm)),
                                    which(colSums(is.na(PXD047528_hihm)) < nrow(PXD047528_hihm))]
PXD047528_himm_na <- PXD047528_himm[which(rowSums(is.na(PXD047528_himm)) < ncol(PXD047528_himm)),
                                    which(colSums(is.na(PXD047528_himm)) < nrow(PXD047528_himm))]
PXD047528_hilm_na <- PXD047528_hilm[which(rowSums(is.na(PXD047528_hilm)) < ncol(PXD047528_hilm)),
                                    which(colSums(is.na(PXD047528_hilm)) < nrow(PXD047528_hilm))]
PXD047528_mihm_na <- PXD047528_mihm[which(rowSums(is.na(PXD047528_mihm)) < ncol(PXD047528_mihm)),
                                    which(colSums(is.na(PXD047528_mihm)) < nrow(PXD047528_mihm))]
PXD047528_mimm_na <- PXD047528_mimm[which(rowSums(is.na(PXD047528_mimm)) < ncol(PXD047528_mimm)),
                                    which(colSums(is.na(PXD047528_mimm)) < nrow(PXD047528_mimm))]
PXD047528_milm_na <- PXD047528_milm[which(rowSums(is.na(PXD047528_milm)) < ncol(PXD047528_milm)),
                                    which(colSums(is.na(PXD047528_milm)) < nrow(PXD047528_milm))]
PXD047528_lihm_na <- PXD047528_lihm[which(rowSums(is.na(PXD047528_lihm)) < ncol(PXD047528_lihm)),
                                    which(colSums(is.na(PXD047528_lihm)) < nrow(PXD047528_lihm))]
PXD047528_limm_na <- PXD047528_limm[which(rowSums(is.na(PXD047528_limm)) < ncol(PXD047528_limm)),
                                    which(colSums(is.na(PXD047528_limm)) < nrow(PXD047528_limm))]
PXD047528_lilm_na <- PXD047528_lilm[which(rowSums(is.na(PXD047528_lilm)) < ncol(PXD047528_lilm)),
                                    which(colSums(is.na(PXD047528_lilm)) < nrow(PXD047528_lilm))]

#real datasets
PXD047528_hihm_real <- PXD047528_hihm_na
PXD047528_himm_real <- PXD047528_himm_na
PXD047528_hilm_real <- PXD047528_hilm_na
PXD047528_mihm_real <- PXD047528_mihm_na
PXD047528_mimm_real <- PXD047528_mimm_na
PXD047528_milm_real <- PXD047528_milm_na
PXD047528_lihm_real <- PXD047528_lihm_na
PXD047528_limm_real <- PXD047528_limm_na
PXD047528_lilm_real <- PXD047528_lilm_na

#generate random sample points, 20%
which(is.na(PXD047528_hihm_na) == FALSE)

na_pos_PXD047528_hihm <- sample(which(is.na(PXD047528_hihm_na) == FALSE), 
                                0.2*length(which(is.na(PXD047528_hihm_na) == FALSE)))
na_pos_PXD047528_himm <- sample(which(is.na(PXD047528_himm_na) == FALSE), 
                                0.2*length(which(is.na(PXD047528_himm_na) == FALSE)))
na_pos_PXD047528_hilm <- sample(which(is.na(PXD047528_hilm_na) == FALSE), 
                                0.2*length(which(is.na(PXD047528_hilm_na) == FALSE)))
na_pos_PXD047528_mihm <- sample(which(is.na(PXD047528_mihm_na) == FALSE), 
                                0.2*length(which(is.na(PXD047528_mihm_na) == FALSE)))
na_pos_PXD047528_mimm <- sample(which(is.na(PXD047528_mimm_na) == FALSE), 
                                0.2*length(which(is.na(PXD047528_mimm_na) == FALSE)))
na_pos_PXD047528_milm <- sample(which(is.na(PXD047528_milm_na) == FALSE), 
                                0.2*length(which(is.na(PXD047528_milm_na) == FALSE)))
na_pos_PXD047528_lihm <- sample(which(is.na(PXD047528_lihm_na) == FALSE), 
                                0.2*length(which(is.na(PXD047528_lihm_na) == FALSE)))
na_pos_PXD047528_limm <- sample(which(is.na(PXD047528_limm_na) == FALSE), 
                                0.2*length(which(is.na(PXD047528_limm_na) == FALSE)))
na_pos_PXD047528_lilm <- sample(which(is.na(PXD047528_lilm_na) == FALSE), 
                                0.2*length(which(is.na(PXD047528_lilm_na) == FALSE)))

#replace true values with NA
PXD047528_hihm_na <- as.matrix(PXD047528_hihm_na)
PXD047528_hihm_na[na_pos_PXD047528_hihm] <- NA
PXD047528_hihm_na <- as.data.frame(PXD047528_hihm_na)
PXD047528_himm_na <- as.matrix(PXD047528_himm_na)
PXD047528_himm_na[na_pos_PXD047528_himm] <- NA
PXD047528_himm_na <- as.data.frame(PXD047528_himm_na)
PXD047528_hilm_na <- as.matrix(PXD047528_hilm_na)
PXD047528_hilm_na[na_pos_PXD047528_hilm] <- NA
PXD047528_hilm_na <- as.data.frame(PXD047528_hilm_na)
PXD047528_mihm_na <- as.matrix(PXD047528_mihm_na)
PXD047528_mihm_na[na_pos_PXD047528_mihm] <- NA
PXD047528_mihm_na <- as.data.frame(PXD047528_mihm_na)
PXD047528_mimm_na <- as.matrix(PXD047528_mimm_na)
PXD047528_mimm_na[na_pos_PXD047528_mimm] <- NA
PXD047528_mimm_na <- as.data.frame(PXD047528_mimm_na)
PXD047528_milm_na <- as.matrix(PXD047528_milm_na)
PXD047528_milm_na[na_pos_PXD047528_milm] <- NA
PXD047528_milm_na <- as.data.frame(PXD047528_milm_na)
PXD047528_lihm_na <- as.matrix(PXD047528_lihm_na)
PXD047528_lihm_na[na_pos_PXD047528_lihm] <- NA
PXD047528_lihm_na <- as.data.frame(PXD047528_lihm_na)
PXD047528_limm_na <- as.matrix(PXD047528_limm_na)
PXD047528_limm_na[na_pos_PXD047528_limm] <- NA
PXD047528_limm_na <- as.data.frame(PXD047528_limm_na)
PXD047528_lilm_na <- as.matrix(PXD047528_lilm_na)
PXD047528_lilm_na[na_pos_PXD047528_lilm] <- NA
PXD047528_lilm_na <- as.data.frame(PXD047528_lilm_na)

#change na pos: hihm
#count na percentage
PXD047528_hihm_na <- PXD047528_hihm_na[-which(rowSums(is.na(PXD047528_hihm_na)) == ncol(PXD047528_hihm_na)),]
PXD047528_mihm_na <- PXD047528_mihm_na[-which(rowSums(is.na(PXD047528_mihm_na)) == ncol(PXD047528_mihm_na)),]
PXD047528_lihm_na <- PXD047528_lihm_na[-which(rowSums(is.na(PXD047528_lihm_na)) == ncol(PXD047528_lihm_na)),]

#impute by intensity and missing rate
#bpca, svd, knn (SEQKNN, TRKNN, KNNIMPUTE), median, missForrest, RSN
#bpca
PXD047528_hihm_bpca <- pca(PXD047528_hihm_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
PXD047528_himm_bpca <- pca(PXD047528_himm_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
PXD047528_hilm_bpca <- pca(PXD047528_hilm_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
PXD047528_mihm_bpca <- pca(PXD047528_mihm_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
PXD047528_mimm_bpca <- pca(PXD047528_mimm_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
PXD047528_milm_bpca <- pca(PXD047528_milm_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
PXD047528_lihm_bpca <- pca(PXD047528_lihm_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
PXD047528_limm_bpca <- pca(PXD047528_limm_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
PXD047528_lilm_bpca <- pca(PXD047528_lilm_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
PXD047528_hihm_bpca <- as.data.frame(PXD047528_hihm_bpca@completeObs)
PXD047528_himm_bpca <- as.data.frame(PXD047528_himm_bpca@completeObs)
PXD047528_hilm_bpca <- as.data.frame(PXD047528_hilm_bpca@completeObs)
PXD047528_mihm_bpca <- as.data.frame(PXD047528_mihm_bpca@completeObs)
PXD047528_mimm_bpca <- as.data.frame(PXD047528_mimm_bpca@completeObs)
PXD047528_milm_bpca <- as.data.frame(PXD047528_milm_bpca@completeObs)
PXD047528_lihm_bpca <- as.data.frame(PXD047528_lihm_bpca@completeObs)
PXD047528_limm_bpca <- as.data.frame(PXD047528_limm_bpca@completeObs)
PXD047528_lilm_bpca <- as.data.frame(PXD047528_lilm_bpca@completeObs)

#svd
PXD047528_hihm_svd <- as.data.frame(impute.wrapper.SVD(PXD047528_hihm_na, K = 10))
PXD047528_himm_svd <- as.data.frame(impute.wrapper.SVD(PXD047528_himm_na, K = 10))
PXD047528_hilm_svd <- as.data.frame(impute.wrapper.SVD(PXD047528_hilm_na, K = 10))
PXD047528_mihm_svd <- as.data.frame(impute.wrapper.SVD(PXD047528_mihm_na, K = 10))
PXD047528_mimm_svd <- as.data.frame(impute.wrapper.SVD(PXD047528_mimm_na, K = 10))
PXD047528_milm_svd <- as.data.frame(impute.wrapper.SVD(PXD047528_milm_na, K = 10))
PXD047528_lihm_svd <- as.data.frame(impute.wrapper.SVD(PXD047528_lihm_na, K = 10))
PXD047528_limm_svd <- as.data.frame(impute.wrapper.SVD(PXD047528_limm_na, K = 10))
PXD047528_lilm_svd <- as.data.frame(impute.wrapper.SVD(PXD047528_lilm_na, K = 10))

#knn
PXD047528_hihm_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(PXD047528_hihm_na), K = 5))
PXD047528_himm_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(PXD047528_himm_na), K = 5))
PXD047528_hilm_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(PXD047528_hilm_na), K = 5))
PXD047528_mihm_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(PXD047528_mihm_na), K = 5))
PXD047528_mimm_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(PXD047528_mimm_na), K = 5))
PXD047528_milm_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(PXD047528_milm_na), K = 5))
PXD047528_lihm_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(PXD047528_lihm_na), K = 5))
PXD047528_limm_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(PXD047528_limm_na), K = 5))
PXD047528_lilm_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(PXD047528_lilm_na), K = 5))

#mle
PXD047528_hihm_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(PXD047528_hihm_na)))
PXD047528_himm_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(PXD047528_himm_na)))
PXD047528_hilm_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(PXD047528_hilm_na)))
PXD047528_mihm_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(PXD047528_mihm_na)))
PXD047528_mimm_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(PXD047528_mimm_na)))
PXD047528_milm_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(PXD047528_milm_na)))
PXD047528_lihm_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(PXD047528_lihm_na)))
PXD047528_limm_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(PXD047528_limm_na)))
PXD047528_lilm_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(PXD047528_lilm_na)))

#lls
PXD047528_hihm_lls <- llsImpute(PXD047528_hihm_na, k = 5, center = TRUE, completeObs = TRUE)
PXD047528_himm_lls <- llsImpute(PXD047528_himm_na, k = 5, center = TRUE, completeObs = TRUE)
PXD047528_hilm_lls <- llsImpute(PXD047528_hilm_na, k = 5, center = TRUE, completeObs = TRUE)
PXD047528_mihm_lls <- llsImpute(PXD047528_mihm_na, k = 5, center = TRUE, completeObs = TRUE)
PXD047528_mimm_lls <- llsImpute(PXD047528_mimm_na, k = 5, center = TRUE, completeObs = TRUE)
PXD047528_milm_lls <- llsImpute(PXD047528_milm_na, k = 5, center = TRUE, completeObs = TRUE)
PXD047528_lihm_lls <- llsImpute(PXD047528_lihm_na, k = 5, center = TRUE, completeObs = TRUE)
PXD047528_limm_lls <- llsImpute(PXD047528_limm_na, k = 5, center = TRUE, completeObs = TRUE)
PXD047528_lilm_lls <- llsImpute(PXD047528_lilm_na, k = 5, center = TRUE, completeObs = TRUE)
PXD047528_hihm_lls <- as.data.frame(PXD047528_hihm_lls@completeObs)
PXD047528_himm_lls <- as.data.frame(PXD047528_himm_lls@completeObs)
PXD047528_hilm_lls <- as.data.frame(PXD047528_hilm_lls@completeObs)
PXD047528_mihm_lls <- as.data.frame(PXD047528_mihm_lls@completeObs)
PXD047528_mimm_lls <- as.data.frame(PXD047528_mimm_lls@completeObs)
PXD047528_milm_lls <- as.data.frame(PXD047528_milm_lls@completeObs)
PXD047528_lihm_lls <- as.data.frame(PXD047528_lihm_lls@completeObs)
PXD047528_limm_lls <- as.data.frame(PXD047528_limm_lls@completeObs)
PXD047528_lilm_lls <- as.data.frame(PXD047528_lilm_lls@completeObs)

#rf
#register number of cores for parallel
registerDoParallel(cores = 15)
PXD047528_hihm_rf <- missForest(PXD047528_hihm_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
PXD047528_himm_rf <- missForest(PXD047528_himm_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
PXD047528_hilm_rf <- missForest(PXD047528_hilm_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
PXD047528_mihm_rf <- missForest(PXD047528_mihm_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
PXD047528_mimm_rf <- missForest(PXD047528_mimm_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
PXD047528_milm_rf <- missForest(PXD047528_milm_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
PXD047528_lihm_rf <- missForest(PXD047528_lihm_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
PXD047528_limm_rf <- missForest(PXD047528_limm_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
PXD047528_lilm_rf <- missForest(PXD047528_lilm_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
PXD047528_hihm_rf <- PXD047528_hihm_rf$ximp
PXD047528_himm_rf <- PXD047528_himm_rf$ximp
PXD047528_hilm_rf <- PXD047528_hilm_rf$ximp
PXD047528_mihm_rf <- PXD047528_mihm_rf$ximp
PXD047528_mimm_rf <- PXD047528_mimm_rf$ximp
PXD047528_milm_rf <- PXD047528_milm_rf$ximp
PXD047528_lihm_rf <- PXD047528_lihm_rf$ximp
PXD047528_limm_rf <- PXD047528_limm_rf$ximp
PXD047528_lilm_rf <- PXD047528_lilm_rf$ximp



#PXD000279
PXD000279 <- read_delim("reference/Pietz.2024/data/datasets/PXD000279_maxlfq_benchmark_human_ecoli_mixture/peptides.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
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

PXD000279_hihm_na <- PXD000279_hihm_na[-which(rowSums(is.na(PXD000279_hihm_na)) == ncol(PXD000279_hihm_na)),]
PXD000279_mihm_na <- PXD000279_mihm_na[-which(rowSums(is.na(PXD000279_mihm_na)) == ncol(PXD000279_mihm_na)),]
PXD000279_lihm_na <- PXD000279_lihm_na[-which(rowSums(is.na(PXD000279_lihm_na)) == ncol(PXD000279_lihm_na)),]
PXD000279_limm_na <- PXD000279_limm_na[-which(rowSums(is.na(PXD000279_limm_na)) == ncol(PXD000279_limm_na)),]

comp_real <- PXD000279_mihm_real[which(rownames(PXD000279_mihm_real) %in% rownames(PXD000279_mihm_na)),] %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "real")
comp_na <- PXD000279_mihm_na %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "na")

comp <- data.frame(sample = comp_real$sample,
                   real = comp_real$real,
                   na = comp_na$na)


#impute by intensity and missing rate
#bpca, svd, knn (SEQKNN, TRKNN, KNNIMPUTE), median, missForrest, RSN
#bpca
PXD000279_hihm_bpca <- pca(PXD000279_hihm_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
PXD000279_himm_bpca <- pca(PXD000279_himm_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
PXD000279_hilm_bpca <- pca(PXD000279_hilm_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
PXD000279_mihm_bpca <- pca(PXD000279_mihm_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
PXD000279_mimm_bpca <- pca(PXD000279_mimm_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
PXD000279_milm_bpca <- pca(PXD000279_milm_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
PXD000279_lihm_bpca <- pca(PXD000279_lihm_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
PXD000279_limm_bpca <- pca(PXD000279_limm_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
PXD000279_lilm_bpca <- pca(PXD000279_lilm_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
PXD000279_hihm_bpca <- as.data.frame(PXD000279_hihm_bpca@completeObs)
PXD000279_himm_bpca <- as.data.frame(PXD000279_himm_bpca@completeObs)
PXD000279_hilm_bpca <- as.data.frame(PXD000279_hilm_bpca@completeObs)
PXD000279_mihm_bpca <- as.data.frame(PXD000279_mihm_bpca@completeObs)
PXD000279_mimm_bpca <- as.data.frame(PXD000279_mimm_bpca@completeObs)
PXD000279_milm_bpca <- as.data.frame(PXD000279_milm_bpca@completeObs)
PXD000279_lihm_bpca <- as.data.frame(PXD000279_lihm_bpca@completeObs)
PXD000279_limm_bpca <- as.data.frame(PXD000279_limm_bpca@completeObs)
PXD000279_lilm_bpca <- as.data.frame(PXD000279_lilm_bpca@completeObs)

#svd
PXD000279_hihm_svd <- as.data.frame(impute.wrapper.SVD(PXD000279_hihm_na, K = 10))
PXD000279_himm_svd <- as.data.frame(impute.wrapper.SVD(PXD000279_himm_na, K = 10))
PXD000279_hilm_svd <- as.data.frame(impute.wrapper.SVD(PXD000279_hilm_na, K = 10))
PXD000279_mihm_svd <- as.data.frame(impute.wrapper.SVD(PXD000279_mihm_na, K = 10))
PXD000279_mimm_svd <- as.data.frame(impute.wrapper.SVD(PXD000279_mimm_na, K = 10))
PXD000279_milm_svd <- as.data.frame(impute.wrapper.SVD(PXD000279_milm_na, K = 10))
PXD000279_lihm_svd <- as.data.frame(impute.wrapper.SVD(PXD000279_lihm_na, K = 10))
PXD000279_limm_svd <- as.data.frame(impute.wrapper.SVD(PXD000279_limm_na, K = 10))
PXD000279_lilm_svd <- as.data.frame(impute.wrapper.SVD(PXD000279_lilm_na, K = 10))

#knn
PXD000279_hihm_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(PXD000279_hihm_na), K = 5))
PXD000279_himm_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(PXD000279_himm_na), K = 5))
PXD000279_hilm_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(PXD000279_hilm_na), K = 5))
PXD000279_mihm_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(PXD000279_mihm_na), K = 5))
PXD000279_mimm_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(PXD000279_mimm_na), K = 5))
PXD000279_milm_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(PXD000279_milm_na), K = 5))
PXD000279_lihm_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(PXD000279_lihm_na), K = 5))
PXD000279_limm_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(PXD000279_limm_na), K = 5))
PXD000279_lilm_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(PXD000279_lilm_na), K = 5))

#mle
PXD000279_hihm_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(PXD000279_hihm_na)))
PXD000279_himm_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(PXD000279_himm_na)))
PXD000279_hilm_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(PXD000279_hilm_na)))
PXD000279_mihm_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(PXD000279_mihm_na)))
PXD000279_mimm_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(PXD000279_mimm_na)))
PXD000279_milm_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(PXD000279_milm_na)))
PXD000279_lihm_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(PXD000279_lihm_na)))
PXD000279_limm_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(PXD000279_limm_na)))
PXD000279_lilm_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(PXD000279_lilm_na)))

#lls
PXD000279_hihm_lls <- llsImpute(PXD000279_hihm_na, k = 5, center = TRUE, completeObs = TRUE)
PXD000279_himm_lls <- llsImpute(PXD000279_himm_na, k = 5, center = TRUE, completeObs = TRUE)
PXD000279_hilm_lls <- llsImpute(PXD000279_hilm_na, k = 5, center = TRUE, completeObs = TRUE)
PXD000279_mihm_lls <- llsImpute(PXD000279_mihm_na, k = 5, center = TRUE, completeObs = TRUE)
PXD000279_mimm_lls <- llsImpute(PXD000279_mimm_na, k = 5, center = TRUE, completeObs = TRUE)
PXD000279_milm_lls <- llsImpute(PXD000279_milm_na, k = 5, center = TRUE, completeObs = TRUE)
PXD000279_lihm_lls <- llsImpute(PXD000279_lihm_na, k = 5, center = TRUE, completeObs = TRUE)
PXD000279_limm_lls <- llsImpute(PXD000279_limm_na, k = 5, center = TRUE, completeObs = TRUE)
PXD000279_lilm_lls <- llsImpute(PXD000279_lilm_na, k = 5, center = TRUE, completeObs = TRUE)
PXD000279_hihm_lls <- as.data.frame(PXD000279_hihm_lls@completeObs)
PXD000279_himm_lls <- as.data.frame(PXD000279_himm_lls@completeObs)
PXD000279_hilm_lls <- as.data.frame(PXD000279_hilm_lls@completeObs)
PXD000279_mihm_lls <- as.data.frame(PXD000279_mihm_lls@completeObs)
PXD000279_mimm_lls <- as.data.frame(PXD000279_mimm_lls@completeObs)
PXD000279_milm_lls <- as.data.frame(PXD000279_milm_lls@completeObs)
PXD000279_lihm_lls <- as.data.frame(PXD000279_lihm_lls@completeObs)
PXD000279_limm_lls <- as.data.frame(PXD000279_limm_lls@completeObs)
PXD000279_lilm_lls <- as.data.frame(PXD000279_lilm_lls@completeObs)

#rf
#register number of cores for parallel
registerDoParallel(cores = 6)
PXD000279_hihm_rf <- missForest(PXD000279_hihm_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
PXD000279_himm_rf <- missForest(PXD000279_himm_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
PXD000279_hilm_rf <- missForest(PXD000279_hilm_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
PXD000279_mihm_rf <- missForest(PXD000279_mihm_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
PXD000279_mimm_rf <- missForest(PXD000279_mimm_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
PXD000279_milm_rf <- missForest(PXD000279_milm_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
PXD000279_lihm_rf <- missForest(PXD000279_lihm_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
PXD000279_limm_rf <- missForest(PXD000279_limm_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
PXD000279_lilm_rf <- missForest(PXD000279_lilm_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
PXD000279_hihm_rf <- PXD000279_hihm_rf$ximp
PXD000279_himm_rf <- PXD000279_himm_rf$ximp
PXD000279_hilm_rf <- PXD000279_hilm_rf$ximp
PXD000279_mihm_rf <- PXD000279_mihm_rf$ximp
PXD000279_mimm_rf <- PXD000279_mimm_rf$ximp
PXD000279_milm_rf <- PXD000279_milm_rf$ximp
PXD000279_lihm_rf <- PXD000279_lihm_rf$ximp
PXD000279_limm_rf <- PXD000279_limm_rf$ximp
PXD000279_lilm_rf <- PXD000279_lilm_rf$ximp

#transpose data frame for pimms
write.csv(t(PXD000279_hihm_na), file = "PXD000279_hihm_na.csv")
write.csv(t(PXD000279_himm_na), file = "PXD000279_himm_na.csv")
write.csv(t(PXD000279_hilm_na), file = "PXD000279_hilm_na.csv")
write.csv(t(PXD000279_mihm_na), file = "PXD000279_mihm_na.csv")
write.csv(t(PXD000279_mimm_na), file = "PXD000279_mimm_na.csv")
write.csv(t(PXD000279_milm_na), file = "PXD000279_milm_na.csv")
write.csv(t(PXD000279_lihm_na), file = "PXD000279_lihm_na.csv")
write.csv(t(PXD000279_limm_na), file = "PXD000279_limm_na.csv")
write.csv(t(PXD000279_lilm_na), file = "PXD000279_lilm_na.csv")

#import pimms result
#cf
PXD000279_hilm_cf <- read.csv("PXD000279_hilm_cf.csv")
PXD000279_milm_cf <- read.csv("PXD000279_milm_cf.csv")
PXD000279_limm_cf <- read.csv("PXD000279_limm_cf.csv")
PXD000279_lilm_cf <- read.csv("PXD000279_lilm_cf.csv")

#vae
PXD000279_hihm_vae <- read.csv("PXD000279_hihm_vae.csv")
PXD000279_himm_vae <- read.csv("PXD000279_himm_vae.csv")
PXD000279_hilm_vae <- read.csv("PXD000279_hilm_vae.csv")
PXD000279_mihm_vae <- read.csv("PXD000279_mihm_vae.csv")
PXD000279_mimm_vae <- read.csv("PXD000279_mimm_vae.csv")
PXD000279_milm_vae <- read.csv("PXD000279_milm_vae.csv")
PXD000279_lihm_vae <- read.csv("PXD000279_lihm_vae.csv")
PXD000279_limm_vae <- read.csv("PXD000279_limm_vae.csv")
PXD000279_lilm_vae <- read.csv("PXD000279_lilm_vae.csv")

#dae
PXD000279_hihm_dae <- read.csv("PXD000279_hihm_dae.csv")
PXD000279_himm_dae <- read.csv("PXD000279_himm_dae.csv")
PXD000279_hilm_dae <- read.csv("PXD000279_hilm_dae.csv")
PXD000279_mihm_dae <- read.csv("PXD000279_mihm_dae.csv")
PXD000279_mimm_dae <- read.csv("PXD000279_mimm_dae.csv")
PXD000279_milm_dae <- read.csv("PXD000279_milm_dae.csv")
PXD000279_lihm_dae <- read.csv("PXD000279_lihm_dae.csv")
PXD000279_limm_dae <- read.csv("PXD000279_limm_dae.csv")
PXD000279_lilm_dae <- read.csv("PXD000279_lilm_dae.csv")

#set colname
colnames(PXD000279_hihm_dae) <- PXD000279_hihm_dae[1,]
colnames(PXD000279_himm_dae) <- PXD000279_himm_dae[1,]
colnames(PXD000279_hilm_dae) <- PXD000279_hilm_dae[1,]
colnames(PXD000279_mihm_dae) <- PXD000279_mihm_dae[1,]
colnames(PXD000279_mimm_dae) <- PXD000279_mimm_dae[1,]
colnames(PXD000279_milm_dae) <- PXD000279_milm_dae[1,]
colnames(PXD000279_lihm_dae) <- PXD000279_lihm_dae[1,]
colnames(PXD000279_limm_dae) <- PXD000279_limm_dae[1,]
colnames(PXD000279_lilm_dae) <- PXD000279_lilm_dae[1,]

colnames(PXD000279_hihm_vae) <- PXD000279_hihm_vae[1,]
colnames(PXD000279_himm_vae) <- PXD000279_himm_vae[1,]
colnames(PXD000279_hilm_vae) <- PXD000279_hilm_vae[1,]
colnames(PXD000279_mihm_vae) <- PXD000279_mihm_vae[1,]
colnames(PXD000279_mimm_vae) <- PXD000279_mimm_vae[1,]
colnames(PXD000279_milm_vae) <- PXD000279_milm_vae[1,]
colnames(PXD000279_lihm_vae) <- PXD000279_lihm_vae[1,]
colnames(PXD000279_limm_vae) <- PXD000279_limm_vae[1,]
colnames(PXD000279_lilm_vae) <- PXD000279_lilm_vae[1,]

#remove excessive rows
PXD000279_hihm_dae <- PXD000279_hihm_dae[-c(1,2),]
PXD000279_himm_dae <- PXD000279_himm_dae[-c(1,2),]
PXD000279_hilm_dae <- PXD000279_hilm_dae[-c(1,2),]
PXD000279_mihm_dae <- PXD000279_mihm_dae[-c(1,2),]
PXD000279_mimm_dae <- PXD000279_mimm_dae[-c(1,2),]
PXD000279_milm_dae <- PXD000279_milm_dae[-c(1,2),]
PXD000279_lihm_dae <- PXD000279_lihm_dae[-c(1,2),]
PXD000279_limm_dae <- PXD000279_limm_dae[-c(1,2),]
PXD000279_lilm_dae <- PXD000279_lilm_dae[-c(1,2),]

PXD000279_hihm_vae <- PXD000279_hihm_vae[-c(1,2),]
PXD000279_himm_vae <- PXD000279_himm_vae[-c(1,2),]
PXD000279_hilm_vae <- PXD000279_hilm_vae[-c(1,2),]
PXD000279_mihm_vae <- PXD000279_mihm_vae[-c(1,2),]
PXD000279_mimm_vae <- PXD000279_mimm_vae[-c(1,2),]
PXD000279_milm_vae <- PXD000279_milm_vae[-c(1,2),]
PXD000279_lihm_vae <- PXD000279_lihm_vae[-c(1,2),]
PXD000279_limm_vae <- PXD000279_limm_vae[-c(1,2),]
PXD000279_lilm_vae <- PXD000279_lilm_vae[-c(1,2),]

#turn chr into num
PXD000279_hihm_dae[,-1] <- sapply(PXD000279_hihm_dae[,-1], as.numeric)
PXD000279_himm_dae[,-1] <- sapply(PXD000279_himm_dae[,-1], as.numeric)
PXD000279_hilm_dae[,-1] <- sapply(PXD000279_hilm_dae[,-1], as.numeric)
PXD000279_mihm_dae[,-1] <- sapply(PXD000279_mihm_dae[,-1], as.numeric)
PXD000279_mimm_dae[,-1] <- sapply(PXD000279_mimm_dae[,-1], as.numeric)
PXD000279_milm_dae[,-1] <- sapply(PXD000279_milm_dae[,-1], as.numeric)
PXD000279_lihm_dae[,-1] <- sapply(PXD000279_lihm_dae[,-1], as.numeric)
PXD000279_limm_dae[,-1] <- sapply(PXD000279_limm_dae[,-1], as.numeric)
PXD000279_lilm_dae[,-1] <- sapply(PXD000279_lilm_dae[,-1], as.numeric)

PXD000279_hihm_vae[,-1] <- sapply(PXD000279_hihm_vae[,-1], as.numeric)
PXD000279_himm_vae[,-1] <- sapply(PXD000279_himm_vae[,-1], as.numeric)
PXD000279_hilm_vae[,-1] <- sapply(PXD000279_hilm_vae[,-1], as.numeric)
PXD000279_mihm_vae[,-1] <- sapply(PXD000279_mihm_vae[,-1], as.numeric)
PXD000279_mimm_vae[,-1] <- sapply(PXD000279_mimm_vae[,-1], as.numeric)
PXD000279_milm_vae[,-1] <- sapply(PXD000279_milm_vae[,-1], as.numeric)
PXD000279_lihm_vae[,-1] <- sapply(PXD000279_lihm_vae[,-1], as.numeric)
PXD000279_limm_vae[,-1] <- sapply(PXD000279_limm_vae[,-1], as.numeric)
PXD000279_lilm_vae[,-1] <- sapply(PXD000279_lilm_vae[,-1], as.numeric)

#set sample id as row name
#cf
rownames(PXD000279_hilm_cf) <- PXD000279_hilm_cf[,1]
rownames(PXD000279_mimm_cf) <- PXD000279_mimm_cf[,1]
rownames(PXD000279_milm_cf) <- PXD000279_milm_cf[,1]
rownames(PXD000279_limm_cf) <- PXD000279_limm_cf[,1]
rownames(PXD000279_lilm_cf) <- PXD000279_lilm_cf[,1]

#vae
rownames(PXD000279_hihm_vae) <- PXD000279_hihm_vae[,1]
rownames(PXD000279_himm_vae) <- PXD000279_himm_vae[,1]
rownames(PXD000279_hilm_vae) <- PXD000279_hilm_vae[,1]
rownames(PXD000279_mihm_vae) <- PXD000279_mihm_vae[,1]
rownames(PXD000279_mimm_vae) <- PXD000279_mimm_vae[,1]
rownames(PXD000279_milm_vae) <- PXD000279_milm_vae[,1]
rownames(PXD000279_lihm_vae) <- PXD000279_lihm_vae[,1]
rownames(PXD000279_limm_vae) <- PXD000279_limm_vae[,1]
rownames(PXD000279_lilm_vae) <- PXD000279_lilm_vae[,1]

#dae
rownames(PXD000279_hihm_dae) <- PXD000279_hihm_dae[,1]
rownames(PXD000279_himm_dae) <- PXD000279_himm_dae[,1]
rownames(PXD000279_hilm_dae) <- PXD000279_hilm_dae[,1]
rownames(PXD000279_mihm_dae) <- PXD000279_mihm_dae[,1]
rownames(PXD000279_mimm_dae) <- PXD000279_mimm_dae[,1]
rownames(PXD000279_milm_dae) <- PXD000279_milm_dae[,1]
rownames(PXD000279_lihm_dae) <- PXD000279_lihm_dae[,1]
rownames(PXD000279_limm_dae) <- PXD000279_limm_dae[,1]
rownames(PXD000279_lilm_dae) <- PXD000279_lilm_dae[,1]

#transpose
#cf
PXD000279_hilm_cf <- as.data.frame(t(PXD000279_hilm_cf[,-1]))
PXD000279_milm_cf <- as.data.frame(t(PXD000279_milm_cf[,-1]))
PXD000279_limm_cf <- as.data.frame(t(PXD000279_limm_cf[,-1]))
PXD000279_lilm_cf <- as.data.frame(t(PXD000279_lilm_cf[,-1]))

#vae
PXD000279_hihm_vae <- as.data.frame(t(PXD000279_hihm_vae[,-1]))
PXD000279_himm_vae <- as.data.frame(t(PXD000279_himm_vae[,-1]))
PXD000279_hilm_vae <- as.data.frame(t(PXD000279_hilm_vae[,-1]))
PXD000279_mihm_vae <- as.data.frame(t(PXD000279_mihm_vae[,-1]))
PXD000279_mimm_vae <- as.data.frame(t(PXD000279_mimm_vae[,-1]))
PXD000279_milm_vae <- as.data.frame(t(PXD000279_milm_vae[,-1]))
PXD000279_lihm_vae <- as.data.frame(t(PXD000279_lihm_vae[,-1]))
PXD000279_limm_vae <- as.data.frame(t(PXD000279_limm_vae[,-1]))
PXD000279_lilm_vae <- as.data.frame(t(PXD000279_lilm_vae[,-1]))

#dae
PXD000279_hihm_dae <- as.data.frame(t(PXD000279_hihm_dae[,-1]))
PXD000279_himm_dae <- as.data.frame(t(PXD000279_himm_dae[,-1]))
PXD000279_hilm_dae <- as.data.frame(t(PXD000279_hilm_dae[,-1]))
PXD000279_mihm_dae <- as.data.frame(t(PXD000279_mihm_dae[,-1]))
PXD000279_mimm_dae <- as.data.frame(t(PXD000279_mimm_dae[,-1]))
PXD000279_milm_dae <- as.data.frame(t(PXD000279_milm_dae[,-1]))
PXD000279_lihm_dae <- as.data.frame(t(PXD000279_lihm_dae[,-1]))
PXD000279_limm_dae <- as.data.frame(t(PXD000279_limm_dae[,-1]))
PXD000279_lilm_dae <- as.data.frame(t(PXD000279_lilm_dae[,-1]))

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

#impute as whole
PXD000279_whole_na <- rbind(PXD000279_hihm_na, PXD000279_himm_na, PXD000279_hilm_na,
                            PXD000279_mihm_na, PXD000279_mimm_na, PXD000279_milm_na,
                            PXD000279_lihm_na, PXD000279_limm_na, PXD000279_lilm_na)
#real
PXD000279_whole_real <- rbind(PXD000279_hihm_real[rownames(PXD000279_hihm_real) %in% rownames(PXD000279_hihm_na),],
                              PXD000279_himm_real[rownames(PXD000279_himm_real) %in% rownames(PXD000279_himm_na),],
                              PXD000279_hilm_real[rownames(PXD000279_hilm_real) %in% rownames(PXD000279_hilm_na),],
                              PXD000279_mihm_real[rownames(PXD000279_mihm_real) %in% rownames(PXD000279_mihm_na),],
                              PXD000279_mimm_real[rownames(PXD000279_mimm_real) %in% rownames(PXD000279_mimm_na),],
                              PXD000279_milm_real[rownames(PXD000279_milm_real) %in% rownames(PXD000279_milm_na),],
                              PXD000279_lihm_real[rownames(PXD000279_lihm_real) %in% rownames(PXD000279_lihm_na),],
                              PXD000279_limm_real[rownames(PXD000279_limm_real) %in% rownames(PXD000279_limm_na),],
                              PXD000279_lilm_real[rownames(PXD000279_lilm_real) %in% rownames(PXD000279_lilm_na),])
#bpca
PXD000279_whole_bpca <- pca(PXD000279_whole_na, method = "bpca", center = TRUE, nPcs = 5, completeObs = TRUE, scale = "vector")
PXD000279_whole_bpca <- as.data.frame(PXD000279_whole_bpca@completeObs)
#svd
PXD000279_whole_svd <- as.data.frame(impute.wrapper.SVD(PXD000279_whole_na, K = 10))
#knn
PXD000279_whole_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(PXD000279_whole_na), K = 5))
#mle
PXD000279_whole_mle <- as.data.frame(impute.wrapper.MLE(as.matrix(PXD000279_whole_na)))
#lls
PXD000279_whole_lls <- llsImpute(PXD000279_whole_na, k = 5, center = TRUE, completeObs = TRUE)
PXD000279_whole_lls <- as.data.frame(PXD000279_whole_lls@completeObs)
#rf
PXD000279_whole_rf <- missForest(PXD000279_whole_na, maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
PXD000279_whole_rf <- PXD000279_whole_rf$ximp

#for pimms
write.csv(t(PXD000279_whole_na), file = "PXD000279_whole_na.csv")
#import pimms result
#fread for fast read
PXD000279_whole_cf <- fread("PXD000279_whole_cf.csv", header = TRUE)
PXD000279_whole_dae <- fread("PXD000279_whole_dae.csv", header = TRUE)
PXD000279_whole_vae <- fread("PXD000279_whole_vae.csv", header = TRUE)

PXD000279_whole_cf <- as.data.frame(PXD000279_whole_cf)
PXD000279_whole_dae <- as.data.frame(PXD000279_whole_dae)
PXD000279_whole_vae <- as.data.frame(PXD000279_whole_vae)

#set colname
colnames(PXD000279_whole_dae) <- PXD000279_whole_dae[1,]
colnames(PXD000279_whole_vae) <- PXD000279_whole_vae[1,]

#remove excessive rows
PXD000279_whole_dae <- PXD000279_whole_dae[-c(1,2),]
PXD000279_whole_vae <- PXD000279_whole_vae[-c(1,2),]

#turn chr into num
PXD000279_whole_dae[,-1] <- sapply(PXD000279_whole_dae[,-1], as.numeric)
PXD000279_whole_vae[,-1] <- sapply(PXD000279_whole_vae[,-1], as.numeric)

#set sample id as row name
rownames(PXD000279_whole_cf) <- PXD000279_whole_cf[,1]
rownames(PXD000279_whole_dae) <- PXD000279_whole_dae[,1]
rownames(PXD000279_whole_vae) <- PXD000279_whole_vae[,1]

PXD000279_whole_cf <- as.data.frame(t(PXD000279_whole_cf[,-1]))
PXD000279_whole_dae <- as.data.frame(t(PXD000279_whole_dae[,-1]))
PXD000279_whole_vae <- as.data.frame(t(PXD000279_whole_vae[,-1]))

#assemble the optimal method matrix
PXD000279_whole_mix <- rbind(PXD000279_hihm_lls, PXD000279_himm_bpca, PXD000279_hilm_rf,
                             PXD000279_mihm_rf, PXD000279_mimm_bpca, PXD000279_milm_rf,
                             PXD000279_lihm_rf, PXD000279_limm_bpca, PXD000279_lilm_rf)

#get na data
comp_na <- PXD000279_whole_na
#get real data
comp_real <- PXD000279_whole_real[which(rownames(PXD000279_whole_real) %in% rownames(comp_na)),]

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

methods <- c("bpca", "knn", "lls", "mle", "svd", "rf", "cf", "dae", "vae", "mix")
id <- "PXD000279"
#get data for different methods
for (j in 1:length(methods)) {
  if (paste(id, "whole", methods[j], sep = "_") %in% ls() == TRUE) {
    temp <- get(paste(id, "whole", methods[j], sep = "_"))
    for (k in 1:length(nrmse_sum$sample)) {
      nrmse_sum[k, grep(methods[j], colnames(nrmse_sum))] <- sqrt(mean((temp[,k] - comp_real[,k])^2, na.rm = TRUE))/sd(comp_real[,k], na.rm = TRUE)
    }
  } else {
    next
  }
  #if highest value > 4 fold of median, remove values higher than 5*(50% values) in a method
  if (max(nrmse_sum[,grep(methods[j], colnames(nrmse_sum))]) > 4*quantile(nrmse_sum[,grep(methods[j], colnames(nrmse_sum))], probs = c(0,0.5,0.9,0.95,1))[2]) {
    print(methods[j])
    print(quantile(nrmse_sum[,grep(methods[j], colnames(nrmse_sum))], probs = c(0,0.5,0.9,0.95,1)))
    nrmse_sum[which(nrmse_sum[,grep(methods[j], colnames(nrmse_sum))] > 3*quantile(nrmse_sum[,grep(methods[j], colnames(nrmse_sum))], probs = c(0,0.5,0.9,0.95,1))[2]),
              grep(methods[j], colnames(nrmse_sum))] <- NA
  }
}


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
  geom_text(aes(label = round(mean, digits = 2)), vjust = -1.5, size = 5) +
  scale_x_discrete(labels = c("RF", "BPCA", "kNN", "LLS", "MLE", "SVD", "CF", "DAE", "VAE", "Mix")) +
  labs(x = "Method",
       y = "NRMSE") +
  ylim(0,0.4) +
  scale_fill_brewer(labels = c("RF", "BPCA", "kNN", "LLS", "MLE", "SVD", "CF", "DAE", "VAE", "Mix"), palette = "Set3") +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 15, color = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.major = element_line(linewidth = 1),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15))

ggsave("figure7B_mix_compare_PXD000279.tiff",
       plot = last_plot(),
       path = "figures/manuscript/",
       device = "tiff",
       width = 2244,
       height = 950,
       units = "px",
       dpi = 300)
ggsave("figure7B_mix_compare_PXD000279.jpeg",
       plot = last_plot(),
       path = "figures/manuscript/",
       device = "jpeg",
       width = 2244,
       height = 950,
       units = "px",
       dpi = 300)

##differential expression analysis
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

#DE by t-test

for (i in 1:length(DE_result$df)) {
  data <- get(DE_result$df[i])
  temp <- data.frame(feature = rownames(data),
                     logFC = NA,
                     t = NA,
                     p = NA,
                     adjp = NA)
  temp$logFC <- rowMeans(data[,grep(" H", colnames(data))]) - rowMeans(data[,grep(" L", colnames(data))])
  for (j in 1:length(temp)) {

    
  }

  
}




ggplot(data = DE_result, aes(x = method, y = sig_DEP, fill = method)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sig_DEP, vjust = -0.5)) +
  scale_fill_discrete(labels = c("Real", "RF", "BPCA", "kNN", "LLS", "MLE", "SVD", "CF", "DAE", "VAE", "Mix")) +
  scale_x_discrete(labels = c("Real", "RF", "BPCA", "kNN", "LLS", "MLE", "SVD", "CF", "DAE", "VAE", "Mix")) +
  labs(x = "Method",
       y = "Number of significant differentially expressed peptides") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 15),
        axis.text = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13))

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
ROC_result <- data.frame(method = DE_result$method,
                         df = DE_result$sig_DF,
                         TP = NA,
                         TN = NA,
                         FP = NA,
                         FN = NA,
                         TPR = NA,
                         FPR = NA)

total <- length(rownames(PXD000279_whole_na))
NP_H <- table(ppmap$species)[2]
NP_E <- table(ppmap$species)[1]


for (i in 1:length(ROC_result$method)) {
  temp <- get(ROC_result$df[i])
  ROC_result$TP[i] <- length(which(rownames(temp) %in% ppmap$peptide[ppmap$species == "E.coli"] == TRUE))
  ROC_result$TN[i] <- NP_H - length(which(rownames(temp) %in% ppmap$peptide[ppmap$species == "Human"] == TRUE))
  ROC_result$FP[i] <- length(which(rownames(temp) %in% ppmap$peptide[ppmap$species == "Human"] == TRUE))
  ROC_result$FN[i] <- NP_E - length(which(rownames(temp) %in% ppmap$peptide[ppmap$species == "E.coli"] == TRUE))
}

ROC_result$TPR <- ROC_result$TP/(ROC_result$TP + ROC_result$FN)
ROC_result$FPR <- ROC_result$FP/(ROC_result$FP + ROC_result$TN)

comp_real <- PXD000279_whole_real %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "real")
comp_na <- PXD000279_whole_na %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "na")
comp_knn <- PXD000279_whole_knn %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "knn")
comp_mle <- PXD000279_whole_mle %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "mle")
comp_svd <- PXD000279_whole_svd %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "svd")
comp_bpca <- PXD000279_whole_bpca %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "bpca")
comp_lls <- PXD000279_whole_lls %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "lls")
comp_rf <- PXD000279_whole_rf %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "rf")
comp_cf <- PXD000279_whole_cf %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "cf")
comp_vae <- PXD000279_whole_vae %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "vae")
comp_dae <- PXD000279_whole_dae %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "dae")
comp_mix <- PXD000279_whole_mix %>% pivot_longer(., 1:length(.), names_to = "sample", values_to = "mix")

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

library(ggplot2)
library(pROC)
#define object to plot
rocobj <- roc(comp$real, comp$mix)
#create ROC plot
ggroc(rocobj)
auc(comp$real, comp$mix)

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


id

sig
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

rocobj <- roc(roc_plot$truth, roc_plot$mix)
#create ROC plot
ggroc(rocobj)



library(plotROC)
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
ggsave("figure6_roc.jpeg",
       plot = last_plot(),
       path = "figures/manuscript/",
       device = "jpeg",
       width = 2244,
       height = 1496,
       units = "px",
       dpi = 300)

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