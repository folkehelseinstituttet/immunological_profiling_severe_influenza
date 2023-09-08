library(ggplot2)
library(ggpubr)
library(RColorBrewer)

source("H:/git/cytof/functions.R")

unsup_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "Cytof unsupperviced")
panel <- "Panel 1"
subData <- "All"
seed <- "seed 1345"

ut_sti <- fs::path(unsup_path, panel, subData, seed)

fig_path <- fs::path(unsup_path, "Endelig_des2022", "Figurer")

temp <- readRDS(fs::path(ut_sti, "dataUMAP_tSNE_5000.RDS"))  #umap = umapAshinfac5$layout,

# dUMAP <- data.frame(temp$umap)
# colnames(dUMAP) <- c("UMAP1", "UMAP2")
dtSNE <- data.frame(temp$res_tsne)
colnames(dtSNE) <- c("tSNE1", "tSNE2")


gates <- read.csv2(fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_panel1_mars2022", "posNeg", "Data", "gater.csv"))
navn <- read.csv2( fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_panel1_mars2022", "posNeg", "Data", "navn.csv" ))
navn$marker <- paste0("X", navn$marker)

gates <- merge(gates, navn, by = "marker_short_name")
rownames(gates) <- gates$marker




fra <- which(colnames(temp$d) == "89Y_CD45")
til <- which(colnames(temp$d) == "151Eu_CD123")
tempMat <- temp$d[,fra:til]






xx <- colnames(tempMat)
xx <- gsub("-", ".", xx)
nyttNavn <-  gates[paste0("X", xx), "marker_short_name"]
tabort <- xx[is.na(nyttNavn)]
tabort


if(length(tabort) > 0){
  tempMat <- tempMat[, colnames(tempMat)[!xx %in% tabort]]
  xx <- colnames(tempMat)
  xx <- gsub("-", ".", xx)
  nyttNavn <-  gates[paste0("X", xx), "marker_short_name"]
}

colnames(tempMat) <-   nyttNavn
posNegMat <- as.data.frame(matrix(NA, ncol = ncol(tempMat), nrow = nrow(tempMat)))

colnames(posNegMat) <- colnames(tempMat)

for(i in 1:ncol(tempMat)){
  x <- gates[gates$marker_short_name == colnames(tempMat)[i],]
  posNegMat[,i] <- tempMat[,i] > x$low
  if(!is.na(x$high)){
    posNegMat[ tempMat[,i] > x$high,i] <- 2
  }
}






posNegMat$cellname_Bagwell <- "unknown"
posNegMat$cellname_Bagwell[posNegMat$CD45 == 1 & posNegMat$CD19 == 0 & posNegMat$CD3 == 1 & posNegMat$CD56 %in% c(2) & posNegMat$TCRgd == 0] <- "NKT"
table(posNegMat$cellname_Bagwell)
posNegMat$cellname_Bagwell[posNegMat$CD45 == 1 & posNegMat$CD19 == 0 & posNegMat$CD3 == 1 & posNegMat$CD56 %in% c(1,2) & posNegMat$TCRgd == 0 & posNegMat$CD57 %in% c(2)] <- "NKT"
table(posNegMat$cellname_Bagwell)
posNegMat$cellname_Bagwell[posNegMat$CD45 == 1 & posNegMat$CD19 == 0 & posNegMat$CD3 == 1 & posNegMat$CD4 == 1 & posNegMat$CD8 == 0] <- "CD4" 
table(posNegMat$cellname_Bagwell)
posNegMat$cellname_Bagwell[posNegMat$CD45 == 1 & posNegMat$CD19 == 0 & posNegMat$CD3 == 1 & posNegMat$CD4 == 0 & posNegMat$CD8 %in% c(1,2)] <- "CD8"
table(posNegMat$cellname_Bagwell)
posNegMat$cellname_Bagwell[posNegMat$CD45 == 1 & posNegMat$CD19 == 0 & posNegMat$CD3 == 1 & posNegMat$TCRVa7.2 %in% c(0) & posNegMat$TCRgd %in% c(1)] <- "gdT"
table(posNegMat$cellname_Bagwell)
posNegMat$cellname_Bagwell[posNegMat$CD45 == 1 & posNegMat$CD19 == 0 & posNegMat$CD3 == 0 & (posNegMat$CD56 %in% c(1,2) ) & posNegMat$CD123 == 0 ] <- "NK"
table(posNegMat$cellname_Bagwell)
posNegMat$cellname_Bagwell[posNegMat$CD45 == 1 & posNegMat$CD19 == 1 & posNegMat$CD3 == 0 & posNegMat$CD56 == 0 & posNegMat$HLADR %in% c(1,2)] <- "B"
table(posNegMat$cellname_Bagwell)
posNegMat$cellname_Bagwell[posNegMat$CD45 == 1 & posNegMat$CD19 == 0 & posNegMat$CD3 == 0 & posNegMat$CD57 == 0 & posNegMat$HLADR %in% c(2) & posNegMat$CD16 == 0 & posNegMat$CD169 == 0] <- "DC"
table(posNegMat$cellname_Bagwell)
posNegMat$cellname_Bagwell[posNegMat$CD45 == 1 & posNegMat$CD19 == 0 & posNegMat$CD3 == 0 & posNegMat$CD56 %in% c(0) & posNegMat$CD57 == 0 & posNegMat$CD11c %in% c(1,2) &  posNegMat$CD11b == 1 &  (posNegMat$CD14 %in% c(1,2)  )] <- "Mo"
table(posNegMat$cellname_Bagwell)
posNegMat$cellname_Bagwell[posNegMat$CD45 == 1 & posNegMat$CD19 == 0 & posNegMat$CD3 == 0 & posNegMat$HLADR %in% c(1,2) & posNegMat$CD11c %in% c(1,2) & posNegMat$CD11b %in% c(1)  &  posNegMat$CD56 == 0 & posNegMat$CD14 == 0 & posNegMat$CD38 == 0 ] <- "Mo" 
table(posNegMat$cellname_Bagwell)
posNegMat$cellname_Bagwell[posNegMat$CD45 == 1 & posNegMat$CD3 == 0 & posNegMat$CD19 == 0 & posNegMat$CD123 == 1 & posNegMat$CD11c == 0 & posNegMat$HLADR %in% c(1,2) & posNegMat$CD14 == 0  & posNegMat$CD169 == 0] <- "DC" 
table(posNegMat$cellname_Bagwell)
posNegMat$cellname_Bagwell[posNegMat$CD45 == 1 & posNegMat$CD3 == 0 & posNegMat$CD19 == 0 & posNegMat$CD123 == 0 & posNegMat$CD11c %in% c(1,2) & posNegMat$HLADR %in% c(1,2) & posNegMat$CD14 == 0 & posNegMat$CD16 %in% c(1,0) & posNegMat$CD38 %in% c(1,0) & posNegMat$CRTH2 == 0 & posNegMat$CD169 == 0] <- "DC" 
table(posNegMat$cellname_Bagwell)
posNegMat$cellname_Bagwell[posNegMat$CD45 == 1 & posNegMat$CD19 == 0 & posNegMat$CD3 == 0 & posNegMat$CD56 == 0 & posNegMat$CD169 == 1] <- "Mo" 
table(posNegMat$cellname_Bagwell)
posNegMat$cellname_Bagwell[posNegMat$CD45 == 1 & posNegMat$CD19 == 0 & posNegMat$CD3 == 1 & posNegMat$TCRVa7.2 %in% c(1) & posNegMat$TCRgd %in% c(0)] <- "MAIT"
table(posNegMat$cellname_Bagwell)


table(posNegMat$cellname_Bagwell)/nrow(posNegMat)

unique(posNegMat$cellname_Bagwell)
#color_tablet_Bagwell <- c("B" = "brown", "CD4" = "yellow",  "Treg" = "orange", "MAIT" = "red", "CD8" = "blue", "Mo" = "light blue",
#                  "NKT" = "darkgreen", "NK" = "green", "DC"  = "purple", 
#                 "gdT" ="cyan")


posNegMat$cellname_Bagwell <- factor(posNegMat$cellname_Bagwell, levels = c("CD4", "CD8", "NK", "DC", "gdT", "MAIT", "B", "Mo", "NKT"))#, "Tregs"))
# 
# temp <- dUMAP
# dUMAP <- cbind(temp, posNegMat$cellname_Bagwell)
# colnames(dUMAP)[3] <- "cellname_Bagwell"



dtSNE_all <- dtSNE
dtSNE <- cbind(dtSNE_all, posNegMat$cellname_Bagwell)
colnames(dtSNE)[3] <- "cellname"

dtSNE <- dtSNE[!is.na(dtSNE$cellname), ]

get

brewer.pal(n = 9, name = "Paired")

# color_group <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6")
# names(color_group) <-  c("CD4", "CD8", "NK", "DC", "gdT", "MAIT", "B", "Mo", "NKT")



color_group <- c("#A6CEE3", "#1F78B4", "#FB9A99", "#E31A1C", "#FDBF6F", "#B2DF8A", "#CAB2D6", "#FF7F00", "#33A02C" )
names(color_group) <-  c("CD4", "CD8",  "gdT", "MAIT", "B", "NK", "NKT", "Mo", "DC")



#alle i en

g <- ggplot(dtSNE, aes(x = tSNE1, y = tSNE2, col = cellname)) + 
  geom_point() +
#  xlim(Min1,Max1) +
#  ylim(Min2,Max2) +
  ggtitle(paste0("")) + 
  theme_classic(base_size = 20) +
  #scale_fill_manual(limits = c("CD4", "CD8",  "gdT", "MAIT", "B", "NK", "NKT", "Mo", "DC")) + #CD4, CD8, gdT, MAIT, B, NK, NKT, Mo, DC
  scale_colour_manual(values = as.character(color_group), breaks = names(color_group), limits = c("CD4", "CD8",  "gdT", "MAIT", "B", "NK", "NKT", "Mo", "DC"))  +
#  scale_color_brewer(palette = "Paired")  #+
  theme(legend.position = "none")
  
  tiff(fs::path(fig_path, paste0("Panel1_tSNE_5000_manuel_manuel_uten_legend.tiff")), width = 700, height = 700)
g
dev.off()


tiff(fs::path(fig_path, "Manuel_ManuelGating", paste0("Panel1_tSNE_5000_manuel_manuel_uten_legend.tiff")), width = 700, height = 700)
g
dev.off()





g2 <- ggplot(dtSNE, aes(x = tSNE1, y = tSNE2, col = cellname)) + 
  geom_point(shape = 15, size = 5) +
  #  xlim(Min1,Max1) +
  #  ylim(Min2,Max2) +
  ggtitle(paste0("")) + 
  theme(legend.text=element_text(size=13))  +
  #scale_fill_manual(limits = c("CD4", "CD8",  "gdT", "MAIT", "B", "NK", "NKT", "Mo", "DC")) + #CD4, CD8, gdT, MAIT, B, NK, NKT, Mo, DC
  scale_colour_manual(values = as.character(color_group), breaks = names(color_group), limits = c("CD4", "CD8",  "gdT", "MAIT", "B", "NK", "NKT", "Mo", "DC"))  +
  #  scale_color_brewer(palette = "Paired")  #+
  theme(legend.title=element_blank())

leg <- get_legend(g2)
gleg <- as_ggplot(leg)
gleg 
  
tiff(fs::path(fig_path, paste0("Panel1_tSNE_5000_manuel_manuel_kun_legend.tiff")), width = 700, height = 700)
gleg
dev.off()


tiff(fs::path(fig_path, "Manuel_ManuelGating", paste0("Panel1_tSNE_5000_manuel_manuel_kun_legend.tiff")), width = 700, height = 700)
gleg
dev.off()



# 
# 
# fig_path_k <- fs::path(fig_path, "per_k")
# 
# k = 10
# 
# dtSNE <- cbind(dtSNE_all, temp$d$"kluster 10")
# colnames(dtSNE)[3] <- "k_10"
# dtSNE$k_10 <- factor(dtSNE$k_10)
# 
# 
# 
# #alle i en
# 
# g <- ggplot(dtSNE, aes(x = tSNE1, y = tSNE2, col = k_10)) + 
#   geom_point() +
#   #  xlim(Min1,Max1) +
#   #  ylim(Min2,Max2) +
#   ggtitle(paste0("")) + 
#   theme_classic(base_size = 20) +
#   scale_colour_manual(values = col50)  
#  
# 
# tiff(fs::path(fig_path_k, paste0("Panel1_tSNE_5000_k_10.tiff")), width = 900, height = 700)
# g
# dev.off()
# 
# 
# 
# k = 20
# 
# dtSNE <- cbind(dtSNE_all, temp$d$"kluster 20")
# colnames(dtSNE)[3] <- "k_20"
# dtSNE$k_20 <- factor(dtSNE$k_20)
# 
# 
# #alle i en
# 
# g <- ggplot(dtSNE, aes(x = tSNE1, y = tSNE2, col = k_20)) + 
#   geom_point() +
#   #  xlim(Min1,Max1) +
#   #  ylim(Min2,Max2) +
#   ggtitle(paste0("")) + 
#   theme_classic(base_size = 20) +
#   scale_colour_manual(values = col50) 
# 
# tiff(fs::path(fig_path_k, paste0("Panel1_tSNE_5000_k_20.tiff")), width = 900, height = 700)
# g
# dev.off()
# 
# 
# 
# 
# 
# k = 30
# 
# dtSNE <- cbind(dtSNE_all, temp$d$"kluster 30")
# colnames(dtSNE)[3] <- "k_30"
# dtSNE$k_30 <- factor(dtSNE$k_30)
# 
# 
# 
# #alle i en
# 
# g <- ggplot(dtSNE, aes(x = tSNE1, y = tSNE2, col = k_30)) + 
#   geom_point() +
#   #  xlim(Min1,Max1) +
#   #  ylim(Min2,Max2) +
#   ggtitle(paste0("")) + 
#   theme_classic(base_size = 20)  +
#   scale_colour_manual(values = col50) 
# 
# 
# tiff(fs::path(fig_path_k, paste0("Panel1_tSNE_5000_k_30.tiff")), width = 900, height = 700)
# g
# dev.off()
# 
# k = 40
# 
# dtSNE <- cbind(dtSNE_all, temp$d$"kluster 40")
# colnames(dtSNE)[3] <- "k_40"
# dtSNE$k_40 <- factor(dtSNE$k_40)
# 
# 
# 
# #alle i en
# 
# g <- ggplot(dtSNE, aes(x = tSNE1, y = tSNE2, col = k_40)) + 
#   geom_point() +
#   #  xlim(Min1,Max1) +
#   #  ylim(Min2,Max2) +
#   ggtitle(paste0("")) + 
#   theme_classic(base_size = 20)  +
#   scale_colour_manual(values = col50) 
# 
# 
# tiff(fs::path(fig_path_k, paste0("Panel1_tSNE_5000_k_40.tiff")), width = 900, height = 700)
# g
# dev.off()
# 
# 
# k = 50
# 
# dtSNE <- cbind(dtSNE_all, temp$d$"kluster 50")
# colnames(dtSNE)[3] <- "k_50"
# dtSNE$k_50 <- factor(dtSNE$k_50)
# 
# 
# 
# #alle i en
# 
# g <- ggplot(dtSNE, aes(x = tSNE1, y = tSNE2, col = k_50)) + 
#   geom_point() +
#   #  xlim(Min1,Max1) +
#   #  ylim(Min2,Max2) +
#   ggtitle(paste0("")) + 
#   theme_classic(base_size = 20)  +
#   scale_colour_manual(values = col50) 
# 
# 
# tiff(fs::path(fig_path_k, paste0("Panel1_tSNE_5000_k_50.tiff")), width = 900, height = 700)
# g
# dev.off()
# 
# 
# k = 60
# 
# dtSNE <- cbind(dtSNE_all, temp$d$"kluster 60")
# colnames(dtSNE)[3] <- "k_60"
# dtSNE$k_60 <- factor(dtSNE$k_60)
# 
# 
# #alle i en
# col100 <- c(col50, col50)
# names(col100) <- as.character(1:100)
# 
# g <- ggplot(dtSNE, aes(x = tSNE1, y = tSNE2, col = k_60)) + 
#   geom_point() +
#   #  xlim(Min1,Max1) +
#   #  ylim(Min2,Max2) +
#   ggtitle(paste0("")) + 
#   theme_classic(base_size = 20)  +
#   scale_colour_manual(values = col100)  +
#   theme(legend.position = "none")
# 
# 
# 
# tiff(fs::path(fig_path_k, paste0("Panel1_tSNE_5000_k_60.tiff")), width = 900, height = 700)
# g
# dev.off()
# 
# 
# 
# 
# # Denne kjC8rer i figUMAPp1.R trenger ikke begge plasser, men sparer hvis fil slettes...
# # 
# # #per markC8r og status
# # 
# # 
# # dUMAP <- data.frame(temp$umap)
# # colnames(dUMAP) <- c("UMAP1", "UMAP2")
# # dtSNE <- data.frame(temp$res_tsne)
# # colnames(dtSNE) <- c("tSNE1", "tSNE2")
# # d <- data.frame(temp$d2_mindre)
# # 
# # 
# # navn <- read.csv2(fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel1_mars2022", "posNeg", "Data", "navn.csv"))
# # 
# # pal <- rev(brewer.pal(10,"Spectral"))  #wes_palette("Zissou1", 100, type = "continuous")
# # 
# # 
# # library(ggplot2)
# # library(ggpubr)
# # 
# # 
# # for(i in 1:ncol(d)){
# #   marker <- navn$V2[i]
# #   
# #   dtSNE$signal <- d[,i]
# #   dtSNE_i <- dtSNE
# #   q0001 <- quantile(dtSNE_i[,3], 0.001)
# #   q0999 <- quantile(dtSNE_i[,3], 0.999)
# #   
# #   dtSNE_i <- dtSNE_i[dtSNE_i[,3] <= q0999, ]
# #   dtSNE_i <- dtSNE_i[dtSNE_i[,3] >= q0001, ]
# #   
# #   maks_signal <- ceiling(max(dtSNE_i[,3]))
# #   min_signal <- floor(min(dtSNE_i[,3]))
# #   
# #   
# #   
# #   #tiff(fs::path(fig_path, "Per marker", paste0("Panel1_tSNE_", marker, "_control.tiff")), width = 500, height = 700)
# #   taMed <- sample(which(posNegMat$status == "Control"), 24909)
# #   g <- ggplot(dtSNE_i[taMed,], aes(x = tSNE1, y = tSNE2, col = signal)) + 
# #     geom_point() +
# #     ggtitle(paste0(marker, ", control")) + 
# #     scale_color_gradientn(colours =  pal,  limits = range(min_signal, maks_signal))+ 
# #     theme(legend.position = "none")
# #   #print(g)
# #   g2 <- g
# #   #dev.off()
# #   
# #   taMed <- sample(which(posNegMat$status == "Control"), 24909)
# #   gC2 <- ggplot(dtSNE_i[taMed,], aes(x = tSNE1, y = tSNE2, col = signal)) + 
# #     geom_point() +
# #     ggtitle(paste0(marker, ", control")) + 
# #     scale_color_gradientn(colours =  pal,  limits = range(min_signal, maks_signal))  +
# #     theme(legend.text=element_text(size=13))  +
# #     guides(col = guide_colourbar(barwidth = 3, barheight = 8)) +
# #     theme(legend.title=element_blank())
# #   # 
# #   leg <- get_legend(gC2)
# #   gleg <- as_ggplot(leg)
# #   
# #   
# #   
# #   
# #   #tiff(fs::path(fig_path, "Per marker", paste0("Panel1_tSNE_", marker, "_Severe T1_.tiff")), width = 500, height = 700)
# #   taMed <- sample(which(posNegMat$status == "Severe T1"), 24909)
# #   g <- ggplot(dtSNE_i[taMed,], aes(x = tSNE1, y = tSNE2, col = signal)) + 
# #     geom_point() + 
# #     ggtitle(paste0(marker, ", severe T1")) + 
# #     scale_color_gradientn(colours =  pal,  limits = range(min_signal, maks_signal))+ 
# #     theme(legend.position = "none")
# #   #print(g)
# #   g3 <- g
# #   #dev.off()
# #   
# #   #tiff(fs::path(fig_path, "Per marker", paste0("Panel1_tSNE_", marker, "_Severe T2.tiff")), width = 500, height = 700)
# #   taMed <- sample(which(posNegMat$status == "Severe T2"), 24909)
# #   g <- ggplot(dtSNE_i[taMed,], aes(x = tSNE1, y = tSNE2, col = signal)) +
# #     geom_point() + 
# #     ggtitle(paste0(marker, ", severe T2")) + 
# #     scale_color_gradientn(colours =  pal,  limits = range(min_signal, maks_signal))+ 
# #     theme(legend.position = "none")
# #   #print(g)
# #   g4 <- g
# #   #dev.off()
# #   
# #   #tiff(fs::path(fig_path, "Per marker", paste0("Panel1_tSNE_", marker, "_Moderate T1.tiff")), width = 500, height = 700)
# #   taMed <- sample(which(posNegMat$status == "Moderate T1"), 24909)
# #   g <- ggplot(dtSNE_i[taMed,], aes(x = tSNE1, y = tSNE2, col = signal)) + 
# #     geom_point() + 
# #     ggtitle(paste0(marker, ", moderate T1")) + 
# #     scale_color_gradientn(colours =  pal,  limits = range(min_signal, maks_signal))+ 
# #     theme(legend.position = "none")
# #   #print(g)
# #   g5 <- g
# #   #dev.off()
# #   
# #   #tiff(fs::path(fig_path, "Per marker", paste0("Panel1_tSNE_", marker, "_Moderate T2.tiff")), width = 500, height = 700)
# #   taMed <- sample(which(posNegMat$status == "Moderate T2"), 24909)
# #   g <- ggplot(dtSNE_i[taMed,], aes(x = tSNE1, y = tSNE2, col = signal)) + 
# #     geom_point() + 
# #     ggtitle(paste0(marker, ", moderate T2")) + 
# #     scale_color_gradientn(colours =  pal,  limits = range(min_signal, maks_signal))+ 
# #     theme(legend.position = "none")
# #   #print(g)
# #   g6 <- g
# #   #dev.off()
# #   
# #   g <- ggarrange(g5, g6, g3, g4, gleg, g2, ncol = 2, nrow = 3)
# #   
# #   
# #   tiff(fs::path(fig_path, "Per marker samlet", paste0("Panel1_tSNE_", marker, ".tiff")), width = 500, height = 700)
# #   print(g)
# #   dev.off()
# #   
# #   dUMAP$signal <- d[,i]
# #   dUMAP_i <- dUMAP
# #   dUMAP_i <- dUMAP_i[dUMAP_i[,3] <= q0999, ]
# #   dUMAP_i <- dUMAP_i[dUMAP_i[,3] >= q0001, ]
# #   
# #   
# #   
# #   #tiff(fs::path(fig_path, "Per marker", paste0("Panel1_UMAP_", marker, "_all.tiff")), width = 500, height = 700)
# #   g1 <- ggplot(dUMAP_i, aes(x = UMAP1, y = UMAP2, col = signal)) + 
# #     geom_point() + 
# #     ggtitle(paste0(marker, ", all")) + 
# #     scale_color_gradientn(colours =  pal,  limits = range(min_signal, maks_signal))+ 
# #     theme(legend.position = "none")
# #   #print(g)
# #   #dev.off()
# #   
# #   
# #   taMed <- sample(which(posNegMat$status == "Control"), 24909)
# #   gC2 <- ggplot(dUMAP_i[taMed,], aes(x = UMAP1, y = UMAP2, col = signal)) + 
# #     geom_point() +
# #     ggtitle(paste0(marker, ", control")) + 
# #     scale_color_gradientn(colours =  pal,  limits = range(min_signal, maks_signal))  +
# #     theme(legend.text=element_text(size=13))  +
# #     guides(col = guide_colourbar(barwidth = 3, barheight = 8)) +
# #     theme(legend.title=element_blank())
# #   # 
# #   leg <- get_legend(gC2)
# #   gleg <- as_ggplot(leg)
# #   
# #   
# #   
# #   #tiff(fs::path(fig_path, "Per marker", paste0("Panel1_UMAP_", marker, "_control.tiff")), width = 500, height = 700)
# #   taMed <- sample(which(posNegMat$status == "Control"), 24909)
# #   g2 <- ggplot(dUMAP_i[taMed,], aes(x = UMAP1, y = UMAP2, col = signal)) + 
# #     geom_point() +
# #     ggtitle(paste0(marker, ", control")) + 
# #     scale_color_gradientn(colours =  pal,  limits = range(min_signal, maks_signal))+ 
# #     theme(legend.position = "none")
# #   #print(g)
# #   #dev.off()
# #   
# #   
# #   #tiff(fs::path(fig_path, "Per marker", paste0("Panel1_UMAP_", marker, "_Severe T1.tiff")), width = 500, height = 700)
# #   taMed <- sample(which(posNegMat$status == "Severe T1"), 24909)
# #   g3 <- ggplot(dUMAP_i[taMed,], aes(x = UMAP1, y = UMAP2, col = signal)) + 
# #     geom_point() + 
# #     ggtitle(paste0(marker, ", severe T1")) + 
# #     scale_color_gradientn(colours =  pal,  limits = range(min_signal, maks_signal))+ 
# #     theme(legend.position = "none")
# #   #print(g)
# #   #dev.off()
# #   
# #   #tiff(fs::path(fig_path, "Per marker", paste0("Panel1_UMAP_", marker, "_Severe T2.tiff")), width = 500, height = 700)
# #   taMed <- sample(which(posNegMat$status == "Severe T2"), 24909)
# #   g4 <- ggplot(dUMAP_i[taMed,], aes(x = UMAP1, y = UMAP2, col = signal)) +
# #     geom_point() + 
# #     ggtitle(paste0(marker, ", severe T2")) + 
# #     scale_color_gradientn(colours =  pal,  limits = range(min_signal, maks_signal))+ 
# #     theme(legend.position = "none")
# #   #print(g)
# #   #dev.off()
# #   
# #   #tiff(fs::path(fig_path, "Per marker", paste0("Panel1_UMAP_", marker, "_Moderate T1.tiff")), width = 500, height = 700)
# #   taMed <- sample(which(posNegMat$status == "Moderate T1"), 24909)
# #   g5 <- ggplot(dUMAP_i[taMed,], aes(x = UMAP1, y = UMAP2, col = signal)) + 
# #     geom_point() + 
# #     ggtitle(paste0(marker, ", moderate T1")) + 
# #     scale_color_gradientn(colours =  pal,  limits = range(min_signal, maks_signal))+ 
# #     theme(legend.position = "none")
# #   #print(g)
# #   #dev.off()
# #   
# #   #tiff(fs::path(fig_path, "Per marker", paste0("Panel1_UMAP_", marker, "_Moderate T2.tiff")), width = 500, height = 700)
# #   taMed <- sample(which(posNegMat$status == "Moderate T2"), 24909)
# #   g6 <- ggplot(dUMAP_i[taMed,], aes(x = UMAP1, y = UMAP2, col = signal)) + 
# #     geom_point() + 
# #     ggtitle(paste0(marker, ", moderate T2")) + 
# #     scale_color_gradientn(colours =  pal,  limits = range(min_signal, maks_signal))+ 
# #     theme(legend.position = "none")
# #   #print(g)
# #   #dev.off()
# #   
# #   g <- ggarrange(g3, g4, g5, g6, gleg, g2, ncol = 2, nrow = 3)
# #   
# #   tiff(fs::path(fig_path, "Per marker samlet", paste0("Panel1_UMAP_", marker, ".tiff")), width = 500, height = 700)
# #   print(g)
# #   dev.off()
# #   
# #   
# # }
