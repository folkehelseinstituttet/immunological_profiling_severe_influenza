
unsup_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "Cytof unsupperviced")
panel <- "Panel 2"
subData <- "All"
seed <- "seed 2234"

ut_sti <- fs::path(unsup_path, panel, subData, seed)

fig_path <- fs::path(unsup_path, "Endelig_des2022", "Figurer")

temp <- readRDS(fs::path(ut_sti, "dataUMAP_tSNE_5000.RDS"))  #umap = umapAshinfac5$layout,
#temp2 <- readRDS(fs::path(ut_sti, "data_med_filnavn1345.RDS"))  #umap = umapAshinfac5$layout,

dtSNE <- data.frame(temp$res_tsne)
colnames(dtSNE) <- c("tSNE1", "tSNE2")

d <- temp$d

signif <- readxl::read_excel(fs::path(unsup_path, "Endelig_des2022", "Figurer", "PANEL2_ST1vsMT1_FINAL.xlsx"))
signif$Navn_tSNE <- signif$Navn_V2
signif <- signif[order(signif$`Kluster nr`),]
Cluster <- signif$`Anjas navnPercent`
groups <- signif$Navn_tSNE

color_group <- c("#CCCCCC", signif$colour)
names(color_group) <- c("Cl0 (Other cells, ns)",  signif$Navn_tSNE)



# 
# 
# d$celler <- "Cl Cl0: (Rest)"
# d$celler[d$`kluster 20` == 10] <- "Cl 1 (CD4 CM)"  
# table(d$celler)
# d$celler[d$`kluster 40` == 31] <- "Cl 2 (CD4 CM PD-L2)"  
# table(d$celler)
# d$celler[d$`kluster 20` == 8] <- "Cl 3 (CD4 CM CD137)"  
# table(d$celler)
# d$celler[d$`kluster 30` == 10] <- "Cl 4 (CD8 EM PD-1)"  
# table(d$celler)
# d$celler[d$`kluster 50` == 19] <-  "Cl 5 (CD8 TEMRA GranzB)"
# table(d$celler)
# d$celler[d$`kluster 30` == 3] <- "Cl 6 (NK GranzB & IFNg)"  
# table(d$celler)
# d$celler[d$`kluster 50` == 5] <-  "Cl 7 (NK GranzB)"
# table(d$celler)
# d$celler[d$`kluster 20` == 4] <-  "Cl 8 (NK GranzB & TIM-3)"
# table(d$celler)
# 


d$celler <- "Cl0 (Other cells, ns)"
#d$celler[d$`kluster 10` == 3] <- "10_3 (CD4)"  
#table(d$celler)
for(i in 1:nrow(signif)){
  klus <- strsplit(signif$`Anjas navnPercent`[i],"_")
  k <- klus[[1]][5]
  nr <- klus[[1]][7]
  col_k <- which(colnames(d) == paste("kluster", k))
  d$celler[d[,col_k] == nr] <- signif$Navn_tSNE[i]
}


dtSNE <- cbind(dtSNE, d$celler)
colnames(dtSNE)[3] <- "clusters"
dtSNE$clusters <- factor(dtSNE$clusters, c("Cl0 (Other cells, ns)", as.character(groups)))


Min1 <- min(dtSNE$tSNE1)
Max1 <- max(dtSNE$tSNE1)
Min2 <- min(dtSNE$tSNE2)
Max2 <- max(dtSNE$tSNE2)


table(d$status.x)

library(ggplot2)
library(ggpubr)
library(RColorBrewer)

taMed <- sample(which(d$status.x == "Control"), 25000)
gC <- ggplot(dtSNE[taMed,], aes(x = tSNE1, y = tSNE2, col = clusters)) + 
  geom_point() +
  xlim(Min1,Max1) +
  ylim(Min2,Max2) +
  ggtitle(paste0("Control")) + 
  scale_colour_manual(values = as.character(color_group), breaks = names(color_group)) + 
  theme_classic(base_size = 10) +
  #scale_color_brewer(palette = "Paired") +
  theme(legend.position = "none")

gC2 <- ggplot(dtSNE[taMed,], aes(x = tSNE1, y = tSNE2, col = clusters)) + 
  geom_point(shape = 15, size = 4) +
  xlim(Min1,Max1) +
  ylim(Min2,Max2) +
  ggtitle(paste0("Control")) + 
  scale_colour_manual(values = as.character(color_group), breaks = names(color_group)) + 
  theme_classic(base_size = 10) +
  #scale_color_brewer(palette = "Paired") +
  theme(legend.text=element_text(size=13))  +
  theme(legend.title=element_blank()) + guides(col=guide_legend(ncol=1))

leg <- get_legend(gC2)
gleg <- as_ggplot(leg)




taMed <- sample(which(d$status.x == "Severe T1"), 25000)
gS1 <- ggplot(dtSNE[taMed,], aes(x = tSNE1, y = tSNE2, col = clusters)) + 
  geom_point() +
  xlim(Min1,Max1) +
  ylim(Min2,Max2) +
  ggtitle(paste0("Severe T1")) + 
  theme_classic(base_size = 10) +
  scale_colour_manual(values = as.character(color_group), breaks = names(color_group)) + 
  #scale_color_brewer(palette = "Paired") +
  theme(legend.position = "none")




taMed <- sample(which(d$status.x == "Severe T2"), 25000)
gS2 <- ggplot(dtSNE[taMed,], aes(x = tSNE1, y = tSNE2, col = clusters)) + 
  geom_point() +
  xlim(Min1,Max1) +
  ylim(Min2,Max2) +
  ggtitle(paste0("Severe T2")) + 
  theme_classic(base_size = 10) +
  scale_colour_manual(values = as.character(color_group), breaks = names(color_group)) + 
  #scale_color_brewer(palette = "Paired") +
  theme(legend.position = "none")



taMed <- sample(which(d$status.x == "Moderate T1"), 25000)
gM1 <- ggplot(dtSNE[taMed,], aes(x = tSNE1, y = tSNE2, col = clusters)) + 
  geom_point() +
  xlim(Min1,Max1) +
  ylim(Min2,Max2) +
  ggtitle(paste0("Moderate T1")) +
  theme_classic(base_size = 10) +
  scale_colour_manual(values = as.character(color_group), breaks = names(color_group)) + 
  #scale_color_brewer(palette = "Paired") +
  theme(legend.position = "none")




taMed <- sample(which(d$status.x == "Moderate T2"), 25000)
gM2 <- ggplot(dtSNE[taMed,], aes(x = tSNE1, y = tSNE2, col = clusters)) + 
  geom_point() +
  xlim(Min1,Max1) +
  ylim(Min2,Max2) +
  ggtitle(paste0("Moderate T2")) +
  theme_classic(base_size = 10) +
  scale_colour_manual(values = as.character(color_group), breaks = names(color_group)) + 
  #scale_color_brewer(palette = "Paired") +
  theme(legend.position = "none")




tiff(fs::path(fig_path, paste0("Panel2_tSNE_unsup_group_per_status_tid.tiff")), width = 500, height = 700)
ggarrange(gS1, gS2, gM1, gM2, gleg, gC, ncol = 2, nrow = 3)
dev.off()




#alle i en

g <- ggplot(dtSNE, aes(x = tSNE1, y = tSNE2, col = clusters)) + 
  geom_point() +
  xlim(Min1,Max1) +
  ylim(Min2,Max2) +
  ggtitle(paste0("")) + 
  theme_classic(base_size = 20) +
  scale_colour_manual(values = as.character(color_group), breaks = names(color_group)) + 
  #scale_color_brewer(palette = "Paired") +
  theme(legend.position = "none")

  tiff(fs::path(fig_path, paste0("Panel2_tSNE_unsup_group_alle_uten_legend.tiff")), width = 700, height = 700)
g
dev.off()

tiff(fs::path(fig_path, paste0("Panel2_tSNE_unsup_group_alle_legend.tiff")), width = 700, height = 700)
gleg
dev.off()



fra <- which(colnames(d) == "89Y_CD45")
til <- which(colnames(d) == "146Nd_TNFa")

signal0 <- d[d$celler == "Cl0 (Other cells, ns)" ,fra:til]

q0 <- apply(signal0, 2, quantile, probs = c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95))

write.csv2(q0, fs::path(ut_sti, "quantilesAllOthers.csv"))


