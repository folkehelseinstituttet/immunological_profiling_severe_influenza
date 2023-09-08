
library(reshape2)
library(ggplot2)
library(ggrepel)
library(ggpubr)


pVerdiStjerne <- function(p){
  if(p < 0.001){
    return("***")
  } else {
    if(p < 0.01){
      return("**")
    } else { 
      if(p < 0.05){
        return("*")
      } else {
        return("-")
      }
    }
  }
}


cytofSignPlot <- function(d, sign_cl, signnavn, main = NULL, ymax = 20){
  d$sign <- FALSE
  d$sign[d$cluster %in% sign_cl] <- TRUE
  d$minuslog10adjp <- -log10(d$adj_p)
  d2 <- merge(d, signnavn, by.x = "cluster", by.y = "navn", all.x = TRUE)
  d2$kort_navn[d2$sign == FALSE] <- NA
  d2 <- d2[!(d2$adj_p < 0.05 & d2$sign == FALSE),]
  d2$diff <- "NO"
  d2$diff[d2$adj_p < 0.05 & d2$IRRmotsatt > 1] <- "UP"
  d2$diff[d2$adj_p < 0.05 & d2$IRRmotsatt < 1] <- "DOWN"
  d2$adj.p.value <- d2$adj_p
  cols <- c("DOWN" = "blue", "NO" = "black", "UP" = "red")
  
  # 
  # gg <- ggplot(d2, aes(x= log(IRR), y = -log10(adj.p.value), col = diff, label = kort_navn)) + 
  #   geom_point() + xlab("log(Incident rate ratio)") +
  #   theme_minimal() +
  #   geom_text_repel() +
  #   theme_classic(base_size = 15) +
  #   scale_color_manual(values=cols) +
  #   geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  #   ggtitle(main) + 
  #   theme(legend.position = "none") + 
  #   ylim(0,20) + xlim(-6,6) 
  
  
  gg <- ggplot(d2, aes(x= IRRmotsatt, y = -log10(adj.p.value), col = diff, label = kort_navn)) + 
    geom_point() + xlab("Incident rate ratio") +
    theme_minimal() +
    geom_text_repel( max.overlaps = Inf, size = 5) +
    theme_classic(base_size = 25) +
    scale_color_manual(values=cols) +
   # geom_vline(xintercept=c(2^-0.6, 2^0.6), col="red") +
    ggtitle(main) + 
    theme(legend.position = "none") + 
    ylim(0,ymax) +  coord_trans(x="log2", xlim = c(0.004,256)) +
    scale_x_continuous(breaks = c(0.004, 0.008 ,0.016, 0.031, 0.063, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256)) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  
  

  
  
  return(list(gg = gg, signcl = d2$kort_navn, x = d2$IRRmotsatt, label = d2$kort_navn, d = d2 ))
}



utSti <- "F:/Forskningsprosjekter/PDB 2794 - Immune responses aga_/Forskningsfiler/JOBO/CyTOF/Analyse i R OUS/CleanUpGatingMarch2022/Analyse/Endelig/Cytof unsupperviced/Endelig_des2022/Figurer/cytofSamlet"

d <- read.csv2("F:/Forskningsprosjekter/PDB 2794 - Immune responses aga_/Forskningsfiler/JOBO/CyTOF/Analyse i R OUS/CleanUpGatingMarch2022/Analyse/Endelig/Cytof unsupperviced/Endelig_des2022/Panel 1/All/seed 1345/IRRadjp_1345.csv")
d$IRRmotsatt <- 1/d$IRR
sign <- as.data.frame(readxl::read_excel("F:/Forskningsprosjekter/PDB 2794 - Immune responses aga_/Forskningsfiler/JOBO/CyTOF/Analyse i R OUS/CleanUpGatingMarch2022/Analyse/Endelig/Cytof unsupperviced/Endelig_des2022/Figurer/151222_Panel1_ST1vsMT1VsCST2vsMT2VsC_FINAL.xlsx"))
sign <- sign[!is.na(sign$Cluster),]
sign <- sign[sign$`ta med` == 1,]
sign_cl_ST1vsMT1 <- sign$`Anjas navnPercent`[!is.na(sign$S1vsM1)]
sign$navn[is.na(sign$navn)] <- sign$Cluster[is.na(sign$navn)]
sign2 <- sign[, c(3,2)]
colnames(sign2) <-  c("kort_navn", "navn")


ggST1vsMT1p1 <- cytofSignPlot(d = d, sign_cl = sign_cl_ST1vsMT1, signnavn = sign2, main =  "Severe T1 vs moderate T1")
ggST1vsMT1p1$gg
tiff(fs::path(utSti, "P1_cytof_sign_plot_ST1vsMT1.tiff"), width = 500, height = 700)
ggST1vsMT1p1$gg 
dev.off()

library(vegan)

d2 <- ggST1vsMT1p1$d
d2 <- d2[!is.na(d2$kort_navn),]
d2$col <- "red"
d2$col[d2$diff == "DOWN"] <- "blue"
d2 <- d2[order(d2$IRRmotsatt),]
d2$stjerner <- lapply(d2$adj_p, pVerdiStjerne)
d2$lab2 <- paste0(d2$kort_navn, " (", d2$stjerner, ")")


tiff(fs::path(utSti, "log2P1_IRR_linestack_ST1_MT1.tiff"), width = 300, height = 700)
par(mar = c(5, 4, 4, 2) + 0.1)
linestack(log2(d2$IRRmotsatt), d2$kort_navn, axis = T, air = 1.1, col = d2$col, main = "Severe T1 vs Moderate T1", ylim = c(-8,8))
dev.off()

d$lab2 <- paste0(d$delabel, " (", d$stjerner, ")")

tiff(fs::path(utSti, "log2P1_IRR_linestack_ST1_MT1_med stjerne.tiff"), width = 340, height = 700)
par(mar = c(5, 4, 4, 2) + 0.1)
linestack(log2(d2$IRRmotsatt), d2$lab2, axis = T, air = 1.1, col = d2$col, main = "Severe T1 vs Moderate T1", ylim = c(-8,8))
dev.off()


tiff(fs::path(utSti, "P1_IRR_linestack_ST1_MT1.tiff"), width = 300, height = 700)
linestack(log2(d2$IRRmotsatt), d2$kort_navn, axis = F, air = 1.1, col = d2$col, main = "Severe T1 vs Moderate T1", ylim = c(-8,8)) 
axis(side = 2, at = (-8):8, label = round(2^((-8):8),3))
dev.off()


tiff(fs::path(utSti, "P1_IRR_linestack_ST1_MT1_med stjerne.tiff"), width = 340, height = 700)
linestack(log2(d2$IRRmotsatt), d2$lab2, axis = F, air = 1.1, col = d2$col, main = "Severe T1 vs Moderate T1", ylim = c(-8,8)) 
axis(side = 2, at = (-8):8, label =round(2^((-8):8),3))
dev.off()



d <- read.csv2("F:/Forskningsprosjekter/PDB 2794 - Immune responses aga_/Forskningsfiler/JOBO/CyTOF/Analyse i R OUS/CleanUpGatingMarch2022/Analyse/Endelig/Cytof unsupperviced/Endelig_des2022/Panel 1/All/seed 1345/ST1vsCont/IRRadjp_1345.csv")
d$IRRmotsatt <- 1/d$IRR

sign_cl_ST1vsC <- sign$`Anjas navnPercent`[!is.na(sign$S1vsC)]
sign$navn[is.na(sign$navn)] <- sign$Cluster[is.na(sign$navn)]
sign2 <- sign[, c(3,2)]
colnames(sign2) <-  c("kort_navn", "navn")


ggST1vsCp1 <- cytofSignPlot(d = d, sign_cl = sign_cl_ST1vsC, signnavn = sign2, main =  "Severe T1 vs control")
ggST1vsCp1$gg


tiff(fs::path(utSti, "P1_cytof_sign_plot_ST1vsC_final.tiff"), width = 500, height = 700)
ggST1vsCp1$gg
dev.off()



d2 <- ggST1vsCp1$d
d2 <- d2[!is.na(d2$kort_navn),]
d2$col <- "red"
d2$col[d2$diff == "DOWN"] <- "blue"
d2 <- d2[order(d2$IRRmotsatt),]
d2$stjerner <- lapply(d2$adj_p, pVerdiStjerne)
d2$lab2 <- paste0(d2$kort_navn, " (", d2$stjerner, ")")


tiff(fs::path(utSti, "log2P1_IRR_linestack_ST1_control.tiff"), width = 360, height = 700)
par(mar = c(5, 4, 4, 2) + 0.1)
linestack(log2(d2$IRRmotsatt), d2$kort_navn, axis = T, air = 1.1, col = d2$col, main = "Severe T1 vs Control", ylim = c(-8,8))
dev.off()

d$lab2 <- paste0(d$delabel, " (", d$stjerner, ")")

tiff(fs::path(utSti, "log2P1_IRR_linestack_ST1_control_med stjerne.tiff"), width = 400, height = 700)
par(mar = c(5, 4, 4, 2) + 0.1)
linestack(log2(d2$IRRmotsatt), d2$lab2, axis = T, air = 1.1, col = d2$col, main = "Severe T1 vs controls", ylim = c(-8,8))
dev.off()


tiff(fs::path(utSti, "P1_IRR_linestack_ST1_control.tiff"), width = 360, height = 700)
linestack(log2(d2$IRRmotsatt), d2$kort_navn, axis = F, air = 1.1, col = d2$col, main = "Severe T1 vs controls", ylim = c(-8,8)) 
axis(side = 2, at = (-8):8, label = round(2^((-8):8),3))
dev.off()


tiff(fs::path(utSti, "P1_IRR_linestack_ST1_control_med stjerne.tiff"), width = 420, height = 700)
linestack(log2(d2$IRRmotsatt), d2$lab2, axis = F, air = 1.1, col = d2$col, main = "Severe T1 vs controls", ylim = c(-8,8)) 
axis(side = 2, at = (-8):8, label =round(2^((-8):8),3))
dev.off()



d <- read.csv2("F:/Forskningsprosjekter/PDB 2794 - Immune responses aga_/Forskningsfiler/JOBO/CyTOF/Analyse i R OUS/CleanUpGatingMarch2022/Analyse/Endelig/Cytof unsupperviced/Endelig_des2022/Panel 1/All/seed 1345/MT1vsCont/IRRadjp_1345.csv")
d$IRRmotsatt <- 1/d$IRR

sign_cl_MT1vsC <- sign$`Anjas navnPercent`[!is.na(sign$M1vsC)]
sign$navn[is.na(sign$navn)] <- sign$Cluster[is.na(sign$navn)]
sign2 <- sign[, c(3,2)]
colnames(sign2) <-  c("kort_navn", "navn")


ggMT1vsCp1 <- cytofSignPlot(d = d, sign_cl = sign_cl_MT1vsC, signnavn = sign2, main =  "Moderate T1 vs control")
ggMT1vsCp1$gg








tiff(fs::path(utSti, "P1_cytof_sign_plot_MT1vsC.tiff"), width = 500, height = 700)
ggMT1vsCp1$gg
dev.off()



d2 <- ggMT1vsCp1$d
d2 <- d2[!is.na(d2$kort_navn),]
d2$col <- "red"
d2$col[d2$diff == "DOWN"] <- "blue"
d2 <- d2[order(d2$IRRmotsatt),]
d2$stjerner <- lapply(d2$adj_p, pVerdiStjerne)
d2$lab2 <- paste0(d2$kort_navn, " (", d2$stjerner, ")")


tiff(fs::path(utSti, "log2P1_IRR_linestack_MT1_control.tiff"), width = 360, height = 700)
par(mar = c(5, 4, 4, 2) + 0.1)
linestack(log2(d2$IRRmotsatt), d2$kort_navn, axis = T, air = 1.1, col = d2$col, main = "Moderate T1 vs Control", ylim = c(-8,8))
dev.off()

d$lab2 <- paste0(d$delabel, " (", d$stjerner, ")")

tiff(fs::path(utSti, "log2P1_IRR_linestack_MT1_control_med stjerne.tiff"), width = 400, height = 700)
par(mar = c(5, 4, 4, 2) + 0.1)
linestack(log2(d2$IRRmotsatt), d2$lab2, axis = T, air = 1.1, col = d2$col, main = "Moderate T1 vs controls", ylim = c(-8,8))
dev.off()


tiff(fs::path(utSti, "P1_IRR_linestack_MT1_control.tiff"), width = 360, height = 700)
linestack(log2(d2$IRRmotsatt), d2$kort_navn, axis = F, air = 1.1, col = d2$col, main = "Moderate T1 vs controls", ylim = c(-8,8)) 
axis(side = 2, at = (-8):8, label = round(2^((-8):8),3))
dev.off()


tiff(fs::path(utSti, "P1_IRR_linestack_MT1_control_med stjerne.tiff"), width = 420, height = 700)
linestack(log2(d2$IRRmotsatt), d2$lab2, axis = F, air = 1.1, col = d2$col, main = "Moderate T1 vs controls", ylim = c(-8,8)) 
axis(side = 2, at = (-8):8, label =round(2^((-8):8),3))
dev.off()






d <- read.csv2("F:/Forskningsprosjekter/PDB 2794 - Immune responses aga_/Forskningsfiler/JOBO/CyTOF/Analyse i R OUS/CleanUpGatingMarch2022/Analyse/Endelig/Cytof unsupperviced/Endelig_des2022/Panel 1/All/seed 1345/ST2vsMT2/IRRadjp_1345.csv")
d$IRRmotsatt <- 1/d$IRR

sign_cl_ST2vsMT2 <- sign$`Anjas navnPercent`[!is.na(sign$ST2vsMT2)]
sign$navn[is.na(sign$navn)] <- sign$Cluster[is.na(sign$navn)]
sign2 <- sign[, c(3,2)]
colnames(sign2) <-  c("kort_navn", "navn")


ggST2vsMT2p1 <- cytofSignPlot(d = d, sign_cl = sign_cl_ST2vsMT2, signnavn = sign2, main =  "Severe T2 vs moderate T2")#, ymax = 3)
ggST2vsMT2p1$gg
tiff(fs::path(utSti, "P1_cytof_sign_plot_ST2vsMT2.tiff"), width = 500, height = 700)
ggST2vsMT2p1$gg
dev.off()


d2 <- ggST2vsMT2p1$d
d2 <- d2[!is.na(d2$kort_navn),]
d2$col <- "red"
d2$col[d2$diff == "DOWN"] <- "blue"
d2 <- d2[order(d2$IRRmotsatt),]
d2$stjerner <- lapply(d2$adj_p, pVerdiStjerne)
d2$lab2 <- paste0(d2$kort_navn, " (", d2$stjerner, ")")


tiff(fs::path(utSti, "log2P1_IRR_linestack_ST2_MT2.tiff"), width = 300, height = 700)
par(mar = c(5, 4, 4, 2) + 0.1)
linestack(log2(d2$IRRmotsatt), d2$kort_navn, axis = T, air = 1.1, col = d2$col, main = "Severe T2 vs Moderate T2", ylim = c(-8,8))
dev.off()

d$lab2 <- paste0(d$delabel, " (", d$stjerner, ")")

tiff(fs::path(utSti, "log2P1_IRR_linestack_ST2_MT2_med stjerne.tiff"), width = 340, height = 700)
par(mar = c(5, 4, 4, 2) + 0.1)
linestack(log2(d2$IRRmotsatt), d2$lab2, axis = T, air = 1.1, col = d2$col, main = "Severe T2 vs Moderate T2", ylim = c(-8,8))
dev.off()


tiff(fs::path(utSti, "P1_IRR_linestack_ST2_MT2.tiff"), width = 300, height = 700)
linestack(log2(d2$IRRmotsatt), d2$kort_navn, axis = F, air = 1.1, col = d2$col, main = "Severe T2 vs Moderate T2", ylim = c(-8,8)) 
axis(side = 2, at = (-8):8, label = round(2^((-8):8),3))
dev.off()


tiff(fs::path(utSti, "P1_IRR_linestack_ST2_MT2_med stjerne.tiff"), width = 340, height = 700)
linestack(log2(d2$IRRmotsatt), d2$lab2, axis = F, air = 1.1, col = d2$col, main = "Severe T2 vs Moderate T2", ylim = c(-8,8)) 
axis(side = 2, at = (-8):8, label =round(2^((-8):8),3))
dev.off()


d <- read.csv2("F:/Forskningsprosjekter/PDB 2794 - Immune responses aga_/Forskningsfiler/JOBO/CyTOF/Analyse i R OUS/CleanUpGatingMarch2022/Analyse/Endelig/Cytof unsupperviced/Endelig_des2022/Panel 1/All/seed 1345/ST2vsCont/IRRadjp_1345.csv")
d$IRRmotsatt <- 1/d$IRR

sign_cl_ST2vsC <- sign$`Anjas navnPercent`[!is.na(sign$ST2vsC)]
sign$navn[is.na(sign$navn)] <- sign$Cluster[is.na(sign$navn)]
sign2 <- sign[, c(3,2)]
colnames(sign2) <-  c("kort_navn", "navn")


ggST2vsCp1 <- cytofSignPlot(d = d, sign_cl = sign_cl_ST2vsC, signnavn = sign2, main =  "Severe T2 vs control")#, ymax = 3)
ggST2vsCp1$gg


tiff(fs::path(utSti, "P1_cytof_sign_plot_ST2vsC.tiff"), width = 500, height = 700)
ggST2vsCp1$gg
dev.off()



d2 <- ggST2vsCp1$d
d2 <- d2[!is.na(d2$kort_navn),]
d2$col <- "red"
d2$col[d2$diff == "DOWN"] <- "blue"
d2 <- d2[order(d2$IRRmotsatt),]
d2$stjerner <- lapply(d2$adj_p, pVerdiStjerne)
d2$lab2 <- paste0(d2$kort_navn, " (", d2$stjerner, ")")


tiff(fs::path(utSti, "log2P1_IRR_linestack_ST2_control.tiff"), width = 360, height = 700)
par(mar = c(5, 4, 4, 2) + 0.1)
linestack(log2(d2$IRRmotsatt), d2$kort_navn, axis = T, air = 1.1, col = d2$col, main = "Severe T2 vs Control", ylim = c(-8,8))
dev.off()

d$lab2 <- paste0(d$delabel, " (", d$stjerner, ")")

tiff(fs::path(utSti, "log2P1_IRR_linestack_ST2_control_med stjerne.tiff"), width = 400, height = 700)
par(mar = c(5, 4, 4, 2) + 0.1)
linestack(log2(d2$IRRmotsatt), d2$lab2, axis = T, air = 1.1, col = d2$col, main = "Severe T2 vs controls", ylim = c(-8,8))
dev.off()


tiff(fs::path(utSti, "P1_IRR_linestack_ST2_control.tiff"), width = 360, height = 700)
linestack(log2(d2$IRRmotsatt), d2$kort_navn, axis = F, air = 1.1, col = d2$col, main = "Severe T2 vs controls", ylim = c(-8,8)) 
axis(side = 2, at = (-8):8, label = round(2^((-8):8),3))
dev.off()


tiff(fs::path(utSti, "P1_IRR_linestack_ST2_control_med stjerne.tiff"), width = 420, height = 700)
linestack(log2(d2$IRRmotsatt), d2$lab2, axis = F, air = 1.1, col = d2$col, main = "Severe T2 vs controls", ylim = c(-8,8)) 
axis(side = 2, at = (-8):8, label =round(2^((-8):8),3))
dev.off()




d <- read.csv2("F:/Forskningsprosjekter/PDB 2794 - Immune responses aga_/Forskningsfiler/JOBO/CyTOF/Analyse i R OUS/CleanUpGatingMarch2022/Analyse/Endelig/Cytof unsupperviced/Endelig_des2022/Panel 1/All/seed 1345/MT2vsCont/IRRadjp_1345.csv")
d$IRRmotsatt <- 1/d$IRR

sign_cl_MT2vsC <- sign$`Anjas navnPercent`[!is.na(sign$MT2vsC)]
sign$navn[is.na(sign$navn)] <- sign$Cluster[is.na(sign$navn)]
sign2 <- sign[, c(3,2)]
colnames(sign2) <-  c("kort_navn", "navn")


ggMT2vsCp1 <- cytofSignPlot(d = d, sign_cl = sign_cl_MT2vsC, signnavn = sign2, main =  "Moderate T2 vs control")#, ymax = 3)
ggMT2vsCp1$gg


tiff(fs::path(utSti, "P1_cytof_sign_plot_MT2vsC.tiff"), width = 500, height = 700)
ggMT2vsCp1$gg
dev.off()


d2 <- ggMT2vsCp1$d
d2 <- d2[!is.na(d2$kort_navn),]
d2$col <- "red"
d2$col[d2$diff == "DOWN"] <- "blue"
d2 <- d2[order(d2$IRRmotsatt),]
d2$stjerner <- lapply(d2$adj_p, pVerdiStjerne)
d2$lab2 <- paste0(d2$kort_navn, " (", d2$stjerner, ")")


tiff(fs::path(utSti, "log2P1_IRR_linestack_MT2_control.tiff"), width = 360, height = 700)
par(mar = c(5, 4, 4, 2) + 0.1)
linestack(log2(d2$IRRmotsatt), d2$kort_navn, axis = T, air = 1.1, col = d2$col, main = "Moderate T2 vs Control", ylim = c(-8,8))
dev.off()

d$lab2 <- paste0(d$delabel, " (", d$stjerner, ")")

tiff(fs::path(utSti, "log2P1_IRR_linestack_MT2_control_med stjerne.tiff"), width = 400, height = 700)
par(mar = c(5, 4, 4, 2) + 0.1)
linestack(log2(d2$IRRmotsatt), d2$lab2, axis = T, air = 1.1, col = d2$col, main = "Moderate T2 vs controls", ylim = c(-8,8))
dev.off()


tiff(fs::path(utSti, "P1_IRR_linestack_MT2_control.tiff"), width = 360, height = 700)
linestack(log2(d2$IRRmotsatt), d2$kort_navn, axis = F, air = 1.1, col = d2$col, main = "Moderate T2 vs controls", ylim = c(-8,8)) 
axis(side = 2, at = (-8):8, label = round(2^((-8):8),3))
dev.off()


tiff(fs::path(utSti, "P1_IRR_linestack_MT2_control_med stjerne.tiff"), width = 420, height = 700)
linestack(log2(d2$IRRmotsatt), d2$lab2, axis = F, air = 1.1, col = d2$col, main = "Moderate T2 vs controls", ylim = c(-8,8)) 
axis(side = 2, at = (-8):8, label =round(2^((-8):8),3))
dev.off()



library(ggVennDiagram)
x <- list("ST1vsC" = as.character(ggST1vsCp1$signcl[!is.na(ggST1vsCp1$signcl)]),
          "ST1vsMT1" = as.character(ggST1vsMT1p1$signcl[!is.na(ggST1vsMT1p1$signcl)]), 
          "MT1vsC" = as.character(ggMT1vsCp1$signcl[!is.na(ggMT1vsCp1$signcl)]))

gvennT1 <- ggVennDiagram(x) + scale_color_brewer(palette = "Paired") + scale_x_continuous(expand = expansion(mult = .2))
gvennT1
tiff(fs::path(utSti, "P1_cytof_venn_ST1_MT1_C.tiff"), width = 700, height = 700)
gvennT1
dev.off()



library(ggVennDiagram)
x <- list("ST2vsC" = as.character(ggST2vsCp1$signcl[!is.na(ggST2vsCp1$signcl)]),
          "ST2vsMT2" = as.character(ggST2vsMT2p1$signcl[!is.na(ggST2vsMT2p1$signcl)]), 
          "MT2vsC" = as.character(ggMT2vsCp1$signcl[!is.na(ggMT2vsCp1$signcl)]))

gvennT2 <- ggVennDiagram(x) + scale_color_brewer(palette = "Paired") + scale_x_continuous(expand = expansion(mult = .2))
gvennT2
tiff(fs::path(utSti, "P1_cytof_venn_ST2_MT2_C.tiff"), width = 700, height = 700)
gvennT2
dev.off()



tiff(fs::path(utSti, "Alle_P1_Cytof_venn_ST1_MT1_C.tiff"), width = 1500, height = 700)
ggarrange(ggST1vsMT1p1$gg, ggST1vsCp1$gg, ggMT1vsCp1$gg, gvennT1 , ncol = 4, nrow = 1)
dev.off()



tiff(fs::path(utSti, "Alle_P1_Cytof_ST1_MT1_C.tiff"), width = 1500, height = 700)
ggarrange(ggST1vsMT1p1$gg, ggST1vsCp1$gg, ggMT1vsCp1$gg , ncol = 3, nrow = 1)
dev.off()



tiff(fs::path(utSti, "Alle_P1_Cytof_venn_ST2_MT2_C.tiff"), width = 1500, height = 700)
ggarrange(ggST2vsMT2p1$gg, ggST2vsCp1$gg, ggMT2vsCp1$gg, gvennT2 , ncol = 4, nrow = 1)
dev.off()


tiff(fs::path(utSti, "Alle_P1_Cytof_ST2_MT2_C.tiff"), width = 1500, height = 700)
ggarrange(ggST2vsMT2p1$gg, ggST2vsCp1$gg, ggMT2vsCp1$gg , ncol = 3, nrow = 1)
dev.off()



#panel 2


d <- read.csv2("F:/Forskningsprosjekter/PDB 2794 - Immune responses aga_/Forskningsfiler/JOBO/CyTOF/Analyse i R OUS/CleanUpGatingMarch2022/Analyse/Endelig/Cytof unsupperviced/Endelig_des2022/Panel 2/All/seed 2234/IRRadjp_2234.csv")
d$IRRmotsatt <- 1/d$IRR
sign <- as.data.frame(readxl::read_excel("F:/Forskningsprosjekter/PDB 2794 - Immune responses aga_/Forskningsfiler/JOBO/CyTOF/Analyse i R OUS/CleanUpGatingMarch2022/Analyse/Endelig/Cytof unsupperviced/Endelig_des2022/Figurer/15122022 PANEL2_ST1vsMT1VsCST2vsMT2VsC_FINAL.xlsx"))
sign <- sign[!is.na(sign$Cluster),]
sign <- sign[sign$`ta med` == 1,]
sign$`Anjas navnPercent` <- gsub("All_", "",sign$`Anjas navnPercent`)
sign_cl_ST1vsMT1 <- sign$`Anjas navnPercent`[!is.na(sign$S1vsM1)]
sign2 <- sign[, c(3,2)]
colnames(sign2) <-  c("kort_navn", "navn")


ggST1vsMT1p2 <- cytofSignPlot(d = d, sign_cl = sign_cl_ST1vsMT1, signnavn = sign2, main =  "Severe T1 vs moderate T1")#, ymax = 7)
ggST1vsMT1p2$gg
tiff(fs::path(utSti, "P2_cytof_sign_plot_ST1vsMT1.tiff"), width = 500, height = 700)
ggST1vsMT1p2$gg
dev.off()


d2 <- ggST1vsMT1p2$d
d2 <- d2[!is.na(d2$kort_navn),]
d2$col <- "red"
d2$col[d2$diff == "DOWN"] <- "blue"
d2 <- d2[order(d2$IRRmotsatt),]
d2$stjerner <- lapply(d2$adj_p, pVerdiStjerne)
d2$lab2 <- paste0(d2$kort_navn, " (", d2$stjerner, ")")


tiff(fs::path(utSti, "log2P2_IRR_linestack_ST1_MT1.tiff"), width = 300, height = 700)
par(mar = c(5, 4, 4, 2) + 0.1)
linestack(log2(d2$IRRmotsatt), d2$kort_navn, axis = T, air = 1.1, col = d2$col, main = "Severe T1 vs Moderate T1", ylim = c(-8,8))
dev.off()

d$lab2 <- paste0(d$delabel, " (", d$stjerner, ")")

tiff(fs::path(utSti, "log2P2_IRR_linestack_ST1_MT1_med stjerne.tiff"), width = 340, height = 700)
par(mar = c(5, 4, 4, 2) + 0.1)
linestack(log2(d2$IRRmotsatt), d2$lab2, axis = T, air = 1.1, col = d2$col, main = "Severe T1 vs Moderate T1", ylim = c(-8,8))
dev.off()


tiff(fs::path(utSti, "P2_IRR_linestack_ST1_MT1.tiff"), width = 300, height = 700)
linestack(log2(d2$IRRmotsatt), d2$kort_navn, axis = F, air = 1.1, col = d2$col, main = "Severe T1 vs Moderate T1", ylim = c(-8,8)) 
axis(side = 2, at = (-8):8, label = round(2^((-8):8),3))
dev.off()


tiff(fs::path(utSti, "P2_IRR_linestack_ST1_MT1_med stjerne.tiff"), width = 340, height = 700)
linestack(log2(d2$IRRmotsatt), d2$lab2, axis = F, air = 1.1, col = d2$col, main = "Severe T1 vs Moderate T1", ylim = c(-8,8)) 
axis(side = 2, at = (-8):8, label =round(2^((-8):8),3))
dev.off()



d <- read.csv2("F:/Forskningsprosjekter/PDB 2794 - Immune responses aga_/Forskningsfiler/JOBO/CyTOF/Analyse i R OUS/CleanUpGatingMarch2022/Analyse/Endelig/Cytof unsupperviced/Endelig_des2022/Panel 2/All/seed 2234/ST1vsCont/IRRadjp_2234.csv")
d$IRRmotsatt <- 1/d$IRR

sign_cl_ST1vsC <- sign$`Anjas navnPercent`[!is.na(sign$S1vsC)]
sign2 <- sign[, c(3,2)]
colnames(sign2) <-  c("kort_navn", "navn")


ggST1vsCp2 <- cytofSignPlot(d = d, sign_cl = sign_cl_ST1vsC, signnavn = sign2, main =  "Severe T1 vs control")#, ymax = 7)
ggST1vsCp2$gg


tiff(fs::path(utSti, "P2_cytof_sign_plot_ST1vsC.tiff"), width = 500, height = 700)
ggST1vsCp2$gg
dev.off()


d2 <- ggST1vsCp2$d
d2 <- d2[!is.na(d2$kort_navn),]
d2$col <- "red"
d2$col[d2$diff == "DOWN"] <- "blue"
d2 <- d2[order(d2$IRRmotsatt),]
d2$stjerner <- lapply(d2$adj_p, pVerdiStjerne)
d2$lab2 <- paste0(d2$kort_navn, " (", d2$stjerner, ")")


tiff(fs::path(utSti, "log2P2_IRR_linestack_ST1_control.tiff"), width = 360, height = 700)
par(mar = c(5, 4, 4, 2) + 0.1)
linestack(log2(d2$IRRmotsatt), d2$kort_navn, axis = T, air = 1.1, col = d2$col, main = "Severe T1 vs Control", ylim = c(-8,8))
dev.off()

d$lab2 <- paste0(d$delabel, " (", d$stjerner, ")")

tiff(fs::path(utSti, "log2P2_IRR_linestack_ST1_control_med stjerne.tiff"), width = 400, height = 700)
par(mar = c(5, 4, 4, 2) + 0.1)
linestack(log2(d2$IRRmotsatt), d2$lab2, axis = T, air = 1.1, col = d2$col, main = "Severe T1 vs controls", ylim = c(-8,8))
dev.off()


tiff(fs::path(utSti, "P2_IRR_linestack_ST1_control.tiff"), width = 360, height = 700)
linestack(log2(d2$IRRmotsatt), d2$kort_navn, axis = F, air = 1.1, col = d2$col, main = "Severe T1 vs controls", ylim = c(-8,8)) 
axis(side = 2, at = (-8):8, label = round(2^((-8):8),3))
dev.off()


tiff(fs::path(utSti, "P2_IRR_linestack_ST1_control_med stjerne.tiff"), width = 420, height = 700)
linestack(log2(d2$IRRmotsatt), d2$lab2, axis = F, air = 1.1, col = d2$col, main = "Severe T1 vs controls", ylim = c(-8,8)) 
axis(side = 2, at = (-8):8, label =round(2^((-8):8),3))
dev.off()


d <- read.csv2("F:/Forskningsprosjekter/PDB 2794 - Immune responses aga_/Forskningsfiler/JOBO/CyTOF/Analyse i R OUS/CleanUpGatingMarch2022/Analyse/Endelig/Cytof unsupperviced/Endelig_des2022/Panel 2/All/seed 2234/MT1vsCont/IRRadjp_2234.csv")
d$IRRmotsatt <- 1/d$IRR

sign_cl_MT1vsC <- sign$`Anjas navnPercent`[!is.na(sign$M1vsC)]
sign2 <- sign[, c(3,2)]
colnames(sign2) <-  c("kort_navn", "navn")


ggMT1vsCp2 <- cytofSignPlot(d = d, sign_cl = sign_cl_MT1vsC, signnavn = sign2, main =  "Moderate T1 vs control")#, ymax = 7)
ggMT1vsCp2$gg


tiff(fs::path(utSti, "P2_cytof_sign_plot_MT1vsC.tiff"), width = 500, height = 700)
ggMT1vsCp2$gg
dev.off()


d2 <- ggMT1vsCp2$d
d2 <- d2[!is.na(d2$kort_navn),]
d2$col <- "red"
d2$col[d2$diff == "DOWN"] <- "blue"
d2 <- d2[order(d2$IRRmotsatt),]
d2$stjerner <- lapply(d2$adj_p, pVerdiStjerne)
d2$lab2 <- paste0(d2$kort_navn, " (", d2$stjerner, ")")


tiff(fs::path(utSti, "log2P2_IRR_linestack_MT1_control.tiff"), width = 360, height = 700)
par(mar = c(5, 4, 4, 2) + 0.1)
linestack(log2(d2$IRRmotsatt), d2$kort_navn, axis = T, air = 1.1, col = d2$col, main = "Moderate T1 vs Control", ylim = c(-8,8))
dev.off()

d$lab2 <- paste0(d$delabel, " (", d$stjerner, ")")

tiff(fs::path(utSti, "log2P2_IRR_linestack_MT1_control_med stjerne.tiff"), width = 400, height = 700)
par(mar = c(5, 4, 4, 2) + 0.1)
linestack(log2(d2$IRRmotsatt), d2$lab2, axis = T, air = 1.1, col = d2$col, main = "Moderate T1 vs controls", ylim = c(-8,8))
dev.off()


tiff(fs::path(utSti, "P2_IRR_linestack_MT1_control.tiff"), width = 360, height = 700)
linestack(log2(d2$IRRmotsatt), d2$kort_navn, axis = F, air = 1.1, col = d2$col, main = "Moderate T1 vs controls", ylim = c(-8,8)) 
axis(side = 2, at = (-8):8, label = round(2^((-8):8),3))
dev.off()


tiff(fs::path(utSti, "P2_IRR_linestack_MT1_control_med stjerne.tiff"), width = 420, height = 700)
linestack(log2(d2$IRRmotsatt), d2$lab2, axis = F, air = 1.1, col = d2$col, main = "Moderate T1 vs controls", ylim = c(-8,8)) 
axis(side = 2, at = (-8):8, label =round(2^((-8):8),3))
dev.off()





d <- read.csv2("F:/Forskningsprosjekter/PDB 2794 - Immune responses aga_/Forskningsfiler/JOBO/CyTOF/Analyse i R OUS/CleanUpGatingMarch2022/Analyse/Endelig/Cytof unsupperviced/Endelig_des2022/Panel 2/All/seed 2234/ST2vsMT2/IRRadjp_2234.csv")
d$IRRmotsatt <- 1/d$IRR

sign_cl_ST2vsMT2 <- sign$`Anjas navnPercent`[!is.na(sign$ST2vsMT2)]
sign2 <- sign[, c(3,2)]
colnames(sign2) <-  c("kort_navn", "navn")


ggST2vsMT2p2 <- cytofSignPlot(d = d, sign_cl = sign_cl_ST2vsMT2, signnavn = sign2, main =  "Severe T2 vs moderate T2")#, ymax = 3)
ggST2vsMT2p2$gg
tiff(fs::path(utSti, "P2_cytof_sign_plot_ST2vsMT2.tiff"), width = 500, height = 700)
ggST2vsMT2p2$gg
dev.off()



d2 <- ggST2vsMT2p2$d
d2 <- d2[!is.na(d2$kort_navn),]
d2$col <- "red"
d2$col[d2$diff == "DOWN"] <- "blue"
d2 <- d2[order(d2$IRRmotsatt),]
d2$stjerner <- lapply(d2$adj_p, pVerdiStjerne)
d2$lab2 <- paste0(d2$kort_navn, " (", d2$stjerner, ")")


tiff(fs::path(utSti, "log2P2_IRR_linestack_ST2_MT2.tiff"), width = 300, height = 700)
par(mar = c(5, 4, 4, 2) + 0.1)
linestack(log2(d2$IRRmotsatt), d2$kort_navn, axis = T, air = 1.1, col = d2$col, main = "Severe T2 vs Moderate T2", ylim = c(-8,8))
dev.off()

d$lab2 <- paste0(d$delabel, " (", d$stjerner, ")")

tiff(fs::path(utSti, "log2P2_IRR_linestack_ST2_MT2_med stjerne.tiff"), width = 340, height = 700)
par(mar = c(5, 4, 4, 2) + 0.1)
linestack(log2(d2$IRRmotsatt), d2$lab2, axis = T, air = 1.1, col = d2$col, main = "Severe T2 vs Moderate T2", ylim = c(-8,8))
dev.off()


tiff(fs::path(utSti, "P2_IRR_linestack_ST2_MT2.tiff"), width = 300, height = 700)
linestack(log2(d2$IRRmotsatt), d2$kort_navn, axis = F, air = 1.1, col = d2$col, main = "Severe T2 vs Moderate T2", ylim = c(-8,8)) 
axis(side = 2, at = (-8):8, label = round(2^((-8):8),3))
dev.off()


tiff(fs::path(utSti, "P2_IRR_linestack_ST2_MT2_med stjerne.tiff"), width = 340, height = 700)
linestack(log2(d2$IRRmotsatt), d2$lab2, axis = F, air = 1.1, col = d2$col, main = "Severe T2 vs Moderate T2", ylim = c(-8,8)) 
axis(side = 2, at = (-8):8, label =round(2^((-8):8),3))
dev.off()


d <- read.csv2("F:/Forskningsprosjekter/PDB 2794 - Immune responses aga_/Forskningsfiler/JOBO/CyTOF/Analyse i R OUS/CleanUpGatingMarch2022/Analyse/Endelig/Cytof unsupperviced/Endelig_des2022/Panel 2/All/seed 2234/ST2vsCont/IRRadjp_2234.csv")
d$IRRmotsatt <- 1/d$IRR

sign_cl_ST2vsC <- sign$`Anjas navnPercent`[!is.na(sign$ST2vsC)]
sign2 <- sign[, c(3,2)]
colnames(sign2) <-  c("kort_navn", "navn")


ggST2vsCp2 <- cytofSignPlot(d = d, sign_cl = sign_cl_ST2vsC, signnavn = sign2, main =  "Severe T2 vs control")#, ymax = 3)
ggST2vsCp2$gg


tiff(fs::path(utSti, "P2_cytof_sign_plot_ST2vsC.tiff"), width = 500, height = 700)
ggST2vsCp2$gg
dev.off()




d2 <- ggST2vsCp2$d
d2 <- d2[!is.na(d2$kort_navn),]
d2$col <- "red"
d2$col[d2$diff == "DOWN"] <- "blue"
d2 <- d2[order(d2$IRRmotsatt),]
d2$stjerner <- lapply(d2$adj_p, pVerdiStjerne)
d2$lab2 <- paste0(d2$kort_navn, " (", d2$stjerner, ")")


tiff(fs::path(utSti, "log2P2_IRR_linestack_ST2_control.tiff"), width = 360, height = 700)
par(mar = c(5, 4, 4, 2) + 0.1)
linestack(log2(d2$IRRmotsatt), d2$kort_navn, axis = T, air = 1.1, col = d2$col, main = "Severe T2 vs Control", ylim = c(-8,8))
dev.off()

d$lab2 <- paste0(d$delabel, " (", d$stjerner, ")")

tiff(fs::path(utSti, "log2P2_IRR_linestack_ST2_control_med stjerne.tiff"), width = 400, height = 700)
par(mar = c(5, 4, 4, 2) + 0.1)
linestack(log2(d2$IRRmotsatt), d2$lab2, axis = T, air = 1.1, col = d2$col, main = "Severe T2 vs controls", ylim = c(-8,8))
dev.off()


tiff(fs::path(utSti, "P2_IRR_linestack_ST2_control.tiff"), width = 360, height = 700)
linestack(log2(d2$IRRmotsatt), d2$kort_navn, axis = F, air = 1.1, col = d2$col, main = "Severe T2 vs controls", ylim = c(-8,8)) 
axis(side = 2, at = (-8):8, label = round(2^((-8):8),3))
dev.off()


tiff(fs::path(utSti, "P2_IRR_linestack_ST2_control_med stjerne.tiff"), width = 420, height = 700)
linestack(log2(d2$IRRmotsatt), d2$lab2, axis = F, air = 1.1, col = d2$col, main = "Severe T2 vs controls", ylim = c(-8,8)) 
axis(side = 2, at = (-8):8, label =round(2^((-8):8),3))
dev.off()




d <- read.csv2("F:/Forskningsprosjekter/PDB 2794 - Immune responses aga_/Forskningsfiler/JOBO/CyTOF/Analyse i R OUS/CleanUpGatingMarch2022/Analyse/Endelig/Cytof unsupperviced/Endelig_des2022/Panel 2/All/seed 2234/MT2vsCont/IRRadjp_2234.csv")
d$IRRmotsatt <- 1/d$IRR

sign_cl_MT2vsC <- sign$`Anjas navnPercent`[!is.na(sign$MT2vsC)]
sign2 <- sign[, c(3,2)]
colnames(sign2) <-  c("kort_navn", "navn")


ggMT2vsCp2 <- cytofSignPlot(d = d, sign_cl = sign_cl_MT2vsC, signnavn = sign2, main =  "Moderate T2 vs control")#, ymax = 3)
ggMT2vsCp2$gg


tiff(fs::path(utSti, "P2_cytof_sign_plot_MT2vsC.tiff"), width = 500, height = 700)
ggMT2vsCp2$gg
dev.off()



d2 <- ggMT2vsCp2$d
d2 <- d2[!is.na(d2$kort_navn),]
d2$col <- "red"
d2$col[d2$diff == "DOWN"] <- "blue"
d2 <- d2[order(d2$IRRmotsatt),]
d2$stjerner <- lapply(d2$adj_p, pVerdiStjerne)
d2$lab2 <- paste0(d2$kort_navn, " (", d2$stjerner, ")")


tiff(fs::path(utSti, "log2P2_IRR_linestack_MT2_control.tiff"), width = 360, height = 700)
par(mar = c(5, 4, 4, 2) + 0.1)
linestack(log2(d2$IRRmotsatt), d2$kort_navn, axis = T, air = 1.1, col = d2$col, main = "Moderate T2 vs Control", ylim = c(-8,8))
dev.off()

d$lab2 <- paste0(d$delabel, " (", d$stjerner, ")")

tiff(fs::path(utSti, "log2P2_IRR_linestack_MT2_control_med stjerne.tiff"), width = 400, height = 700)
par(mar = c(5, 4, 4, 2) + 0.1)
linestack(log2(d2$IRRmotsatt), d2$lab2, axis = T, air = 1.1, col = d2$col, main = "Moderate T2 vs controls", ylim = c(-8,8))
dev.off()


tiff(fs::path(utSti, "P2_IRR_linestack_MT2_control.tiff"), width = 360, height = 700)
linestack(log2(d2$IRRmotsatt), d2$kort_navn, axis = F, air = 1.1, col = d2$col, main = "Moderate T2 vs controls", ylim = c(-8,8)) 
axis(side = 2, at = (-8):8, label = round(2^((-8):8),3))
dev.off()


tiff(fs::path(utSti, "P2_IRR_linestack_MT2_control_med stjerne.tiff"), width = 420, height = 700)
linestack(log2(d2$IRRmotsatt), d2$lab2, axis = F, air = 1.1, col = d2$col, main = "Moderate T2 vs controls", ylim = c(-8,8)) 
axis(side = 2, at = (-8):8, label =round(2^((-8):8),3))
dev.off()




library(ggVennDiagram)
x <- list("ST1vsC" = as.character(ggST1vsCp2$signcl[!is.na(ggST1vsCp2$signcl)]),
          "ST1vsMT1" = as.character(ggST1vsMT1p2$signcl[!is.na(ggST1vsMT1p2$signcl)]), 
          "MT1vsC" = as.character(ggMT1vsCp2$signcl[!is.na(ggMT1vsCp2$signcl)]))

gvennT1 <- ggVennDiagram(x) + scale_color_brewer(palette = "Paired") + scale_x_continuous(expand = expansion(mult = .2))
gvennT1
tiff(fs::path(utSti, "P2_cytof_venn_ST1_MT1_C.tiff"), width = 700, height = 700)
gvennT1
dev.off()



library(ggVennDiagram)
x <- list("ST2vsC" = as.character(ggST2vsCp2$signcl[!is.na(ggST2vsCp2$signcl)]),
          "ST2vsMT2" = as.character(ggST2vsMT2p2$signcl[!is.na(ggST2vsMT2p2$signcl)]), 
          "MT2vsC" = as.character(ggMT2vsCp2$signcl[!is.na(ggMT2vsCp2$signcl)]))

gvennT2 <- ggVennDiagram(x) + scale_color_brewer(palette = "Paired") + scale_x_continuous(expand = expansion(mult = .2))
gvennT2
tiff(fs::path(utSti, "P2_cytof_venn_ST2_MT2_C.tiff"), width = 700, height = 700)
gvennT2
dev.off()



tiff(fs::path(utSti, "Alle_P2_Cytof_venn_ST1_MT1_C.tiff"), width = 1500, height = 700)
ggarrange(ggST1vsMT1p2$gg, ggST1vsCp2$gg, ggMT1vsCp2$gg , ncol = 3, nrow = 1)
dev.off()


tiff(fs::path(utSti, "Alle_P2_Cytof_ST1_MT1_C.tiff"), width = 1500, height = 700)
ggarrange(ggST1vsMT1p2$gg, ggST1vsCp2$gg, ggMT1vsCp2$gg , ncol = 3, nrow = 1)
dev.off()



tiff(fs::path(utSti, "Alle_P2_Cytof_ST2_MT2_C.tiff"), width = 1500, height = 700)
ggarrange(ggST2vsMT2p2$gg, ggST2vsCp2$gg, ggMT2vsCp2$gg , ncol = 3, nrow = 1)
dev.off()



tiff(fs::path(utSti, "Alle_P2_Cytof_venn_ST2_MT2_C.tiff"), width = 1500, height = 700)
ggarrange(ggST2vsMT2p2$gg, ggST2vsCp2$gg, ggMT2vsCp2$gg, gvennT2 , ncol = 4, nrow = 1)
dev.off()


