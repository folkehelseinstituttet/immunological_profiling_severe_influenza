
analyseSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "gen")
utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "Figurer")
utSti2 <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "Cytof unsupperviced", "Endelig_des2022", "Figurer")

d <- read.csv2(fs::path(analyseSti, "dGener.csv"))


fra <- which(colnames(d) == "AIRE")
til <- which(colnames(d) == "ZNF532")

for(i in fra:til){
  x <- d[,i]
  max_value <- mean(x) + 2*sd(x)
  d[d[,i] > max_value, i] <- max_value
}

sign_funn_gen <- d[,c("ID_tid",  "CD4")]
colnames(sign_funn_gen)[2] <- "CD4, Gene"


#panel 2

scriptPath_UnsupAnalysis <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "Cytof unsupperviced")
utSti <- fs::path(scriptPath_UnsupAnalysis,  "resultat_panel2_cor095", "figurer")

colorlegend <- c("Severe T1" = "red", "Moderate T1" = "blue", "Severe T2" = "pink", "Moderate T2" = "light blue", "Control" = "light green")

params <- list()
params$panel <- "Panel 2"
params$seed <- 2234#nb mC% endres vil man vil gjC8re et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$utSti <- fs::path(scriptPath_UnsupAnalysis, "Endelig_des2022", params$panel, params$ext_name, paste("seed", params$seed))
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c("Severe T1", "Moderate T1", "Severe T2", "Moderate T2", "Control" )


d <- read.csv2(fs::path(params$utSti, paste0("unike_klustre", params$seed, ".csv")))
d$statusTid <- factor(d$statusTid, levels = params$nivaa)


# d$celler <- "Gr 0 No difference"
# d$celler[d$`kluster 30` == 10] <- "Gr 4: CD8 EM IFNg"  # "gr 3 (CD 8 EM)"
# d$celler[d$`kluster 20` == 8] <- "Gr 2: CD4 CM CD137"  #"gr 5 (CD 4 CM)"
# d$celler[d$`kluster 20` == 4] <-  "Gr 6: NK GranzB & TIM3" #"gr 9 (NK)"
# d$celler[d$`kluster 50` == 19] <-  "Gr 5: CD8 TEMRA GranzB" #"gr 10 (CD 8)"
# d$celler[d$`kluster 20` == 10] <- "Gr 1: CD4 CM"  #"gr 11 (CD 4 CM)"
# d$celler[d$`kluster 40` == 31] <- "Gr 3: CD4 CM PD-L2"  #"gr 8 (CD 4 CM)"



clusters <- colnames(d)[grepl("kluster", colnames(d))]
length(clusters)
# med <- c("seed_2234_k_30_kluster_10", "seed_2234_k_20_kluster_8", "seed_2234_k_40_kluster_31", 
#              "seed_2234_k_20_kluster_4", "seed_2234_k_50_kluster_19", "seed_2234_k_20_kluster_10")
# groups <- c("P2 Cl 4 (CD8 EM IFNg)", "P2 Cl 2 (CD4 CM CD137)", "P2 Cl 3 (CD4 CM PD-L2)", 
#             "P2 Cl 6 (NK GranzB & TIM3)", "P2 Cl 5 (CD8 TEMRA GranzB)", "P2 Cl 1 (CD4 CM)")


unsup_p2_tid1 <- readxl::read_excel(fs::path(scriptPath_UnsupAnalysis, "Endelig_des2022", "Figurer", "PANEL2_ST1vsMT1_FINAL.xlsx"))

med <- gsub("All_", "", unsup_p2_tid1$`Anjas navnPercent`)
groups <- paste0(unsup_p2_tid1$Navn_boxplot, ", Cytof Stim")

sign_P2 <- d[, c(med, "kort_filnavn", "status", "age", "sex", "statusTid", "tid") ]
sign_P2[,med] <- log(sign_P2[,med]/matrix(rep(d$antall_cells_totalt, length(med)), ncol = length(med), byrow = F) + 0.0001)
colnames(sign_P2) <- c(groups, "kort_filnavn", "status", "age", "sex", "statusTid", "tid")


d_alle <- merge(sign_P2, sign_funn_gen, by.x = "kort_filnavn", by.y = "ID_tid")



#panel 1

scriptPath_UnsupAnalysis <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "Cytof unsupperviced")
utSti <- fs::path(scriptPath_UnsupAnalysis,  "resultat_panel1_cor095", "figurer")

colorlegend <- c("Severe T1" = "red", "Moderate T1" = "blue", "Severe T2" = "pink", "Moderate T2" = "light blue", "Control" = "light green")

params <- list()
params$panel <- "Panel 1"
params$seed <- 1345#nb mC% endres vil man vil gjC8re et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$utSti <- fs::path(scriptPath_UnsupAnalysis, "Endelig_des2022", params$panel, params$ext_name, paste("seed", params$seed))
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c("Severe T1", "Moderate T1", "Severe T2", "Moderate T2", "Control" )


d <- read.csv2(fs::path(params$utSti, paste0("unike_klustre", params$seed, ".csv")))
d$statusTid <- factor(d$statusTid, levels = params$nivaa)


clusters <- colnames(d)[grepl("kluster", colnames(d))]
length(clusters)


# d$celler[d$`kluster 10` == 1] <- "Gr 6 NK CD38 & (NKG2A)" #"gr 1 (NK)"
# d$celler[d$`kluster 30` == 26] <- "Gr 2 CD8 CM CXCR3 & TIGIT & PD-1" # "gr 7 (CD8 CM)"
# d$celler[d$`kluster 60` == 51] <- "Gr 1 CD4 CM CXCR5 & ICOS"  # "gr 9 (CD4 CM)"
# d$celler[d$`kluster 20` == 7] <-  "Gr 7 Mo CD169 CD85j & CD95" #"gr 14 (mo CD169)"
# d$celler[d$`kluster 10` == 4] <- "Gr 8 DC myeloid CD85j" #"gr 15 (DC myeloid)"
# d$celler[d$`kluster 30` == 4] <- "Gr 9 DC plasmacytoid CD85j" #"gr 17 (DC plasmacytoid)"
# d$celler[d$`kluster 50` == 18] <- "Gr 4 B CXCR5 CD85j" #"gr 18 (B)"
# d$celler[d$`kluster 20` == 14] <- "Gr 3 MAIT (PD-1) & CD95" #"gr 23 (MAIT)"
# d$celler[d$`kluster 50` == 17] <- "Gr 5 plasmacells CD95 & CD38" #"gr 27 (plasmacell)"



# med <- c("seed_1345_k_10_kluster_1", "seed_1345_k_30_kluster_26", "seed_1345_k_60_kluster_51", 
#              "seed_1345_k_20_kluster_7", "seed_1345_k_10_kluster_4", "seed_1345_k_30_kluster_4", 
#              "seed_1345_k_50_kluster_18", "seed_1345_k_20_kluster_14", "seed_1345_k_50_kluster_17")
# 
# 
# groups <- c("P1 Cl 6 (NK CD38 & (NKG2A))", "P1 Cl 2 (CD8 CM CXCR3 & TIGIT & PD-1)", "P1 Cl 1 (CD4 CM CXCR5 & ICOS)", 
#             "P1 Cl 7 (Mo CD169 CD85j & CD95)", "P1 Cl 8 (DC myeloid CD85j)", "P1 Cl 9 (DC plasmacytoid CD85j)",
#             "P1 Cl 4 (B CXCR5 CD85j)", "P1 Cl 3 (MAIT (PD-1) & CD95)" , "P1 Cl 5 (plasmacells CD95 & CD38)")



unsup_p1_tid1 <- readxl::read_excel(fs::path(scriptPath_UnsupAnalysis, "Endelig_des2022", "Figurer", "151222_Panel1_ST1vsMT1_FINAL.xlsx"))

med <- unsup_p1_tid1$`Anjas navnPercent`
groups <- paste0(unsup_p1_tid1$Navn_boxplot, ", Cytof Unstim")

sign_p1 <- d[, c(med, "kort_filnavn") ]
sign_p1[,med] <- log(sign_p1[,med]/matrix(rep(d$antall_cells_totalt, length(med)), ncol = length(med), byrow = F) + 0.0001)
colnames(sign_p1) <-  c(groups, "kort_filnavn")

#colnames(sign_p1) <- gsub("Cl", "P1 Cl", colnames(sign_p1))


d_alle <- merge(sign_p1, d_alle, by = "kort_filnavn")



# manuel 


scriptPath_manuel <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig")
params$utSti <- fs::path(scriptPath_manuel,  "figUMAP")

d <- read.csv2(fs::path(params$utSti, paste0("cells_per_manuel_cluster_and_sample.csv")))

med <- c("CD4", "MAIT", "Mo", "NK", "NKT")
navn <- c("CD4, Cytof Unstim manuel", "MAIT, Cytof Unstim manuel", "Mo, Cytof Unstim manuel", "NK, Cytof Unstim manuel", "NKT, Cytof Unstim manuel")

sign_p1_man <- d[,c(med, "X")]
sign_p1_man[,med] <- log(sign_p1_man[,med]/matrix(rep(d$totalt_antall, length(med)), ncol = length(med), byrow = F) + 0.0001)
colnames(sign_p1_man) <- c(navn, "kort_filnavn")

d_alle <- merge(sign_p1_man, d_alle, by = "kort_filnavn")


#CLCX13

analyseSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "gen")


d <- read.csv2(fs::path(analyseSti, "elisa_mm_status"))
d <- d[,c("Pasientnr", "CXCL13")]
d$CXCL13 <- log(d$CXCL13)
colnames(d) <- c("Pasientnr", "CXCL13, plasma")

d_alle <- merge(d_alle, d, by.x = "kort_filnavn", by.y = "Pasientnr")





# d_alle <- d_alle[,c("P1 manuel: CD4",  "P1 manuel: MAIT",  "P1 manuel: Mo", "P1 manuel: NK", "P1 manuel: NKT",
#                     colnames(sign_p1)[grep("P1", colnames(sign_p1))], colnames(sign_P2)[grep("P2", colnames(sign_P2))],
#                     "Gene CD4", "CXCL13", 
#                     "kort_filnavn", "status", "age", "sex", "statusTid", "tid")]


rekkefolge <- c("CD4, Gene","CD4, Cytof Unstim manuel",  "Cl1 (CD4 Tfh CM), Cytof Unstim", "Cl1 (CD4 CM), Cytof Stim",
                "Cl2 (CD4 CM PD-L2), Cytof Stim", "Cl3 (CD4 CM CD137), Cytof Stim", 
                "Cl2 (CD8 CM), Cytof Unstim", "Cl4 (CD8 EM PD-1), Cytof Stim", "Cl5 (CD8 TEMRA GranzB), Cytof Stim", 
                "MAIT, Cytof Unstim manuel",  "Cl3 (MAIT CD4), Cytof Unstim", "Cl4 (MAIT CD8), Cytof Unstim",
                "Cl5 (B CXCR5 CD85j), Cytof Unstim", 
                "Cl6 (plasma cells), Cytof Unstim", 
                "NK, Cytof Unstim manuel", "Cl7 (NK), Cytof Unstim", "Cl6 (NK GranzB & IFNg), Cytof Stim", 
                "Cl7 (NK GranzB), Cytof Stim", "Cl8 (NK GranzB & TIM-3), Cytof Stim", 
                "NKT, Cytof Unstim manuel",
                "Mo, Cytof Unstim manuel", "Cl8 (Mo CD169 CD123neg), Cytof Unstim", "Cl9 (Mo CD169 CD123), Cytof Unstim",
                "Cl10 (DC myeloid CD16), Cytof Unstim", "Cl11 (DC plasmacytoid), Cytof Unstim", 
                "CXCL13, plasma")




d_alle <- d_alle[,c("CD4, Gene","CD4, Cytof Unstim manuel",  "Cl1 (CD4 Tfh CM), Cytof Unstim", "Cl1 (CD4 CM), Cytof Stim",
                    "Cl2 (CD4 CM PD-L2), Cytof Stim", "Cl3 (CD4 CM CD137), Cytof Stim", 
                    "Cl2 (CD8 CM), Cytof Unstim", "Cl4 (CD8 EM PD-1), Cytof Stim", "Cl5 (CD8 TEMRA GranzB), Cytof Stim", 
                    "MAIT, Cytof Unstim manuel",  "Cl3 (MAIT CD4), Cytof Unstim", "Cl4 (MAIT CD8), Cytof Unstim",
                    "Cl5 (B CXCR5 CD85j), Cytof Unstim", 
                    "Cl6 (plasma cells), Cytof Unstim", 
                    "NK, Cytof Unstim manuel", "Cl7 (NK), Cytof Unstim", "Cl6 (NK GranzB & IFNg), Cytof Stim", 
                    "Cl7 (NK GranzB), Cytof Stim", "Cl8 (NK GranzB & TIM-3), Cytof Stim", 
                    "NKT, Cytof Unstim manuel",
                    "Mo, Cytof Unstim manuel", "Cl8 (Mo CD169 CD123neg), Cytof Unstim", "Cl9 (Mo CD169 CD123), Cytof Unstim",
                    "Cl10 (DC myeloid CD16), Cytof Unstim", "Cl11 (DC plasmacytoid), Cytof Unstim", 
                     "CXCL13, plasma", 
                    "kort_filnavn", "status", "age", "sex", "statusTid", "tid")]

#rekkefC8lge   c("CD4", "CD8",  "gdT", "MAIT", "B", "NK", "NKT", "Mo", "DC")


  
fra <- 1
til <- which(colnames(d_alle) == "CXCL13, plasma")

d_alle_raa <- d_alle
rownames(d_alle) <- d_alle$kort_filnavn
d_alle_scaled <- d_alle
d_alle_scaled[, fra:til] <- apply(d_alle_scaled[, fra:til], 2, scale)

summary(d_alle_scaled)


library(ComplexHeatmap)

col = list(statusTid = c("Severe T1" = "#990000", "Moderate T1" = "#000099", "Severe T2" = "#FF6666", "Moderate T2" = "#6699FF", "Control" = "#00CC00"),  
           sex = c("F" = "purple", "M" = "green"))
  

# Create the heatmap annotation
ha <- HeatmapAnnotation(
  statusTid = d_alle_scaled$statusTid, sex = d_alle_scaled$sex,
  col = col
)


set.seed(2)
#pa = pam(t(dat1.gener), k = 3)
#Heatmap(t(dat1.gener), split = paste0("pam", pa$clustering),  top_annotation = ha)

Heatmap(t(d_alle_scaled[,fra:til]), top_annotation = ha) 



#heatmap kun T1

d_alle_T1 <- d_alle[d_alle$tid == "T1", ]



d_alle_T1_scaled <- d_alle_T1
d_alle_T1_scaled[, fra:til] <- apply(d_alle_T1_scaled[, fra:til], 2, scale)



col = list(statusTid = c("Severe T1" = "#990000", "Moderate T1" = "#000099"),  
           sex = c("F" = "purple", "M" = "green"))


# Create the heatmap annotation
ha <- HeatmapAnnotation(
  statusTid = d_alle_T1_scaled$statusTid, sex = d_alle_T1_scaled$sex,
  col = col
)


set.seed(2)
#pa = pam(t(dat1.gener), k = 3)
#Heatmap(t(dat1.gener), split = paste0("pam", pa$clustering),  top_annotation = ha)

Heatmap(t(d_alle_T1_scaled[,fra:til]), top_annotation = ha) 




library(corrplot)
library(RColorBrewer)
M <-cor(d_alle_scaled[,fra:til])

corrplot(M, type="upper",  order="hclust",
         col= colorRampPalette(c("blue","yellow", "red"))(10), 
         tl.col = "black") #brewer.pal(n=8, name="RdYlBu"))


corrplot(M[rekkefolge,rekkefolge], type="upper", order="original",
         col= colorRampPalette(c("blue","yellow", "red"))(10), 
         tl.col = "black") #brewer.pal(n=8, name="RdYlBu"))



utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "Figurer")
utSti2 <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "Cytof unsupperviced", "Endelig_des2022", "Figurer")

tiff(fs::path(utSti, "corrplot.tiff"), width = 600, height = 600)
corrplot(M[rekkefolge, rekkefolge], type="upper",  order="original",
         col= colorRampPalette(c("blue","yellow", "red"))(10), tl.col = "black") #brewer.pal(n=8, name="RdYlBu"))
dev.off()
write.csv2(M[rekkefolge, rekkefolge], fs::path(utSti, "corrplot.csv"))
                  
tiff(fs::path(utSti2, "corrplot.tiff"), width = 600, height = 600)
corrplot(M[rekkefolge, rekkefolge], type="upper",  order="original",
         col= colorRampPalette(c("blue","yellow", "red"))(10), tl.col = "black") #brewer.pal(n=8, name="RdYlBu"))
dev.off()
write.csv2(M[rekkefolge, rekkefolge], fs::path(utSti2, "corrplot.csv"))


M_ST1 <-cor(d_alle_scaled[d_alle_scaled$statusTid == "Severe T1",fra:til])
corrplot(M_ST1[rekkefolge, rekkefolge], type="upper",  order="original",
         col= colorRampPalette(c("blue","yellow", "red"))(10), 
         tl.col = "black") #brewer.pal(n=8, name="RdYlBu"))

tiff(fs::path(utSti, "corrplotST1.tiff"), width = 600, height = 600)
corrplot(M_ST1[rekkefolge, rekkefolge], type="upper",  order="original",
         col= colorRampPalette(c("blue","yellow", "red"))(10), tl.col = "black") #brewer.pal(n=8, name="RdYlBu"))
dev.off()
write.csv2(M_ST1[rekkefolge, rekkefolge], fs::path(utSti, "corrplotST1.csv"))

tiff(fs::path(utSti2, "corrplotST1.tiff"), width = 600, height = 600)
corrplot(M_ST1[rekkefolge, rekkefolge], type="upper",  order="original",
         col= colorRampPalette(c("blue","yellow", "red"))(10), tl.col = "black") #brewer.pal(n=8, name="RdYlBu"))
dev.off()
write.csv2(M_ST1[rekkefolge, rekkefolge], fs::path(utSti2, "corrplotST1.csv"))





M_MT1 <-cor(d_alle_scaled[d_alle_scaled$statusTid == "Moderate T1",fra:til])
corrplot(M_MT1[rekkefolge, rekkefolge], type="upper",  order="original",
         col= colorRampPalette(c("blue","yellow", "red"))(10), 
         tl.col = "black") #brewer.pal(n=8, name="RdYlBu"))

tiff(fs::path(utSti, "corrplotMT1.tiff"), width = 600, height = 600)
corrplot(M_MT1[rekkefolge, rekkefolge], type="upper",  order="original",
         col= colorRampPalette(c("blue","yellow", "red"))(10), tl.col = "black") #brewer.pal(n=8, name="RdYlBu"))
dev.off()
write.csv2(M_MT1[rekkefolge, rekkefolge], fs::path(utSti, "corrplotMT1.csv"))

tiff(fs::path(utSti2, "corrplotMT1.tiff"), width = 600, height = 600)
corrplot(M_MT1[rekkefolge, rekkefolge], type="upper",  order="original",
         col= colorRampPalette(c("blue","yellow", "red"))(10), tl.col = "black") #brewer.pal(n=8, name="RdYlBu"))
dev.off()
write.csv2(M_MT1[rekkefolge, rekkefolge], fs::path(utSti2, "corrplotMT1.csv"))





M_ST2 <-cor(d_alle_scaled[d_alle_scaled$statusTid == "Severe T2",fra:til])
corrplot(M_ST2[rekkefolge, rekkefolge], type="upper",  order="original",
         col= colorRampPalette(c("blue","yellow", "red"))(10), 
         tl.col = "black") #brewer.pal(n=8, name="RdYlBu"))

tiff(fs::path(utSti, "corrplotST2.tiff"), width = 600, height = 600)
corrplot(M_ST2[rekkefolge, rekkefolge], type="upper",  order="original",
         col= colorRampPalette(c("blue","yellow", "red"))(10), tl.col = "black") #brewer.pal(n=8, name="RdYlBu"))
dev.off()
write.csv2(M_ST2[rekkefolge, rekkefolge], fs::path(utSti, "corrplotST2.csv"))

tiff(fs::path(utSti2, "corrplotST2.tiff"), width = 600, height = 600)
corrplot(M_ST2[rekkefolge, rekkefolge], type="upper",  order="original",
         col= colorRampPalette(c("blue","yellow", "red"))(10), tl.col = "black") #brewer.pal(n=8, name="RdYlBu"))
dev.off()
write.csv2(M_ST2[rekkefolge, rekkefolge], fs::path(utSti2, "corrplotST2.csv"))





M_MT2 <-cor(d_alle_scaled[d_alle_scaled$statusTid == "Moderate T2",fra:til])
corrplot(M_MT2[rekkefolge, rekkefolge], type="upper",  order="original",
         col= colorRampPalette(c("blue","yellow", "red"))(10), 
         tl.col = "black") #brewer.pal(n=8, name="RdYlBu"))

tiff(fs::path(utSti, "corrplotMT2.tiff"), width = 600, height = 600)
corrplot(M_MT2[rekkefolge, rekkefolge], type="upper",  order="original",
         col= colorRampPalette(c("blue","yellow", "red"))(10), tl.col = "black") #brewer.pal(n=8, name="RdYlBu"))
dev.off()
write.csv2(M_MT2[rekkefolge, rekkefolge], fs::path(utSti, "corrplotMT2.csv"))

tiff(fs::path(utSti2, "corrplotMT2.tiff"), width = 600, height = 600)
corrplot(M_MT2[rekkefolge, rekkefolge], type="upper",  order="original",
         col= colorRampPalette(c("blue","yellow", "red"))(10), tl.col = "black") #brewer.pal(n=8, name="RdYlBu"))
dev.off()
write.csv2(M_MT2[rekkefolge, rekkefolge], fs::path(utSti2, "corrplotMT2.csv"))



M_C <-cor(d_alle_scaled[d_alle_scaled$statusTid == "Control",fra:til])
corrplot(M_C[rekkefolge, rekkefolge], type="upper",  order="original",
         col= colorRampPalette(c("blue","yellow", "red"))(10), 
         tl.col = "black") #brewer.pal(n=8, name="RdYlBu"))

tiff(fs::path(utSti, "corrplotC.tiff"), width = 600, height = 600)
corrplot(M_C[rekkefolge, rekkefolge], type="upper",  order="original",
         col= colorRampPalette(c("blue","yellow", "red"))(10), tl.col = "black") #brewer.pal(n=8, name="RdYlBu"))
dev.off()
write.csv2(M_C[rekkefolge, rekkefolge], fs::path(utSti, "corrplotC.csv"))

tiff(fs::path(utSti2, "corrplotC.tiff"), width = 600, height = 600)
corrplot(M_C[rekkefolge, rekkefolge], type="upper",  order="original",
         col= colorRampPalette(c("blue","yellow", "red"))(10), tl.col = "black") #brewer.pal(n=8, name="RdYlBu"))
dev.off()
write.csv2(M_C[rekkefolge, rekkefolge], fs::path(utSti2, "corrplotC.csv"))




dStatus <- matrix(NA, ncol = 5, nrow = ncol(d_alle_scaled_scaled[,fra:til] ))

for(i in fra:til){
  dStatus[i,] <- tapply(d_alle_scaled_scaled[,i], d_alle_scaled_scaled$statusTid, median)
}

colnames(dStatus) <- names(tapply(d_alle_scaled_scaled[,1], d_alle_scaled_scaled$statusTid, median))
rownames(dStatus) <- colnames(d_alle_scaled_scaled)[fra:til]

h <- Heatmap(dStatus, heatmap_legend_param = list(title = ""), row_names_max_width = unit(10, "cm")) #,show_heatmap_legend = FALSE
h

tiff(fs::path(utSti, "heatmap_status.tiff"), width = 700, height = 600)
h
dev.off()

tiff(fs::path(utSti2, "heatmap_status.tiff"), width = 700, height = 600)
h
dev.off()







dStatusT1 <- dStatus[, c("Severe T1", "Moderate T1", "Control")]

h <- Heatmap(dStatusT1, heatmap_legend_param = list(title = ""), row_names_max_width = unit(10, "cm"))
h

tiff(fs::path(utSti, "heatmap_status_T1.tiff"), width = 700, height = 600)
h
dev.off()

tiff(fs::path(utSti2, "heatmap_status_T1.tiff"), width = 700, height = 600)
h
dev.off()





dStatusT2 <- dStatus[, c("Severe T2", "Moderate T2", "Control")]

h <- Heatmap(dStatusT2, heatmap_legend_param = list(title = ""), row_names_max_width = unit(10, "cm"))
h

tiff(fs::path(utSti, "heatmap_status_T2.tiff"), width = 700, height = 600)
h
dev.off()

tiff(fs::path(utSti2, "heatmap_status_T2.tiff"), width = 700, height = 600)
h
dev.off()



ggplot(data = d_alle_scaled, aes(x = `CXCL13, plasma`, y = `Cl6 (plasma cells), Cytof Unstim`, col = statusTid)) +
  
  geom_point(size = 4) + scale_fill_manual(values = colorlegend)



ggplot(data = d_alle_raa, aes(x = 100 *(exp(`CXCL13, plasma`) -0.000001), y = exp(`Cl6 (plasma cells), Cytof Unstim`) -0.000001, col = statusTid)) +
  
  geom_point(size = 4) + scale_fill_manual(values = colorlegend)
