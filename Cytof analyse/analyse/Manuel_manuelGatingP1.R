posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel1_mars2022", "posNeg", "Data")
fig_path <- fs::path("F:" ,"Forskningsprosjekter" ,"PDB 2794 - Immune responses aga_" ,"Forskningsfiler" ,"JOBO" ,"CyTOF" ,"Analyse i R OUS" ,"CleanUpGatingMarch2022" ,"Analyse" ,"Endelig" , "Cytof unsupperviced", "Endelig_des2022", "Figurer", "Manuel_ManuelGating")
functionPath <- fs::path("H:", "git", "cytof")
source(fs::path(functionPath, "functions.R"))

posNegFilnavn <- as.character(readRDS(fs::path(posNeg_path,  "posNegFilnavn.rds")))
posNeg <- readRDS(fs::path(posNeg_path,  "posNeg.rds"))
kort_fil_navn <- read.csv2(fs::path(posNeg_path, "kortFilNavn.csv"))






for(i in 1:length(posNeg)){

  
  
  posNeg[[i]]$cellname_Bagwell <- "unknown"
  # posNeg[[i]]$cellname_Bagwell[posNeg[[i]]$CD45 == 0] <- "CD45 neg" 
  # posNeg[[i]]$cellname_Bagwell[posNeg[[i]]$CD45 == 1 & posNeg[[i]]$CD19 == 0 & posNeg[[i]]$CD3 == 1 ] <- "CD3" 
#  posNeg[[i]]$cellname_Bagwell[posNeg[[i]]$CD45 == 1 & posNeg[[i]]$CD19 == 0 & posNeg[[i]]$CD3 == 1 & posNeg[[i]]$CD4 == 1 & posNeg[[i]]$CD8 == 0 & posNeg[[i]]$CD25 %in% c(1,2) & posNeg[[i]]$CD127 %in% c(0)] <- "Tregs" 
  posNeg[[i]]$cellname_Bagwell[posNeg[[i]]$CD45 == 1 & posNeg[[i]]$CD19 == 0 & posNeg[[i]]$CD3 == 1 & posNeg[[i]]$CD56 %in% c(2) & posNeg[[i]]$TCRgd == 0] <- "NKT"
  posNeg[[i]]$cellname_Bagwell[posNeg[[i]]$CD45 == 1 & posNeg[[i]]$CD19 == 0 & posNeg[[i]]$CD3 == 1 & posNeg[[i]]$CD56 %in% c(1,2) & posNeg[[i]]$TCRgd == 0 & posNeg[[i]]$CD57 %in% c(2)] <- "NKT"
  posNeg[[i]]$cellname_Bagwell[posNeg[[i]]$CD45 == 1 & posNeg[[i]]$CD19 == 0 & posNeg[[i]]$CD3 == 1 & posNeg[[i]]$CD4 == 1 & posNeg[[i]]$CD8 == 0] <- "CD4" 
  posNeg[[i]]$cellname_Bagwell[posNeg[[i]]$CD45 == 1 & posNeg[[i]]$CD19 == 0 & posNeg[[i]]$CD3 == 1 & posNeg[[i]]$CD4 == 0 & posNeg[[i]]$CD8 %in% c(1,2)] <- "CD8"
  posNeg[[i]]$cellname_Bagwell[posNeg[[i]]$CD45 == 1 & posNeg[[i]]$CD19 == 0 & posNeg[[i]]$CD3 == 1 & posNeg[[i]]$TCRVa7.2 %in% c(0) & posNeg[[i]]$TCRgd %in% c(1)] <- "gdT"
  posNeg[[i]]$cellname_Bagwell[posNeg[[i]]$CD45 == 1 & posNeg[[i]]$CD19 == 0 & posNeg[[i]]$CD3 == 0 & (posNeg[[i]]$CD56 %in% c(1,2) | posNeg[[i]]$CD57 %in% c(2)) & posNeg[[i]]$CD123 == 0] <- "NK"
  posNeg[[i]]$cellname_Bagwell[posNeg[[i]]$CD45 == 1 & posNeg[[i]]$CD19 == 1 & posNeg[[i]]$CD3 == 0 & posNeg[[i]]$CD56 == 0 & posNeg[[i]]$HLADR %in% c(1,2)] <- "B"
  posNeg[[i]]$cellname_Bagwell[posNeg[[i]]$CD45 == 1 & posNeg[[i]]$CD19 == 0 & posNeg[[i]]$CD3 == 0 & posNeg[[i]]$CD57 == 0 & posNeg[[i]]$HLADR %in% c(2) & posNeg[[i]]$CD16 == 0 & posNeg[[i]]$CD169 == 0] <- "DC"
  posNeg[[i]]$cellname_Bagwell[posNeg[[i]]$CD45 == 1 & posNeg[[i]]$CD19 == 0 & posNeg[[i]]$CD3 == 0 & posNeg[[i]]$CD56 %in% c(0) & posNeg[[i]]$CD11c %in% c(1,2) &  posNeg[[i]]$CD11b == 1 &  posNeg[[i]]$CD57 == 0 & (posNeg[[i]]$CD14 == c(1,2))] <- "Mo" 
  posNeg[[i]]$cellname_Bagwell[posNeg[[i]]$CD45 == 1 & posNeg[[i]]$CD19 == 0 & posNeg[[i]]$CD3 == 0 & posNeg[[i]]$HLADR %in% c(1,2) & posNeg[[i]]$CD11c %in% c(1,2) &  posNeg[[i]]$CD11b == 1 &  posNeg[[i]]$CD56 == 0 & posNeg[[i]]$CD14 == 0 & posNeg[[i]]$CD38 == 0 ] <- "Mo" 
  posNeg[[i]]$cellname_Bagwell[posNeg[[i]]$CD45 == 1 & posNeg[[i]]$CD3 == 0 & posNeg[[i]]$CD19 == 0 & posNeg[[i]]$CD123 == 1 & posNeg[[i]]$CD11c == 0 & posNeg[[i]]$HLADR == c(1,2) & posNeg[[i]]$CD14 == 0 & posNeg[[i]]$CD169 == 0] <- "DC" 
  posNeg[[i]]$cellname_Bagwell[posNeg[[i]]$CD45 == 1 & posNeg[[i]]$CD3 == 0 & posNeg[[i]]$CD19 == 0 & posNeg[[i]]$CD123 == 0 & posNeg[[i]]$CD11c %in% c(1,2) & posNeg[[i]]$HLADR %in% c(1,2) & posNeg[[i]]$CD14 == 0 & posNeg[[i]]$CD16 %in% c(1,0) & posNeg[[i]]$CD38 %in% c(1,0) & posNeg[[i]]$CRTH2 == 0 & posNeg[[i]]$CD169 == 0] <- "DC" 
  posNeg[[i]]$cellname_Bagwell[posNeg[[i]]$CD45 == 1 & posNeg[[i]]$CD19 == 0 & posNeg[[i]]$CD3 == 0 & posNeg[[i]]$CD56 == 0 & posNeg[[i]]$CD169 == 1] <- "Mo" 
  posNeg[[i]]$cellname_Bagwell[posNeg[[i]]$CD45 == 1 & posNeg[[i]]$CD19 == 0 & posNeg[[i]]$CD3 == 1 & posNeg[[i]]$TCRVa7.2 %in% c(1) & posNeg[[i]]$TCRgd %in% c(0)] <- "MAIT"
  
}



                                                                                                     

mulige_celler <- sort(unique(posNeg[[1]]$cellname_Bagwell))

cells <- as.data.frame(matrix(ncol = length(mulige_celler), nrow = length(posNeg) ))
colnames(cells) <- mulige_celler
rownames(cells) <- kort_fil_navn$kort_filnavn

for(i in 1:nrow(cells)){
  z <- table(posNeg[[i]]$cellname_Bagwell)
  cells[i,names(z)] <- z
}

cells[is.na(cells)] <- 0

cells$totalt_antall <- apply(cells, 1, sum)


cells$StatusTid <- kort_fil_navn$status

cells$age <- kort_fil_navn$age

cells$sex <- kort_fil_navn$sex


cells <- cells[cells$StatusTid %in% c("Control", "Severe T1", "Severe T2", "Moderate T1", "Moderate T2"),]

cells$StatusTidSex <- paste(cells$StatusTid, cells$sex)

write.csv2(cells, fs::path(fig_path, "P1_cells_per_manuel_cluster_and_sample.csv"))



