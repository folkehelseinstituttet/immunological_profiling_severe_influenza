scriptPath_UnsupAnalysis <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "Cytof unsupperviced")

library(ggplot2)

stiDes2022 <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "Cytof unsupperviced", "Endelig_des2022",  "Panel 1")
sti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "Cytof unsupperviced",   "Panel 1")


gates <- read.csv2(fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_panel1_mars2022", "posNeg", "Data", "gater.csv"))
navn <- read.csv2( fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_panel1_mars2022", "posNeg", "Data", "navn.csv" ))
navn$marker <- paste0("X", navn$marker)


mapper <- "All"




seeds <- "seed 1345"
temp_unike <- read.csv2(fs::path(sti,mapper[1], seeds[1], paste0("unike_klustre", gsub("seed ", "", seeds[1]), ".csv")))
info <- temp_unike[,c("kort_filnavn", "status", "age", "sex", "statusTid", "tid")]

klustre <- data.frame(id = info[,1])

qTemp <- read.csv2(fs::path(sti,mapper[1], seeds[1], paste0("q5_per_kluster_k_", 10, "_seed", gsub("seed ", "", seeds[1]), mapper[1], ".csv")))
qTom2 <- data.frame(matrix(NA, ncol = ncol(qTemp), nrow = 2))
colnames(qTom2) <- colnames(qTemp)
rownames(qTom2) <- c("slett 1", "slett 2")
qTom2 <- qTom2[,-1]
qTom <- t(qTom2)
qTom <- data.frame(qTom)
qTom$marker <- rownames(qTom)

q5_trans <- qTom
q10_trans <- qTom
q25_trans <- qTom
q50_trans <- qTom
q75_trans <- qTom
q90_trans <- qTom
q95_trans <- qTom


strsplit4 <- function(x,y){
  return(strsplit(x,y)[[1]][4])
}

strsplit6 <- function(x,y){
  return(strsplit(x,y)[[1]][6])
}



signif <- readxl::read_excel(fs::path(scriptPath_UnsupAnalysis, "Endelig_des2022", "Figurer", "151222_Panel1_ST1vsMT1_FINAL.xlsx"))
tall <- as.numeric(signif$`Kluster nr`)
signif <- signif[order(tall),]
Cluster <- signif$`Anjas navnPercent`
signif$Navn_markerplot <- signif$Navn_tSNE
groups <- signif$Navn_markerplot


i <- 1
j <- 1
    temp_unike <- read.csv2(fs::path(stiDes2022,mapper[i], seeds[j], paste0("unike_klustre", gsub("seed ", "", seeds[j]), ".csv")))
    #if(file.info(fs::path(stiDes2022,mapper[i], seeds[j], paste0("sign_klustre_", gsub("seed ", "", seeds[j]), ".csv")))$size > 5){
      temp_sign <- read.csv2(fs::path(stiDes2022,mapper[i], seeds[j], paste0("sign_klustre_", gsub("seed ", "", seeds[j]), ".csv")))
      colnames(temp_sign) <- c("X", "sign")
      temp <-  temp_unike[, as.character(temp_sign$sign)]
      if(!is.null(dim(temp))){
        colnames(temp) <- paste0(mapper[i], "_", colnames(temp))
        antall <- matrix(rep(temp_unike$antall_cells_totalt, ncol(temp_sign)), byrow = FALSE, ncol = ncol(temp_sign))
        klustre <- cbind(klustre, temp/antall)
      } else {
        antall <- temp_unike$antall_cells_totalt
        klustre <- cbind(klustre, temp/antall)
        colnames(klustre)[ncol(klustre)] <- paste0(mapper[i], "_", as.character(temp_sign$sign))
      }


      
      unike_k <-  unique(unlist(lapply(as.character(temp_sign$sign),strsplit4, y = "_")))
      for(k in unike_k){
        #q5
        temp_q5 <- read.csv2(fs::path(sti,mapper[i], seeds[j], paste0("q5_per_kluster_k_", k, "_seed", gsub("seed ", "", seeds[j]), mapper[i], ".csv")))
        temp_q5 <- temp_q5[,-1]
        unike_temp_sign_k <- as.character(temp_sign$sign)[grepl(paste0("k_", k), as.character(temp_sign$sign))]
        unike_k_k <- as.integer(unique(unlist(lapply(as.character(unike_temp_sign_k),strsplit6, y = "_"))))
        q5_temp <- t(temp_q5)
        q5_temp <- data.frame(q5_temp)
        colnames(q5_temp) <- paste0(mapper[i], "_seed_", gsub("seed ", "", seeds[j]), "_k_", k, "_kluster_", gsub("X", "", colnames(q5_temp)))
        q5_temp$marker <- rownames(q5_temp)
        n_marker <- which(colnames(q5_temp) == "marker")
        q5_trans <- merge(q5_trans, q5_temp[, c(unike_k_k, n_marker)], by = "marker",all = T)
        #q10
        temp_q10 <- read.csv2(fs::path(sti,mapper[i], seeds[j], paste0("q10_per_kluster_k_", k, "_seed", gsub("seed ", "", seeds[j]), mapper[i], ".csv")))
        temp_q10 <- temp_q10[,-1]
        unike_temp_sign_k <- as.character(temp_sign$sign)[grepl(paste0("k_", k), as.character(temp_sign$sign))]
        unike_k_k <- as.integer(unique(unlist(lapply(as.character(unike_temp_sign_k),strsplit6, y = "_"))))
        q10_temp <- t(temp_q10)
        q10_temp <- data.frame(q10_temp)
        colnames(q10_temp) <- paste0(mapper[i], "_seed_", gsub("seed ", "", seeds[j]), "_k_", k, "_kluster_", gsub("X", "", colnames(q10_temp)))
        q10_temp$marker <- rownames(q10_temp)
        n_marker <- which(colnames(q10_temp) == "marker")
        q10_trans <- merge(q10_trans, q10_temp[, c(unike_k_k, n_marker)], by = "marker",all = T)
        #q25
        temp_q25 <- read.csv2(fs::path(sti,mapper[i], seeds[j], paste0("q25_per_kluster_k_", k, "_seed", gsub("seed ", "", seeds[j]), mapper[i], ".csv")))
        temp_q25 <- temp_q25[,-1]
        unike_temp_sign_k <- as.character(temp_sign$sign)[grepl(paste0("k_", k), as.character(temp_sign$sign))]
        unike_k_k <- as.integer(unique(unlist(lapply(as.character(unike_temp_sign_k),strsplit6, y = "_"))))
        q25_temp <- t(temp_q25)
        q25_temp <- data.frame(q25_temp)
        colnames(q25_temp) <- paste0(mapper[i], "_seed_", gsub("seed ", "", seeds[j]), "_k_", k, "_kluster_", gsub("X", "", colnames(q25_temp)))
        q25_temp$marker <- rownames(q25_temp)
        n_marker <- which(colnames(q25_temp) == "marker")
        q25_trans <- merge(q25_trans, q25_temp[, c(unike_k_k, n_marker)], by = "marker",all = T)
        #q50
        temp_q50 <- read.csv2(fs::path(sti,mapper[i], seeds[j], paste0("medians_per_kluster_k_", k, "_seed", gsub("seed ", "", seeds[j]), mapper[i], ".csv")))
        temp_q50 <- temp_q50[,-1]
        unike_temp_sign_k <- as.character(temp_sign$sign)[grepl(paste0("k_", k), as.character(temp_sign$sign))]
        unike_k_k <- as.integer(unique(unlist(lapply(as.character(unike_temp_sign_k),strsplit6, y = "_"))))
        q50_temp <- t(temp_q50)
        q50_temp <- data.frame(q50_temp)
        colnames(q50_temp) <- paste0(mapper[i], "_seed_", gsub("seed ", "", seeds[j]), "_k_", k, "_kluster_", gsub("X", "", colnames(q50_temp)))
        q50_temp$marker <- rownames(q50_temp)
        n_marker <- which(colnames(q50_temp) == "marker")
        q50_trans <- merge(q50_trans, q50_temp[, c(unike_k_k, n_marker)], by = "marker",all = T)
        #q5
        temp_q75 <- read.csv2(fs::path(sti,mapper[i], seeds[j], paste0("q75_per_kluster_k_", k, "_seed", gsub("seed ", "", seeds[j]), mapper[i], ".csv")))
        temp_q75 <- temp_q75[,-1]
        unike_temp_sign_k <- as.character(temp_sign$sign)[grepl(paste0("k_", k), as.character(temp_sign$sign))]
        unike_k_k <- as.integer(unique(unlist(lapply(as.character(unike_temp_sign_k),strsplit6, y = "_"))))
        q75_temp <- t(temp_q75)
        q75_temp <- data.frame(q75_temp)
        colnames(q75_temp) <- paste0(mapper[i], "_seed_", gsub("seed ", "", seeds[j]), "_k_", k, "_kluster_", gsub("X", "", colnames(q75_temp)))
        q75_temp$marker <- rownames(q75_temp)
        n_marker <- which(colnames(q75_temp) == "marker")
        q75_trans <- merge(q75_trans, q75_temp[, c(unike_k_k, n_marker)], by = "marker",all = T)
        #q5
        temp_q90 <- read.csv2(fs::path(sti,mapper[i], seeds[j], paste0("q90_per_kluster_k_", k, "_seed", gsub("seed ", "", seeds[j]), mapper[i], ".csv")))
        temp_q90 <- temp_q90[,-1]
        unike_temp_sign_k <- as.character(temp_sign$sign)[grepl(paste0("k_", k), as.character(temp_sign$sign))]
        unike_k_k <- as.integer(unique(unlist(lapply(as.character(unike_temp_sign_k),strsplit6, y = "_"))))
        q90_temp <- t(temp_q90)
        q90_temp <- data.frame(q90_temp)
        colnames(q90_temp) <- paste0(mapper[i], "_seed_", gsub("seed ", "", seeds[j]), "_k_", k, "_kluster_", gsub("X", "", colnames(q90_temp)))
        q90_temp$marker <- rownames(q90_temp)
        n_marker <- which(colnames(q90_temp) == "marker")
        q90_trans <- merge(q90_trans, q90_temp[, c(unike_k_k, n_marker)], by = "marker",all = T)
        #q5
        temp_q95 <- read.csv2(fs::path(sti,mapper[i], seeds[j], paste0("q95_per_kluster_k_", k, "_seed", gsub("seed ", "", seeds[j]), mapper[i], ".csv")))
        temp_q95 <- temp_q95[,-1]
        unike_temp_sign_k <- as.character(temp_sign$sign)[grepl(paste0("k_", k), as.character(temp_sign$sign))]
        unike_k_k <- as.integer(unique(unlist(lapply(as.character(unike_temp_sign_k),strsplit6, y = "_"))))
        q95_temp <- t(temp_q95)
        q95_temp <- data.frame(q95_temp)
        colnames(q95_temp) <- paste0(mapper[i], "_seed_", gsub("seed ", "", seeds[j]), "_k_", k, "_kluster_", gsub("X", "", colnames(q95_temp)))
        q95_temp$marker <- rownames(q95_temp)
        n_marker <- which(colnames(q95_temp) == "marker")
        q95_trans <- merge(q95_trans, q95_temp[, c(unike_k_k, n_marker)], by = "marker",all = T)
      }
    
 
rownames(q5_trans) <- q5_trans$marker
q5_trans <- q5_trans[,-c(1,2,3)]
q5_all <- t(q5_trans)

rownames(q10_trans) <- q10_trans$marker
q10_trans <- q10_trans[,-c(1,2,3)]
q10_all <- t(q10_trans)

rownames(q25_trans) <- q25_trans$marker
q25_trans <- q25_trans[,-c(1,2,3)]
q25_all <- t(q25_trans)

rownames(q50_trans) <- q50_trans$marker
q50_trans <- q50_trans[,-c(1,2,3)]
q50_all <- t(q50_trans)

rownames(q75_trans) <- q75_trans$marker
q75_trans <- q75_trans[,-c(1,2,3)]
q75_all <- t(q75_trans)

rownames(q90_trans) <- q90_trans$marker
q90_trans <- q90_trans[,-c(1,2,3)]
q90_all <- t(q90_trans)

rownames(q95_trans) <- q95_trans$marker
q95_trans <- q95_trans[,-c(1,2,3)]
q95_all <- t(q95_trans)

groups <- factor(groups, levels = groups)
q5 <- data.frame(Group = groups,  q5_all[paste0("All_", Cluster),])
q10 <- data.frame(Group = groups,  q10_all[paste0("All_", Cluster),])
q25 <- data.frame(Group = groups,  q25_all[paste0("All_", Cluster),])
q50 <- data.frame(Group = groups,  q50_all[paste0("All_", Cluster),])
q75 <- data.frame(Group = groups,  q75_all[paste0("All_", Cluster),])
q90 <- data.frame(Group = groups,  q90_all[paste0("All_", Cluster),])
q95 <- data.frame(Group = groups,  q95_all[paste0("All_", Cluster),])


q0_temp <- read.csv2(fs::path(sti,mapper[1], seeds[1],  "quantilesAllOthers.csv"))


q0 <- data.frame(Group = "Cl0 (Other cells, ns)", q0_temp[,colnames(q5)[2:ncol(q5)]])
  
q5 <- data.frame(Group = c("Cl0 (Other cells, ns)", as.character(q5$Group)), rbind(q0[1,-1], q5[,-1]))
q10 <- data.frame(Group = c("Cl0 (Other cells, ns)", as.character(q10$Group)), rbind(q0[2,-1], q10[,-1]))
q25 <- data.frame(Group = c("Cl0 (Other cells, ns)", as.character(q25$Group)), rbind(q0[3,-1], q25[,-1]))
q50 <- data.frame(Group = c("Cl0 (Other cells, ns)", as.character(q50$Group)), rbind(q0[4,-1], q50[,-1]))
q75 <- data.frame(Group = c("Cl0 (Other cells, ns)", as.character(q75$Group)), rbind(q0[5,-1], q75[,-1]))
q90 <- data.frame(Group = c("Cl0 (Other cells, ns)", as.character(q90$Group)), rbind(q0[6,-1], q90[,-1]))
q95 <- data.frame(Group = c("Cl0 (Other cells, ns)", as.character(q95$Group)), rbind(q0[7,-1], q95[,-1]))

 

  q5long <- data.table::melt(data.table::setDT(q5),id.vars = "Group", variable.names = "marker")
  colnames(q5long)[3] <- "q5"
  q10long <- data.table::melt(data.table::setDT(q10),id.vars = "Group", variable.names = "marker")
  colnames(q10long)[3] <- "q10"
  q25long <- data.table::melt(data.table::setDT(q25),id.vars = "Group", variable.names = "marker")
  colnames(q25long)[3] <- "q25"
  q50long <- data.table::melt(data.table::setDT(q50),id.vars = "Group", variable.names = "marker")
  colnames(q50long)[3] <- "q50"
  q75long <- data.table::melt(data.table::setDT(q75),id.vars = "Group", variable.names = "marker")
  colnames(q75long)[3] <- "q75"
  q90long <- data.table::melt(data.table::setDT(q90),id.vars = "Group", variable.names = "marker")
  colnames(q90long)[3] <- "q90"
  q95long <- data.table::melt(data.table::setDT(q95),id.vars = "Group", variable.names = "marker")
  colnames(q95long)[3] <- "q95"
  
  d <-merge(q5long, q10long)
  d <-merge(d, q25long)
  d <-merge(d, q50long)
  d <-merge(d, q75long)
  d <- merge(d, q90long)
  d <- merge(d, q95long)
  
  
  
  d <- merge(d, navn, by.x= "variable", by.y = "marker")
  if(!is.null(gates)){
    d <- merge(d, gates, by = "marker_short_name")
  }
  
   d$marker_short_name <- factor(d$marker_short_name)
  
  
 
  color_group <- c("#CCCCCC", signif$colour)
  names(color_group) <- c("Cl0 (Other cells, ns)",  signif$Navn_markerplot)

  
  
  d$col <- color_group[d$Group]
   
#  d$Group <- factor(d$Group, levels = rev(names(color_group)))
  
                     
d$marker_short_name <- factor(d$marker_short_name, levels = c("CCR4", "CCR6", "CCR7", "CD3", "CD4", "CD5", "CD8",   
                                                              "CD11b", "CD11c", "CD14", "CD15", "CD16", "CD19", "CD25", 
                                                              "CD27", "CD28", "CD38", "CD45", "CD45RA", "CD56", "CD57",  "CD85j", "CD95", 
                                                              "CD123", "CD127", "CD134", "CD141", "CD160", "CD161", "CD169",     
                                                              "CRTH2", "CXCR3", "CXCR5", "HLADR", "ICOS", "IgD", "IgG", "KLRG1", "NKG2A",
                                                              "PD-1", "TCRVa7.2", "TCRgd", "TIGIT"))                  
               
  
  unique_col <- levels(d$col)
  
  d$Group <- factor(d$Group, rev(c("Cl0 (Other cells, ns)", as.character(groups))))

  g <- ggplot(d, aes(x = q50, y = Group, xmin = q10, xmax = q90, col = Group)) +
    geom_point(size = 2) +
    geom_errorbar() +
    geom_errorbar(data = d, aes(x = q50, y = Group, xmin = q5, xmax = q95, col = Group))+
    geom_errorbar(data = d, aes(x = q50, y = Group, xmin = q25, xmax = q75, col = Group))+
    facet_wrap(vars(marker_short_name), ncol = 11, scale = "free_x") + theme_classic(base_size = 18) +
    scale_color_manual(values = color_group) + xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + theme(axis.title = element_text(face="bold")) +  theme(legend.position = "none")
  
  
  if(!is.null(gates)){
    g <- g + geom_vline(aes(xintercept  = low), linetype = "dashed", linewidth = 1.1) +
      geom_vline(aes(xintercept  = high), linetype = "dashed", linewidth = 1.1)
  }
  
  
  scriptPath_UnsupAnalysisUt <- fs::path(scriptPath_UnsupAnalysis, "Endelig_des2022",  "Figurer")
  tiff( fs::path(scriptPath_UnsupAnalysisUt, "P1_marker_plot.tiff"), width = 1350, height = 900)
  print(g)
  dev.off()
  
  
  
  
  
  
  