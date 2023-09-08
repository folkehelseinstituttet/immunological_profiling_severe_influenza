
library(ggplot2)
library(cowplot)
library(ggpubr)

# d: data frame with one row per cell, coloum

markerplot <- function(d, markers = "all", colnamesCluster = "Clusters"  , clusters = "all", color_clusters = NA, with_percent = FALSE, gates = NULL, returnDF = FALSE, ncol = 11){
    
  
  
  q5to95 <- function(x){
    quantile(x, c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95))
  }
   
  n <- nrow(d)
     
    
    if(colnamesCluster %in% colnames(d)){
      i_clusters <- which(colnames(d) == colnamesCluster)
    } else {
      stop("colnamesClusters has to be the name of the column with clusternames")
    } 
    
    if(!clusters[1] == "all"){
      d <- d[d[,i_clusters] %in% clusters,]
      d[,i_clusters] <- factor(d[,i_clusters], levels = clusters)
    } else {
      clusters <- sort(unique(d[,i])) 
    }

    
    if(is.na(color_clusters)[1]){
      color_clusters <- rep("grey", length(clusters))
      names(color_clusters) <- clusters
    }
    


    tmp <- data.frame(expand.grid(markers, clusters))
    colnames(tmp) <- c("markers", "clusters")
    tmp$q5 <- NA
    tmp$q10 <- NA
    tmp$q25 <- NA
    tmp$q50 <- NA
    tmp$q75 <- NA
    tmp$q90 <- NA
    tmp$q95 <- NA
    rownames(tmp) <- paste(tmp[,1], tmp[,2], sep = "_")
    
    for(i in 1:length(clusters)){
      di <- d[d[, i_clusters] == clusters[i],]
      for(j in 1:length(markers)){
        mi <- markers[j]
        name <- paste(mi, clusters[i], sep = "_")
        tmp[name,3:9] <- q5to95(di[,mi])  
      }
    }
    
    
    if(with_percent){
      percent <- round(table(d[,i_clusters])/n * 100, 1)
      percent <- data.frame(clusters = names(percent), percent = as.numeric(percent))
      tmp <- merge(tmp, percent, by = "clusters")
      tmp$clusters <- paste0(tmp$clusters, " (", tmp$percent, "%)")
      percent$name <- paste0(percent$clusters, " (", percent$percent, "%)")
      names(color_clusters) <- percent$name
    }
    
    if(!is.null(gates)[1]){
      tmp <- merge(tmp, gates, by = "markers")
    }
    
    tmp$clusters <- factor(tmp$clusters, levels = names(color_clusters))    
    
    g <- ggplot(tmp, aes(x = q50, y = clusters, xmin = q10, xmax = q90, col = clusters)) +
      geom_point(size = 2) +
      geom_errorbar() +
      geom_errorbar(data = tmp, aes(x = q50, y = clusters, xmin = q5, xmax = q95, col = clusters))+
      geom_errorbar(data = tmp, aes(x = q50, y = clusters, xmin = q25, xmax = q75, col = clusters))+
      facet_wrap(vars(markers), ncol = ncol, scale = "free_x") + theme_classic(base_size = 18) +
      scale_color_manual(values = color_clusters) + xlab("") + ylab("") +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + theme(axis.title = element_text(face="bold")) +  theme(legend.position = "none")
    
    
    if(!is.null(gates)){
      g <- g + geom_vline(aes(xintercept  = low), linetype = "dashed", linewidth = 1.1) +
        geom_vline(aes(xintercept  = high), linetype = "dashed", linewidth = 1.1)
    }
    
    if(returnDF){
      return(list(g = g, tmp = tmp))
    } else {
      return(g)
    }
}



### user example


dataSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "Cytof unsupperviced", "Panel 1", "All", "seed 1345")
d2 <- readRDS(fs::path(dataSti, "data_med_filnavn1345.RDS"))
#d2 should look like this:
# > d2[1:2,]
# CD45       CD3       CD4      CD8      CD19 TCRgd  TCRVa7.2      CD25     CD27     CD28     CD127     CD38
# 1 2.643928 0.4873854 0.7138038 0.000000 0.0000000     0 0.0000000 0.4213454 0.000000 0.000000 0.2400642 2.138743
# 2 4.290144 3.5603067 2.5441205 0.550012 0.2687022     0 0.5780855 0.8425878 4.386372 2.360449 1.0967351 1.177777
# CCR7   CD45RA     HLADR       IgD       IgG     CXCR3      CCR4      CCR6     CXCR5      ICOS     KLRG1 CD134
# 1 0.000000 1.646066 2.3992607 0.1496770 4.6169348 1.5275158 0.8226483 0.1837734 0.0000000 0.4830149 0.4444857     0
# 2 4.522755 3.575275 0.4251954 0.5085317 0.2055259 0.1017059 0.6107758 0.0000000 0.4050297 2.4103504 0.0000000     0
# TIGIT     CD160      CD161      CD95     CRTH2     CD85j NKG2A    CD141 PD-1       CD14      CD16       CD56
# 1 0.08092619 0.4186461 0.15365474 3.1055616 0.8770397 3.0586106     0 1.263599    0 0.92448137 0.6989533 0.15952566
# 2 0.26905544 0.3056073 0.05781532 0.8693487 0.2068264 0.3444777     0 0.000000    0 0.05107892 0.1654434 0.02069772
# CD57     CD11b    CD11c       CD5 CD15 CD169     CD123 kluster 10 kluster 20 kluster 30 kluster 40 kluster 50
# 1 0.0000000 2.2297068 3.093126 0.7875798    0     0 0.4270527          1          4          4          6          6
# 2 0.6645122 0.1799299 0.000000 3.3464267    0     0 0.0000000          3         13         19         24         30
# kluster 60                 dataset
# 1          6 Ko_Fe_A_66_542K1_Panel1
# 2         32 Ko_Fe_A_66_542K1_Panel1
#NB you only need one column with clustering, here I have 6. 

gates <- read.csv2(fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_panel1_mars2022", "posNeg", "Data", "gater.csv"))
colnames(gates) <- c("markers", "low", "high")
#> head(gates)
# markers      low     high
# 1               CD3 1.447716       NA
# 2               CD4 2.431271       NA
# 3               CD5 1.676759       NA
# 4               CD8 1.000000 3.327063
# 5              CD19 2.000613       NA
# 6              CD45 1.638044       NA
#


# #this is specific for my data. don't look at it.
# navn <- read.csv2( fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_panel1_mars2022", "posNeg", "Data", "navn.csv" ))
# navn$marker[navn$marker == "174Yb_CD279_PD.1"] <- "174Yb_CD279_PD-1"
# rownames(navn) <- navn$marker
# 
# nyeKolNavn <- rep(NA, 43)
# for(i in 1:43){
#  nyeKolNavn[i] <- navn[colnames(d2)[i], "marker_short_name"]
# }
# colnames(d2)[1:43] <- nyeKolNavn
# #until here...
# 
colnamesCluster <- "kluster 60" #column name of the column which have the clustering in it 
clusters <- c("1", "3","5","10","55") #name of the clusters you want to have in the marker plot. If all you might just write "all" defalut.
color_clusters <- c("1" ="#33CCCC","3" = "#1F78B4", "5" = "#E31A1C", "10" = "#FB9A99", "55" = "#FDBF6F") #default is all grey.


markers <- c("CCR4", "CCR6", "CCR7", "CD3", "CD4", "CD5", "CD8",   
             "CD11b", "CD11c", "CD14", "CD15", "CD16", "CD19", "CD25", 
             "CD27", "CD28", "CD38", "CD45", "CD45RA", "CD56", "CD57",  "CD85j", "CD95", 
             "CD123", "CD127", "CD134", "CD141", "CD160", "CD161", "CD169",     
             "CRTH2", "CXCR3", "CXCR5", "HLADR", "ICOS", "IgD", "IgG", "KLRG1", "NKG2A",
             "PD-1", "TCRVa7.2", "TCRgd", "TIGIT")  #the order you want the plot in . Have to include all the markers you want to plot
# 
# # colnames(d)[1:43] #vector of the markers you want in the marker plot with same name as column names in d

#gates a matrix with the columns: markers, low and high (is only made for two possible gates, can have only low). if included will give vertical lines in the plot
d <- d2

#without vertical gates
markerplot(d = d, markers = markers, colnamesCluster = "kluster 60", clusters = clusters, color_clusters = color_clusters, with_percent = TRUE, gates = NULL)


#with vertical gates
markerplot(d = d, markers = markers, colnamesCluster = "kluster 60", clusters = clusters, color_clusters = color_clusters, with_percent = TRUE, gates = gates)


markerplot_multipleclusteringColumns <- function(d, markers = "all", colnamesCluster = "Clusters"  , clusters = "all",  with_percent = FALSE, gates = NULL, returnDF = FALSE, ncol = 11){


  q5to95 <- function(x){
    quantile(x, c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95))
  }

  n <- nrow(d)
  d_all <- d



 
  if(all(colnamesCluster %in% colnames(d))){
    i_clusters <- which(colnames(d) %in% colnamesCluster)
  } else {
    stop("colnamesClusters has to be the name of columns with clusternames")
  }

d <- NULL



  if(is.list(clusters)){
    if(all(names(clusters) %in% colnamesCluster)){
      clusters <- clusters[colnamesCluster[colnamesCluster %in% names(clusters)]]
      colnamesCluster <- colnamesCluster[colnamesCluster %in% names(clusters)]
      i_clusters <- which(colnames(d_all) %in% colnamesCluster)
    } else {
      stop("the name of each element in the list clusters has to correspond to names in colnamesCluster")
    }
    for(i_clunr in 1:length(i_clusters)){
      i_clu <- i_clusters[i_clunr]
      clus <- clusters[[i_clunr]]
      temp <-  d_all[d_all[,i_clu] %in% clus, ]
      temp$clusters <- paste(colnames(temp)[i_clu], temp[,i_clu], sep = "_")
      d <- rbind(d, temp)
    }
  } else {
    if(clusters[1] == "all"){
      i_clusters <- which(colnames(d) %in% colnamesCluster)
      for(i_clu in i_clusters){
        temp <- d_all
        temp$clusters <- paste(colnames(temp)[i_clu], temp[,i_clu], sep = "_")
        d <- rbind(d, temp)  
      }
    }
  }

clusters <- unique(d$clusters)



  # if(is.na(color_clusters)[1]){
  #   color_clusters <- rep("grey", length(clusters))
  #   names( color_clusters) <- clusters
  # }
  # 


  tmp <- data.frame(expand.grid(markers, clusters))
  colnames(tmp) <- c("markers", "clusters")
  tmp$q5 <- NA
  tmp$q10 <- NA
  tmp$q25 <- NA
  tmp$q50 <- NA
  tmp$q75 <- NA
  tmp$q90 <- NA
  tmp$q95 <- NA
  rownames(tmp) <- paste(tmp[,1], tmp[,2], sep = "_")
  
  
 

  for(i in 1:length(clusters)){
    di <- d[d[, "clusters"] == clusters[i],]
    for(j in 1:length(markers)){
      mi <- markers[j]
      name <- paste(mi, clusters[i], sep = "_")
      tmp[name,3:9] <- q5to95(di[,mi])
    }
  }

  if(with_percent){
    percent <- round(table(d[,"clusters"])/n * 100, 2)
    percent <- data.frame(clusters = names(percent), percent = as.numeric(percent))
    tmp <- merge(tmp, percent, by = "clusters")
    tmp$clusters <- paste0(tmp$clusters, " (", tmp$percent, "%)")
    percent$name <- paste0(percent$clusters, " (", percent$percent, "%)")
  #  names(color_clusters) <- percent$name
  }

  
  
  if(!is.null(gates)[1]){
    tmp <- merge(tmp, gates, by = "markers")
  }



  g <- ggplot(tmp, aes(x = q50, y = clusters, xmin = q10, xmax = q90, col = clusters)) +
    geom_point(size = 2) +
    geom_errorbar() +
    geom_errorbar(data = tmp, aes(x = q50, y = clusters, xmin = q5, xmax = q95, col = clusters))+
    geom_errorbar(data = tmp, aes(x = q50, y = clusters, xmin = q25, xmax = q75, col = clusters))+
    facet_wrap(vars(markers), ncol = ncol, scale = "free_x") + theme_classic(base_size = 18) +
    # scale_color_manual(values = color_clusters) +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + theme(axis.title = element_text(face="bold")) +  theme(legend.position = "none")


  if(!is.null(gates)){
    g <- g + geom_vline(aes(xintercept  = low), linetype = "dashed", linewidth = 1.1) +
      geom_vline(aes(xintercept  = high), linetype = "dashed", linewidth = 1.1)
  }
  
  if(returnDF){
    return(list(g = g, qTable = tmp))
  } else {
    return(g)
  }
}


#example


dataSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "Cytof unsupperviced", "Endelig_des2022", "Figurer")
#only Anja specific
# cluNavn <- read.csv2(fs::path(dataSti, "klusternavnP1_JB.csv"))
# sort(as.numeric(gsub("20_", "", cluNavn$Anja_navn[grep("20_", cluNavn$Anja_navn)])))
# sort(as.numeric(gsub("30_", "", cluNavn$Anja_navn[grep("30_", cluNavn$Anja_navn)])))
# sort(as.numeric(gsub("40_", "", cluNavn$Anja_navn[grep("40_", cluNavn$Anja_navn)])))
# sort(as.numeric(gsub("50_", "", cluNavn$Anja_navn[grep("50_", cluNavn$Anja_navn)])))
# sort(as.numeric(gsub("60_", "", cluNavn$Anja_navn[grep("60_", cluNavn$Anja_navn)])))
#until here


colnamesCluster <- c("kluster 10", "kluster 20", "kluster 30", "kluster 40", "kluster 50", "kluster 60")
clusters <- list( "kluster 10" = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10) ,  
                  "kluster 20" = c(1, 3, 4, 6, 7, 8, 9, 13, 14, 15, 16, 17, 19), 
                  "kluster 30" = c(1, 3, 4, 5, 6, 9, 13, 14, 15, 16, 17, 18, 19, 20, 21, 23, 25, 26),
                  "kluster 40" = c(1,  3,  4,  5,  7,  8, 10, 15, 19, 25, 27, 28, 31, 32, 33, 37, 39, 40), 
                  "kluster 50" = c(14, 15, 17, 18, 21, 22, 24, 25, 26, 29, 34, 36, 39, 42, 43, 44, 46), 
                  "kluster 60" = c(11, 22, 29, 32, 36, 38, 40, 41, 44, 48, 49, 51, 54, 55, 56, 58))
d <- d2

markers <- c("CCR4", "CCR6", "CCR7", "CD3", "CD4", "CD5", "CD8",   
             "CD11b", "CD11c", "CD14", "CD15", "CD16", "CD19", "CD25", 
             "CD27", "CD28", "CD38", "CD45", "CD45RA", "CD56", "CD57",  "CD85j", "CD95", 
             "CD123", "CD127", "CD134", "CD141", "CD160", "CD161", "CD169",     
             "CRTH2", "CXCR3", "CXCR5", "HLADR", "ICOS", "IgD", "IgG", "KLRG1", "NKG2A",
             "PD-1", "TCRVa7.2", "TCRgd", "TIGIT")  #the order you want the plot in . 
#or just in the same order as in the matrix d:
#colnames(d)[1:43] #vector of the markers you want in the marker plot with same name as column names in d

tmp <- markerplot_multipleclusteringColumns(d = d, markers = markers, colnamesCluster = colnamesCluster, clusters = clusters,  with_percent = TRUE, gates = gates, returnDF = T)

qTable <- tmp$qTable

tmp$g #this might be best as a tiff where you can make the hight large. 


#this is to specific color in markerplots with clusters from more than one clustering ...


  
markerplot_fromQuantileTable <- function(qTable, colnamesCluster = "Clusters", color_group = color_group, color_group_mindre = color_group_mindre, ncol = 11){
  
  
  
  tmp <- qTable    
  i_col <- which(colnames(tmp) == colnamesCluster)
  tmp$clustersToPlot <- tmp[,i_col]
  tmp$clustersToPlot <- factor(tmp$clustersToPlot, levels = unique(tmp$clustersToPlot))
  
 
  
  g <- ggplot(tmp, aes(x = q50, y = clustersToPlot, xmin = q10, xmax = q90, col = clustersToPlot)) +
    geom_point(size = 2) +
    geom_errorbar() +
    geom_errorbar(data = tmp, aes(x = q50, y = clustersToPlot, xmin = q5, xmax = q95, col = clustersToPlot))+
    geom_errorbar(data = tmp, aes(x = q50, y = clustersToPlot, xmin = q25, xmax = q75, col = clustersToPlot))+
    facet_wrap(vars(markers), ncol = ncol, scale = "free_x") + theme_classic(base_size = 18) +
    scale_color_manual(values = color_group) + xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + theme(axis.title = element_text(face="bold")) +  theme(legend.position = "none")
  
  
  if(!is.null(tmp$low)){
    g <- g + geom_vline(aes(xintercept  = low), linetype = "dashed", linewidth = 1.1) +
      geom_vline(aes(xintercept  = high), linetype = "dashed", linewidth = 1.1)
  }
  
  dtemp <- data.frame(col = color_group_mindre, cellname = names(color_group_mindre), x = 1:length(color_group_mindre),  y = 1:length(color_group_mindre))
  
  g2 <- ggplot(dtemp, aes(x = x, y = y, col = cellname)) +
    geom_point(shape = 15, size = 5) +
    # #  xlim(Min1,Max1) +
    # #  ylim(Min2,Max2) +
    # ggtitle(paste0("")) +
    theme(legend.text=element_text(size=15))  +
   
    scale_colour_manual(values = as.character(color_group_mindre), breaks = names(color_group_mindre),
                        limits = names(color_group_mindre))  +
   
    theme(legend.title=element_blank())
  
  leg <- get_legend(g2)
  gleg <- as_ggplot(leg)
  gleg
  
  
  return(list(g = g, gleg = gleg))

}


#example, have to use markerplot_multipleclusteringColumns with "returnDF = T"
#then change name of columns if wanted and give color per cluster, remember same name

qTable$navn <-gsub("kluster ", "", qTable$clusters )
qTable$navn <- substr(qTable$navn, 1,5)
qTable$navn <-gsub(" ", "", qTable$navn )
qTable <- merge(qTable, cluNavn[, c("Anja_navn", "kort_navn")], by.x = "navn", by.y = "Anja_navn")

color_group_mindre <- c("#33CCCC", "#1F78B4", "#E31A1C",  "#FB9A99", "#FDBF6F", "#9933CC", "#B2DF8A", "#CAB2D6", "#f708db", "#fd028f",  "#FF7F00", "#33A02C", "#999999", "#CCCCCC")
names(color_group_mindre) <- c("CD4", "CD8", "MAIT", "gdT", "B", "plasma cells", "NK", "NKT", "Th17", "Treg", "Mo","DC", "Neutrophil", "NA")

colmat <- data.frame(color_group_mindre )
colmat$kort_navn <- rownames(colmat)
qTable <- merge(qTable, colmat, by = "kort_navn")

qTable <- qTable[order(qTable$kort_navn),]
qTable$navn_prosent <- paste(qTable$navn, qTable$percent)
colmat2 <- qTable[, c("nnavn_prosent", "color_group_mindre")]
colmat2 <- unique(colmat2)

color_group <- colmat2$color_group_mindre
names(color_group) <- colmat2$navn_prosent


#qTable should now look something like this
# > head(qTable)
# kort_navn  navn markers             clusters        q5      q10      q25       q50       q75       q90       q95
# 1         B  10_7     CD8 kluster 10_7 (8.62%) 0.0000000 0.000000 0.000000 0.0000000 0.3085980 0.6445363 0.8594836
# 2         B  10_7   TCRgd kluster 10_7 (8.62%) 0.0000000 0.000000 0.000000 0.0000000 0.0000000 0.1170456 0.2491127
# 3         B 30_13     CD5 kluster 30_13 (0.3%) 0.0000000 0.000000 0.000000 0.1320241 0.4389344 0.7889156 1.0131819
# 4         B  10_7     IgD kluster 10_7 (8.62%) 0.0000000 0.000000 1.513912 2.9024391 3.6835968 4.1829866 4.4398471
# 5         B  10_7   CD11c kluster 10_7 (8.62%) 0.0000000 0.000000 0.000000 0.0000000 0.2623267 0.8179904 1.4544024
# 6         B  10_7   CD85j kluster 10_7 (8.62%) 0.3783907 1.012516 2.015136 2.6912822 3.1494285 3.4896913 3.6859091
# percent        low     high color_group_mindre
# 1    8.62 1.00000000 3.327063            #FDBF6F
# 2    8.62 1.71768555       NA            #FDBF6F
# 3    0.30 1.67675889       NA            #FDBF6F
# 4    8.62 2.80059886       NA            #FDBF6F
# 5    8.62 0.02229601 3.222714            #FDBF6F
# 6    8.62 1.82053191       NA            #FDBF6F
#if it do not contain the columns low it will not make vertical lines.


mplot <- markerplot_fromQuantileTable(qTable = qTable, colnamesCluster = "navn_prosent", color_group = color_group, color_group_mindre = color_group_mindre)
  

mplot$g
mplot$gleg
  