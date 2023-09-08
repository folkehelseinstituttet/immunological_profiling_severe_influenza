

# d: data frame with one row per cell, coloum

markerplot <- function(d, markers = "all", colnamesCluster = "Clusters"  , clusters = "all", color_clusters = NA, with_percent = FALSE, gates = NULL){
    
  
  q5to95 <- function(x){
    quantile(x, c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95))
  }
   
  n <- nrow(d)
     
    
    if(colnamesCluster %in% colnames(d)){
      i_clusters <- which(colnames(d) == colnamesCluster)
    } else {
      stop("colnamesClusters has to be the name of the cluster")
    } 
    
    if(!clusters[1] == "all"){
      d <- d[d[,i_clusters] %in% clusters,]
      d[,i_clusters] <- factor(d[,i_clusters], levels = clusters)
    } else {
      clusters <- sort(unique(d[,i])) 
    }

    
    if(is.na(color_clusters)[1]){
      color_clusters <- rep("grey", length(clusters))
      names( color_clusters) <- clusters
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
    }
    
    if(!is.null(gates)[1]){
      tmp <- merge(tmp, gates, by = "markers")
    }
    
    
    
    g <- ggplot(tmp, aes(x = q50, y = clusters, xmin = q10, xmax = q90, col = clusters)) +
      geom_point(size = 2) +
      geom_errorbar() +
      geom_errorbar(data = tmp, aes(x = q50, y = clusters, xmin = q5, xmax = q95, col = clusters))+
      geom_errorbar(data = tmp, aes(x = q50, y = clusters, xmin = q25, xmax = q75, col = clusters))+
      facet_wrap(vars(markers), ncol = 11, scale = "free_x") + theme_classic(base_size = 18) +
      scale_color_manual(values = color_clusters) + xlab("") + ylab("") +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + theme(axis.title = element_text(face="bold")) +  theme(legend.position = "none")
    
    
    if(!is.null(gates)){
      g <- g + geom_vline(aes(xintercept  = low), linetype = "dashed", linewidth = 1.1) +
        geom_vline(aes(xintercept  = high), linetype = "dashed", linewidth = 1.1)
    }
    
    return(g)
}



### user example


dataSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "Cytof unsupperviced", "Panel 1", "All", "seed 1345")
d <- readRDS(fs::path(dataSti, "data_med_filnavn1345.RDS"))

colnamesCluster <- "kluster 60"
clusters <- c("1", "3","5","10","55") #name of the clusters you want to have in the marker plot. If all you might just write "all" defalut. 
color_clusters <- c("1" ="#33CCCC","3" = "#1F78B4", "5" = "#E31A1C", "10" = "#FB9A99", "55" = "#FDBF6F") #default is all grey. 


markers <- colnames(d)[1:43] #vector of the markers you want in the marker plot with same name as column names in d

#gates a matrix with the columns: markers, low and high (is only made for two possible gates, can have only low). if included will give vertical lines in the plot

markerplot(d = d, markers = markers, colnamesCluster = "kluster 60", clusters = clusters, color_clusters = color_clusters, with_percent = TRUE, gates = NULL)
