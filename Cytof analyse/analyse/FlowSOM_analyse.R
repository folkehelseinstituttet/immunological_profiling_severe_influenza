set.seed(params$seed) # to ensure same plot every time
out <- FlowSOM::ReadInput(as.matrix(params$data[,params$kanaler]), transform = F, scale = params$scaling)
print(1)
out <- FlowSOM::BuildSOM(out, colsToUse = 1:(ncol(params$data[,params$kanaler])), xdim = params$xdim, ydim = params$ydim)
print(2)
out <- FlowSOM::BuildMST(out)
print(3)
cluster_FlowSOM_pre <- out$map$mapping[, 1]

clusterFraFlowSOM <- matrix(NA, nrow = nrow(params$data), ncol <- length(params$ks) )
colnames(clusterFraFlowSOM) <- paste("kluster", params$ks)

for(i in 1:length(params$ks)){
  k <- params$ks[i]
  print(k)
  set.seed(params$seed)
  out_k <- FlowSOM::metaClustering_consensus(out$map$codes, k = k, seed = params$seed)
  cluster_FlowSOM_k <- out_k[cluster_FlowSOM_pre]
  clusterFraFlowSOM[,i] <- cluster_FlowSOM_k
  cluster_FlowSOM_k_factor <- factor(cluster_FlowSOM_k, levels = 1:k)
print("cluster")
# ved å inkludere disse ble ikke datasettet likt. 
  q5_k <- q_per_cluster_marker(data = params$data[,params$kanaler], kluster = cluster_FlowSOM_k, probs = 0.05)
  write.csv2(q5_k, fs::path(params$utSti, paste0("q5_per_kluster_k_", k, "_seed", params$seed, params$ext_name, ".csv")))
  q10_k <- q_per_cluster_marker(data = params$data[,params$kanaler], kluster = cluster_FlowSOM_k, probs = 0.1)
  write.csv2(q10_k, fs::path(params$utSti, paste0("q10_per_kluster_k_", k, "_seed", params$seed, params$ext_name, ".csv")))
  q25_k <- q_per_cluster_marker(data = params$data[,params$kanaler], kluster = cluster_FlowSOM_k, probs = 0.25)
  write.csv2(q25_k, fs::path(params$utSti, paste0("q25_per_kluster_k_", k, "_seed", params$seed, params$ext_name, ".csv")))
  q75_k <- q_per_cluster_marker(data = params$data[,params$kanaler], kluster = cluster_FlowSOM_k, probs = 0.75)
  write.csv2(q75_k, fs::path(params$utSti, paste0("q75_per_kluster_k_", k, "_seed", params$seed, params$ext_name, ".csv")))
  q90_k <- q_per_cluster_marker(data = params$data[,params$kanaler], kluster = cluster_FlowSOM_k, probs = 0.9)
  write.csv2(q90_k, fs::path(params$utSti, paste0("q90_per_kluster_k_", k, "_seed", params$seed, params$ext_name, ".csv")))
  q95_k <- q_per_cluster_marker(data = params$data[,params$kanaler], kluster = cluster_FlowSOM_k, probs = 0.95)
  write.csv2(q95_k, fs::path(params$utSti, paste0("q95_per_kluster_k_", k, "_seed", params$seed, params$ext_name, ".csv")))
  
print("q")

  medians_k <- median_per_cluster_marker(data = params$data[,params$kanaler], kluster = cluster_FlowSOM_k)
  write.csv2(medians_k, fs::path(params$utSti, paste0("medians_per_kluster_k_", k, "_seed", params$seed, params$ext_name, ".csv")))
print("median")



  tiff(fs::path(params$utSti, paste0("heatmap_median_k_", k, "_cluster_seed", params$seed, params$ext_name, ".tiff")), width = 1000, height = 800)
    print(heatmap(as.matrix(medians_k[,params$kanaler], cluster_columns = params$column_cluster)))
  dev.off()
print("heatmap")


  totalt_antall_i_kluster <- table(cluster_FlowSOM_k_factor)
  

  antallPerKluster <- cbind(totalt_antall_i_kluster, totalt_antall_i_kluster/nrow(params$data))
  colnames(antallPerKluster) <- c("antall", "prosent")
  write.csv2(antallPerKluster, fs::path(params$utSti, paste0("antallPerCluster_k_", k, "_seed", params$seed, params$ext_name, ".csv")))
print("totalt antall")

  
  perSample <- table(params$data$dataset, cluster_FlowSOM_k_factor)
  colnames(perSample) <- paste0("kluster_", 1:k)

  write.csv2(perSample, fs::path(params$utSti, paste0("antallPerClusterOgPat_k_", k, "_seed", params$seed, params$ext_name, ".csv")))
print("antall per fil")



}

#dataBrukt <- cbind(params$data[,params$kanaler], clusterFraFlowSOM)
dataBrukt <- cbind(params$data, clusterFraFlowSOM)
saveRDS(dataBrukt, fs::path(params$utSti, paste0("data", params$seed, ".RDS")))
