
#Manuel_manuel P1

scriptPath_UnsupAnalysis <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "Cytof unsupperviced")
scriptPath_UnsupAnalysisNy <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "Cytof unsupperviced", "Endelig_des2022")

params <- list()
params$panel <- "Panel 1"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c("Severe T1", "Moderate T1", "Severe T2", "Moderate T2", "Control" )
params$dataSti <- fs::path(scriptPath_UnsupAnalysisNy, "Figurer", "Manuel_ManuelGating", "P1_cells_per_manuel_cluster_and_sample.csv")
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, "Figurer", "Manuel_ManuelGating")
params$klustre <- c("B", "CD4", "CD8", "DC", "MAIT", "Mo", "NK", "NKT", "gdT")
params$test <- "S1vsM1"
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_manuel_manuel.Rmd"), output_file = fs::path(params$resultSti, paste0("P", gsub("Panel ", "", params$panel), "S1vsM1.docx")))





scriptPath_UnsupAnalysis <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "Cytof unsupperviced")
scriptPath_UnsupAnalysisNy <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "Cytof unsupperviced", "Endelig_des2022")

params <- list()
params$panel <- "Panel 1"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c("Severe T1", "Control" , "Moderate T1", "Severe T2", "Moderate T2")
params$dataSti <- fs::path(scriptPath_UnsupAnalysisNy, "Figurer", "Manuel_ManuelGating", "P1_cells_per_manuel_cluster_and_sample.csv")
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, "Figurer", "Manuel_ManuelGating")
params$klustre <- c("B", "CD4", "CD8", "DC", "MAIT", "Mo", "NK", "NKT", "gdT")
params$test <- "S1vsC"
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_manuel_manuel.Rmd"), output_file = fs::path(params$resultSti, paste0("P", gsub("Panel ", "", params$panel), "S1vsC.docx")))




scriptPath_UnsupAnalysis <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "Cytof unsupperviced")
scriptPath_UnsupAnalysisNy <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "Cytof unsupperviced", "Endelig_des2022")

params <- list()
params$panel <- "Panel 1"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c( "Moderate T1", "Control" , "Severe T1", "Severe T2", "Moderate T2")
params$dataSti <- fs::path(scriptPath_UnsupAnalysisNy, "Figurer", "Manuel_ManuelGating", "P1_cells_per_manuel_cluster_and_sample.csv")
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, "Figurer", "Manuel_ManuelGating")
params$klustre <- c("B", "CD4", "CD8", "DC", "MAIT", "Mo", "NK", "NKT", "gdT")
params$test <- "M1vsC"
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_manuel_manuel.Rmd"), output_file = fs::path(params$resultSti, paste0("P", gsub("Panel ", "", params$panel), "M1vsC.docx")))


scriptPath_UnsupAnalysis <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "Cytof unsupperviced")
scriptPath_UnsupAnalysisNy <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "Cytof unsupperviced", "Endelig_des2022")

params <- list()
params$panel <- "Panel 1"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c("Severe T2", "Moderate T2", "Severe T1", "Moderate T1", "Control" )
params$dataSti <- fs::path(scriptPath_UnsupAnalysisNy, "Figurer", "Manuel_ManuelGating", "P1_cells_per_manuel_cluster_and_sample.csv")
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, "Figurer", "Manuel_ManuelGating")
params$klustre <- c("B", "CD4", "CD8", "DC", "MAIT", "Mo", "NK", "NKT", "gdT")
params$test <- "S2vsM2"
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_manuel_manuel.Rmd"), output_file = fs::path(params$resultSti, paste0("P", gsub("Panel ", "", params$panel), "S2vsM2.docx")))





scriptPath_UnsupAnalysis <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "Cytof unsupperviced")
scriptPath_UnsupAnalysisNy <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "Cytof unsupperviced", "Endelig_des2022")

params <- list()
params$panel <- "Panel 1"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c("Severe T2", "Control" , "Moderate T2", "Severe T1", "Moderate T1")
params$dataSti <- fs::path(scriptPath_UnsupAnalysisNy, "Figurer", "Manuel_ManuelGating", "P1_cells_per_manuel_cluster_and_sample.csv")
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, "Figurer", "Manuel_ManuelGating")
params$klustre <- c("B", "CD4", "CD8", "DC", "MAIT", "Mo", "NK", "NKT", "gdT")
params$test <- "S2vsC"
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_manuel_manuel.Rmd"), output_file = fs::path(params$resultSti, paste0("P", gsub("Panel ", "", params$panel), "S2vsC.docx")))




scriptPath_UnsupAnalysis <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "Cytof unsupperviced")
scriptPath_UnsupAnalysisNy <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "Cytof unsupperviced", "Endelig_des2022")

params <- list()
params$panel <- "Panel 1"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c( "Moderate T2", "Control" , "Severe T1", "Severe T2", "Moderate T1")
params$dataSti <- fs::path(scriptPath_UnsupAnalysisNy, "Figurer", "Manuel_ManuelGating", "P1_cells_per_manuel_cluster_and_sample.csv")
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, "Figurer", "Manuel_ManuelGating")
params$klustre <- c("B", "CD4", "CD8", "DC", "MAIT", "Mo", "NK", "NKT", "gdT")
params$test <- "M2vsC"
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_manuel_manuel.Rmd"), output_file = fs::path(params$resultSti, paste0("P", gsub("Panel ", "", params$panel), "M2vsC.docx")))







#panel 1 alle ----

scriptPath_UnsupAnalysis <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "Cytof unsupperviced")
scriptPath_UnsupAnalysisNy <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "Cytof unsupperviced", "Endelig_des2022")

params <- list()
params$panel <- "Panel 1"
params$seed <- 1234 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c("Severe T1", "Moderate T1", "Severe T2", "Moderate T2", "Control" )
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, paste0("Result_", params$seed, "ST1motMT1.docx")))



params <- list()
params$panel <- "Panel 1"
params$seed <- 1345 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c("Severe T1", "Moderate T1", "Severe T2", "Moderate T2", "Control" )
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, paste0("Result_", params$seed, "ST1motMT1.docx")))

params <- list()
params$panel <- "Panel 1"
params$seed <- 1456 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c("Severe T1", "Moderate T1", "Severe T2", "Moderate T2", "Control" )
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, paste0("Result_", params$seed, "ST1motMT1.docx")))

params <- list()
params$panel <- "Panel 1"
params$seed <- 1567 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$utSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c("Severe T1", "Moderate T1", "Severe T2", "Moderate T2", "Control" )
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, paste0("Result_", params$seed, "ST1motMT1.docx")))




#panel 1 sev T1 mot kontrol

#panel 1 alle ----

params <- list()
params$panel <- "Panel 1"
params$seed <- 1234 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c("Severe T1", "Control", "Moderate T1", "Severe T2", "Moderate T2")
params$ekstraSti <-  "ST1vsCont"
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, "ST1vsCont", paste0("Res.docx")))




params <- list()
params$panel <- "Panel 1"
params$seed <- 1345 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c("Severe T1", "Control", "Moderate T1", "Severe T2", "Moderate T2" )
params$ekstraSti <-  "ST1vsCont"
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, "ST1vsCont", paste0("Res.docx")))

params <- list()
params$panel <- "Panel 1"
params$seed <- 1456 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c("Severe T1", "Control", "Moderate T1", "Severe T2", "Moderate T2" )
params$ekstraSti <-  "ST1vsCont"
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, "ST1vsCont", paste0("Res.docx")))

params <- list()
params$panel <- "Panel 1"
params$seed <- 1567 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c("Severe T1", "Control", "Moderate T1", "Severe T2", "Moderate T2" )
params$ekstraSti <-  "ST1vsCont"
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, "ST1vsCont", paste0("Res.docx")))




#panel 1 mod T1 mot kontrol

#panel 1 alle ----

params <- list()
params$panel <- "Panel 1"
params$seed <- 1234 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c( "Moderate T1", "Control", "Severe T1", "Severe T2", "Moderate T2" )
params$ekstraSti <-  "MT1vsCont"
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, "MT1vsCont", paste0("Res.docx")))



params <- list()
params$panel <- "Panel 1"
params$seed <- 1345 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c( "Moderate T1", "Control", "Severe T1", "Severe T2", "Moderate T2" )
params$ekstraSti <-  "MT1vsCont"
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, "MT1vsCont", paste0("Res.docx")))

params <- list()
params$panel <- "Panel 1"
params$seed <- 1456 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c( "Moderate T1", "Control", "Severe T1", "Severe T2", "Moderate T2" )
params$ekstraSti <-  "MT1vsCont"
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, "MT1vsCont", paste0("Res.docx")))

params <- list()
params$panel <- "Panel 1"
params$seed <- 1567 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c( "Moderate T1", "Control", "Severe T1", "Severe T2", "Moderate T2" )
params$ekstraSti <-  "MT1vsCont"
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, "MT1vsCont", paste0("Res.docx")))









#panel 2 alle ----

params <- list()
params$panel <- "Panel 2"
params$seed <- 2234 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c("Severe T1", "Moderate T1", "Severe T2", "Moderate T2", "Control" )
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, paste0("Result_", params$seed, "ST1motMT1.docx")))




params <- list()
params$panel <- "Panel 2"
params$seed <- 2345 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c("Severe T1", "Moderate T1", "Severe T2", "Moderate T2", "Control" )
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, paste0("Result_", params$seed, "ST1motMT1.docx")))


params <- list()
params$panel <- "Panel 2"
params$seed <- 2456 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c("Severe T1", "Moderate T1", "Severe T2", "Moderate T2", "Control" )
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, paste0("Result_", params$seed, "ST1motMT1.docx")))



params <- list()
params$panel <- "Panel 2"
params$seed <- 2567 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c("Severe T1", "Moderate T1", "Severe T2", "Moderate T2", "Control" )
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, paste0("Result_", params$seed, "ST1motMT1.docx")))






#Panel 2 sev T1 mot kontrol

#Panel 2 alle ----


params <- list()
params$panel <- "Panel 2"
params$seed <- 2234 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c("Severe T1", "Control", "Moderate T1", "Severe T2", "Moderate T2")
params$ekstraSti <-  "ST1vsCont"
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, "ST1vsCont", paste0("Res.docx")))




params <- list()
params$panel <- "Panel 2"
params$seed <- 2345 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c("Severe T1", "Control", "Moderate T1", "Severe T2", "Moderate T2")
params$ekstraSti <-  "ST1vsCont"
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, "ST1vsCont", paste0("Res.docx")))

params <- list()
params$panel <- "Panel 2"
params$seed <- 2456 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c("Severe T1", "Control", "Moderate T1", "Severe T2", "Moderate T2")
params$ekstraSti <-  "ST1vsCont"
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, "ST1vsCont", paste0("Res.docx")))

params <- list()
params$panel <- "Panel 2"
params$seed <- 2567 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c("Severe T1", "Control", "Moderate T1", "Severe T2", "Moderate T2")
params$ekstraSti <-  "ST1vsCont"
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, "ST1vsCont", paste0("Res.docx")))



#Panel 2 mod T1 mot kontrol

#Panel 2 alle ----


params <- list()
params$panel <- "Panel 2"
params$seed <- 2234 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c( "Moderate T1", "Control", "Severe T1", "Severe T2", "Moderate T2" )
params$ekstraSti <-  "MT1vsCont"
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, "MT1vsCont", paste0("Res.docx")))



params <- list()
params$panel <- "Panel 2"
params$seed <- 2345 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c( "Moderate T1", "Control", "Severe T1", "Severe T2", "Moderate T2" )
params$ekstraSti <-  "MT1vsCont"
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, "MT1vsCont", paste0("Res.docx")))

params <- list()
params$panel <- "Panel 2"
params$seed <- 2456 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c( "Moderate T1", "Control", "Severe T1", "Severe T2", "Moderate T2" )
params$ekstraSti <-  "MT1vsCont"
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, "MT1vsCont", paste0("Res.docx")))

params <- list()
params$panel <- "Panel 2"
params$seed <- 2567 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c( "Moderate T1", "Control", "Severe T1", "Severe T2", "Moderate T2" )
params$ekstraSti <-  "MT1vsCont"
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, "MT1vsCont", paste0("Res.docx")))




#ST2 MT2 kont



scriptPath_UnsupAnalysis <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "Cytof unsupperviced")
scriptPath_UnsupAnalysisNy <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "Cytof unsupperviced", "Endelig_des2022")

params <- list()
params$panel <- "Panel 1"
params$seed <- 1234 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c("Severe T2", "Moderate T2", "Severe T1", "Moderate T1", "Control" )
params$ekstraSti <-  "ST2vsMT2"
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, "ST2vsMT2", paste0("Res.docx")))



params <- list()
params$panel <- "Panel 1"
params$seed <- 1345 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c("Severe T2", "Moderate T2", "Severe T1", "Moderate T1", "Control" )
params$ekstraSti <-  "ST2vsMT2"
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, "ST2vsMT2", paste0("Res.docx")))

params <- list()
params$panel <- "Panel 1"
params$seed <- 1456 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c("Severe T2", "Moderate T2", "Severe T1", "Moderate T1", "Control" )
params$ekstraSti <-  "ST2vsMT2"
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, "ST2vsMT2", paste0("Res.docx")))

params <- list()
params$panel <- "Panel 1"
params$seed <- 1567 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$utSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c("Severe T2", "Moderate T2", "Severe T1", "Moderate T1", "Control" )
params$ekstraSti <-  "ST2vsMT2"
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, "ST2vsMT2", paste0("Res.docx")))




#panel 1 sev T2 mot kontrol

#panel 1 alle ----

params <- list()
params$panel <- "Panel 1"
params$seed <- 1234 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c("Severe T2", "Control", "Moderate T2", "Severe T1", "Moderate T1")
params$ekstraSti <-  "ST2vsCont"
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, "ST2vsCont", paste0("Res.docx")))




params <- list()
params$panel <- "Panel 1"
params$seed <- 1345 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c("Severe T2", "Control", "Moderate T2", "Severe T1", "Moderate T1")
params$ekstraSti <-  "ST2vsCont"
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, "ST2vsCont", paste0("Res.docx")))

params <- list()
params$panel <- "Panel 1"
params$seed <- 1456 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c("Severe T2", "Control", "Moderate T2", "Severe T1", "Moderate T1")
params$ekstraSti <-  "ST2vsCont"
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, "ST2vsCont", paste0("Res.docx")))

params <- list()
params$panel <- "Panel 1"
params$seed <- 1567 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c("Severe T2", "Control", "Moderate T2", "Severe T1", "Moderate T1")
params$ekstraSti <-  "ST2vsCont"
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, "ST2vsCont", paste0("Res.docx")))




#panel 1 mod T1 mot kontrol

#panel 1 alle ----

params <- list()
params$panel <- "Panel 1"
params$seed <- 1234 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c("Moderate T2", "Control", "Severe T2", "Severe T1", "Moderate T1")
params$ekstraSti <-  "MT2vsCont"
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, "MT2vsCont", paste0("Res.docx")))



params <- list()
params$panel <- "Panel 1"
params$seed <- 1345 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c("Moderate T2", "Control", "Severe T2", "Severe T1", "Moderate T1")
params$ekstraSti <-  "MT2vsCont"
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, "MT2vsCont", paste0("Res.docx")))

params <- list()
params$panel <- "Panel 1"
params$seed <- 1456 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c("Moderate T2", "Control", "Severe T2", "Severe T1", "Moderate T1")
params$ekstraSti <-  "MT2vsCont"
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, "MT2vsCont", paste0("Res.docx")))

params <- list()
params$panel <- "Panel 1"
params$seed <- 1567 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c("Moderate T2", "Control", "Severe T2", "Severe T1", "Moderate T1")
params$ekstraSti <-  "MT2vsCont"
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, "MT2vsCont", paste0("Res.docx")))









#panel 2 alle ----

params <- list()
params$panel <- "Panel 2"
params$seed <- 2234 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c("Severe T2", "Moderate T2", "Severe T1", "Moderate T1", "Control" )
params$ekstraSti <-  "ST2vsMT2"
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, "ST2vsMT2", paste0("Res.docx")))




params <- list()
params$panel <- "Panel 2"
params$seed <- 2345 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c("Severe T2", "Moderate T2", "Severe T1", "Moderate T1", "Control" )
params$ekstraSti <-  "ST2vsMT2"
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, "ST2vsMT2", paste0("Res.docx")))


params <- list()
params$panel <- "Panel 2"
params$seed <- 2456 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c("Severe T2", "Moderate T2", "Severe T1", "Moderate T1", "Control" )
params$ekstraSti <-  "ST2vsMT2"
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, "ST2vsMT2", paste0("Res.docx")))



params <- list()
params$panel <- "Panel 2"
params$seed <- 2567 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c("Severe T2", "Moderate T2", "Severe T1", "Moderate T1", "Control" )
params$ekstraSti <-  "ST2vsMT2"
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, "ST2vsMT2", paste0("Res.docx")))






#Panel 2 sev T2 mot kontrol

#Panel 2 alle ----


params <- list()
params$panel <- "Panel 2"
params$seed <- 2234 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c("Severe T2", "Control", "Moderate T2", "Severe T1", "Moderate T1")
params$ekstraSti <-  "ST2vsCont"
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, "ST2vsCont", paste0("Res.docx")))




params <- list()
params$panel <- "Panel 2"
params$seed <- 2345 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c("Severe T2", "Control", "Moderate T2", "Severe T1", "Moderate T1")
params$ekstraSti <-  "ST2vsCont"
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, "ST2vsCont", paste0("Res.docx")))

params <- list()
params$panel <- "Panel 2"
params$seed <- 2456 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c("Severe T2", "Control", "Moderate T2", "Severe T1", "Moderate T1")
params$ekstraSti <-  "ST2vsCont"
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, "ST2vsCont", paste0("Res.docx")))

params <- list()
params$panel <- "Panel 2"
params$seed <- 2567 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c("Severe T2", "Control", "Moderate T2", "Severe T1", "Moderate T1")
params$ekstraSti <-  "ST2vsCont"
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, "ST2vsCont", paste0("Res.docx")))



#Panel 2 mod T2 mot kontrol

#Panel 2 alle ----


params <- list()
params$panel <- "Panel 2"
params$seed <- 2234 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c("Moderate T2", "Control", "Severe T2", "Severe T1", "Moderate T1")
params$ekstraSti <-  "MT2vsCont"
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, "MT2vsCont", paste0("Res.docx")))



params <- list()
params$panel <- "Panel 2"
params$seed <- 2345 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c("Moderate T2", "Control", "Severe T2", "Severe T1", "Moderate T1")
params$ekstraSti <-  "MT2vsCont"
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, "MT2vsCont", paste0("Res.docx")))

params <- list()
params$panel <- "Panel 2"
params$seed <- 2456 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c("Moderate T2", "Control", "Severe T2", "Severe T1", "Moderate T1")
params$ekstraSti <-  "MT2vsCont"
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, "MT2vsCont", paste0("Res.docx")))

params <- list()
params$panel <- "Panel 2"
params$seed <- 2567 #nb må endres vil man vil gjøre et annet uttrekk
params$ks <- c(10, 20, 30 ,40, 50, 60)
params$n_per_file <- 25000
params$ext_name <- "All"
params$adj_p <- 0.05
params$adj_p_methods <- "fdr" #"bonferroni" # "fdr" # see help(p.adjust) for other methods
params$tidspunkt <- TRUE
params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
params$nivaa <- c("Moderate T2", "Control", "Severe T2", "Severe T1", "Moderate T1")
params$ekstraSti <-  "MT2vsCont"
params$dataSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
params$resultSti <-  fs::path(scriptPath_UnsupAnalysisNy, params$panel, params$ext_name, paste("seed", params$seed))
rmarkdown::render(fs:::path(scriptPath_UnsupAnalysisNy , "resultater_Regression_negbin_fraFlowSOM.Rmd"), output_file = fs::path(params$resultSti, "MT2vsCont", paste0("Res.docx")))

