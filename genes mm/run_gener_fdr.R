scriptPath <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig")
# gener regresjon ----

params <- list()
params$data <- "dGener.csv"
params$fra <- "AIRE"
params$til <- "ZNF532"
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "gen_result_fdr")
params$adj_p <- 0.05
params$adj_p_methods <-  "fdr"
params$utFil <- "gene_analyse_T2_severe_mot_mod.csv"
params$tidspunkt <- TRUE  # vil sammenlikne med 2 level av Status_tid!
params$outliers <- TRUE
#params$subgr <- "H3_control"
params$nivaaer <-  c("Severe T2", "Moderate T2",  "Moderate T1", "Severe T1", "Control")
params$nivaaerSex <-  c("Severe T2 M",  "Severe T2 F",  "Moderate T2 M", "Moderate T2 F",  "Severe T1 M", "Severe T1 F", "Moderate T1 M",  "Moderate T1 F", "Control M", "Control F")
rmarkdown::render(fs:::path(scriptPath, "gen_analyse.Rmd"), output_file = fs::path(params$utSti, gsub(".csv", ".docx", params$utFil)))



params <- list()
params$data <- "dGener.csv"
params$fra <- "AIRE"
params$til <- "ZNF532"
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "gen_result_fdr")
params$adj_p <- 0.05
params$adj_p_methods <-  "fdr"
params$utFil <- "gene_analyse_T1_severe_mot_mod.csv"
params$tidspunkt <- TRUE  # vil sammenlikne med 2 level av Status_tid!
params$outliers <- TRUE
#params$subgr <- "H3_control"
params$nivaaer <-  c("Severe T1", "Moderate T1", "Severe T2", "Moderate T2", "Control")
params$nivaaerSex <-  c( "Severe T1 M", "Severe T1 F", "Moderate T1 M",  "Moderate T1 F", "Severe T2 M",  "Severe T2 F",  "Moderate T2 M", "Moderate T2 F", "Control M", "Control F")
rmarkdown::render(fs:::path(scriptPath, "gen_analyse.Rmd"), output_file = fs::path(params$utSti, gsub(".csv", ".docx", params$utFil)))





params <- list()
params$data <- "dGener.csv"
params$fra <- "AIRE"
params$til <- "ZNF532"
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "gen_result_fdr")
params$adj_p <- 0.05
params$adj_p_methods <-  "fdr"
params$utFil <- "gene_analyse_severeT1_mot_control.csv"
params$tidspunkt <- TRUE  # vil sammenlikne med 2 level av Status_tid!
params$outliers <- TRUE
#params$subgr <- "H3_control"
params$nivaaer <-  c("Severe T1", "Control", "Moderate T1", "Severe T2", "Moderate T2")
params$nivaaerSex <-  c( "Severe T1 M", "Severe T1 F", "Control M", "Control F", "Moderate T1 M",  "Moderate T1 F", "Severe T2 M",  "Severe T2 F",  "Moderate T2 M", "Moderate T2 F")
rmarkdown::render(fs:::path(scriptPath, "gen_analyse.Rmd"), output_file = fs::path(params$utSti, gsub(".csv", ".docx", params$utFil)))




params <- list()
params$data <- "dGener.csv"
params$fra <- "AIRE"
params$til <- "ZNF532"
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "gen_result_fdr")
params$adj_p <- 0.05
params$adj_p_methods <-  "fdr"
params$utFil <- "gene_analyse_moderateT1_mot_control.csv"
params$tidspunkt <- TRUE  # vil sammenlikne med 2 level av Status_tid!
params$outliers <- TRUE
#params$subgr <- "H3_control"
params$nivaaer <-  c( "Moderate T1", "Control", "Severe T1","Severe T2", "Moderate T2")
params$nivaaerSex <-  c("Moderate T1 M",  "Moderate T1 F",  "Control M", "Control F", "Severe T1 M", "Severe T1 F", "Severe T2 M",  "Severe T2 F",  "Moderate T2 M", "Moderate T2 F")
rmarkdown::render(fs:::path(scriptPath, "gen_analyse.Rmd"), output_file = fs::path(params$utSti, gsub(".csv", ".docx", params$utFil)))




params <- list()
params$data <- "dGener.csv"
params$fra <- "AIRE"
params$til <- "ZNF532"
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "gen_result_fdr")
params$adj_p <- 0.05
params$adj_p_methods <-  "fdr"
params$utFil <- "gene_analyse_severeT2_mot_control.csv"
params$tidspunkt <- TRUE  # vil sammenlikne med 2 level av Status_tid!
params$outliers <- TRUE
#params$subgr <- "H3_control"
params$nivaaer <-  c("Severe T2", "Control", "Moderate T2", "Severe T1", "Moderate T1")
params$nivaaerSex <-  c( "Severe T2 M", "Severe T2 F", "Control M", "Control F", "Moderate T2 M",  "Moderate T2 F", "Severe T1 M",  "Severe T1 F",  "Moderate T1 M", "Moderate T1 F")
rmarkdown::render(fs:::path(scriptPath, "gen_analyse.Rmd"), output_file = fs::path(params$utSti, gsub(".csv", ".docx", params$utFil)))




params <- list()
params$data <- "dGener.csv"
params$fra <- "AIRE"
params$til <- "ZNF532"
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "gen_result_fdr")
params$adj_p <- 0.05
params$adj_p_methods <-  "fdr"
params$utFil <- "gene_analyse_moderateT1_mot_control.csv"
params$tidspunkt <- TRUE  # vil sammenlikne med 2 level av Status_tid!
params$outliers <- TRUE
#params$subgr <- "H3_control"
params$nivaaer <-  c( "Moderate T2", "Control", "Severe T2","Severe T1", "Moderate T1")
params$nivaaerSex <-  c("Moderate T2 M",  "Moderate T2 F",  "Control M", "Control F", "Severe T2 M", "Severe T2 F", "Severe T1 M",  "Severe T1 F",  "Moderate T1 M", "Moderate T1 F")
rmarkdown::render(fs:::path(scriptPath, "gen_analyse.Rmd"), output_file = fs::path(params$utSti, gsub(".csv", ".docx", params$utFil)))





params <- list()
params$data <- "dGener.csv"
params$fra <- "AIRE"
params$til <- "ZNF532"
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "gen_result_fdr")
params$adj_p <- 0.05
params$adj_p_methods <-  "fdr"
params$utFil <- "gene_analyse_T1_mot_T2.csv"
params$tidspunkt <- TRUE  # vil sammenlikne med 2 level av Status_tid!
params$outliers <- TRUE
#params$subgr <- "H3_control"
rmarkdown::render(fs:::path(scriptPath, "gen_analyse_kunTid.Rmd"), output_file = fs::path(params$utSti, gsub(".csv", ".docx", params$utFil)))



params <- list()
params$data <- "dGener.csv"
params$fra <- "AIRE"
params$til <- "ZNF532"
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "gen_result_fdr")
params$adj_p <- 0.05
params$adj_p_methods <-  "fdr"
params$utFil <- "gene_analyse_T1_mot_Control.csv"
params$tidspunkt <- TRUE  # vil sammenlikne med 2 level av Status_tid!
params$outliers <- TRUE
#params$subgr <- "H3_control"
params$nivaaer <-  c("T1", "Control", "T2")
params$nivaaerSex <-  c( "T1 M", "T1 F", "Control M", "Control F", "T2 M",  "T2 F")
rmarkdown::render(fs:::path(scriptPath, "gen_analyse_kunTid.Rmd"), output_file = fs::path(params$utSti, gsub(".csv", ".docx", params$utFil)))


params <- list()
params$data <- "dGener.csv"
params$fra <- "AIRE"
params$til <- "ZNF532"
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "gen_result_fdr")
params$adj_p <- 0.05
params$adj_p_methods <-  "fdr"
params$utFil <- "gene_analyse_T2_mot_Control.csv"
params$tidspunkt <- TRUE  # vil sammenlikne med 2 level av Status_tid!
params$outliers <- TRUE
#params$subgr <- "H3_control"
params$nivaaer <-  c("T2", "Control", "T1")
params$nivaaerSex <-  c( "T2 M", "T2 F", "Control M", "Control F", "T1 M",  "T1 F")
rmarkdown::render(fs:::path(scriptPath, "gen_analyse_kunTid.Rmd"), output_file = fs::path(params$utSti, gsub(".csv", ".docx", params$utFil)))



#H3 og kontroll
# gener regresjon ----

params <- list()
params$data <- "dGener.csv"
params$fra <- "AIRE"
params$til <- "ZNF532"
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "gen_result_fdr")
params$adj_p <- 0.05
params$adj_p_methods <-  "fdr"
params$utFil <- "gene_analyse_T2_severe_mot_mod_H3.csv"
params$tidspunkt <- TRUE  # vil sammenlikne med 2 level av Status_tid!
params$outliers <- TRUE
params$subgr <- "H3_control"
params$nivaaer <-  c("Severe T2", "Moderate T2",  "Moderate T1", "Severe T1", "Control")
params$nivaaerSex <-  c("Severe T2 M",  "Severe T2 F",  "Moderate T2 M", "Moderate T2 F",  "Severe T1 M", "Severe T1 F", "Moderate T1 M",  "Moderate T1 F", "Control M", "Control F")
rmarkdown::render(fs:::path(scriptPath, "gen_analyse.Rmd"), output_file = fs::path(params$utSti, gsub(".csv", ".docx", params$utFil)))



params <- list()
params$data <- "dGener.csv"
params$fra <- "AIRE"
params$til <- "ZNF532"
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "gen_result_fdr")
params$adj_p <- 0.05
params$adj_p_methods <-  "fdr"
params$utFil <- "gene_analyse_T1_severe_mot_mod_H3.csv"
params$tidspunkt <- TRUE  # vil sammenlikne med 2 level av Status_tid!
params$outliers <- TRUE
params$subgr <- "H3_control"
params$nivaaer <-  c("Severe T1", "Moderate T1", "Severe T2", "Moderate T2", "Control")
params$nivaaerSex <-  c( "Severe T1 M", "Severe T1 F", "Moderate T1 M",  "Moderate T1 F", "Severe T2 M",  "Severe T2 F",  "Moderate T2 M", "Moderate T2 F", "Control M", "Control F")
rmarkdown::render(fs:::path(scriptPath, "gen_analyse.Rmd"), output_file = fs::path(params$utSti, gsub(".csv", ".docx", params$utFil)))





params <- list()
params$data <- "dGener.csv"
params$fra <- "AIRE"
params$til <- "ZNF532"
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "gen_result_fdr")
params$adj_p <- 0.05
params$adj_p_methods <-  "fdr"
params$utFil <- "gene_analyse_severeT1_mot_control_H3.csv"
params$tidspunkt <- TRUE  # vil sammenlikne med 2 level av Status_tid!
params$outliers <- TRUE
params$subgr <- "H3_control"
params$nivaaer <-  c("Severe T1", "Control", "Moderate T1", "Severe T2", "Moderate T2")
params$nivaaerSex <-  c( "Severe T1 M", "Severe T1 F", "Control M", "Control F", "Moderate T1 M",  "Moderate T1 F", "Severe T2 M",  "Severe T2 F",  "Moderate T2 M", "Moderate T2 F")
rmarkdown::render(fs:::path(scriptPath, "gen_analyse.Rmd"), output_file = fs::path(params$utSti, gsub(".csv", ".docx", params$utFil)))




params <- list()
params$data <- "dGener.csv"
params$fra <- "AIRE"
params$til <- "ZNF532"
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "gen_result_fdr")
params$adj_p <- 0.05
params$adj_p_methods <-  "fdr"
params$utFil <- "gene_analyse_moderateT1_mot_control_H3.csv"
params$tidspunkt <- TRUE  # vil sammenlikne med 2 level av Status_tid!
params$outliers <- TRUE
params$subgr <- "H3_control"
params$nivaaer <-  c( "Moderate T1", "Control", "Severe T1","Severe T2", "Moderate T2")
params$nivaaerSex <-  c("Moderate T1 M",  "Moderate T1 F",  "Control M", "Control F", "Severe T1 M", "Severe T1 F", "Severe T2 M",  "Severe T2 F",  "Moderate T2 M", "Moderate T2 F")
rmarkdown::render(fs:::path(scriptPath, "gen_analyse.Rmd"), output_file = fs::path(params$utSti, gsub(".csv", ".docx", params$utFil)))




params <- list()
params$data <- "dGener.csv"
params$fra <- "AIRE"
params$til <- "ZNF532"
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "gen_result_fdr")
params$adj_p <- 0.05
params$adj_p_methods <-  "fdr"
params$utFil <- "gene_analyse_severeT2_mot_control_H3.csv"
params$tidspunkt <- TRUE  # vil sammenlikne med 2 level av Status_tid!
params$outliers <- TRUE
params$subgr <- "H3_control"
params$nivaaer <-  c("Severe T2", "Control", "Moderate T2", "Severe T1", "Moderate T1")
params$nivaaerSex <-  c( "Severe T2 M", "Severe T2 F", "Control M", "Control F", "Moderate T2 M",  "Moderate T2 F", "Severe T1 M",  "Severe T1 F",  "Moderate T1 M", "Moderate T1 F")
rmarkdown::render(fs:::path(scriptPath, "gen_analyse.Rmd"), output_file = fs::path(params$utSti, gsub(".csv", ".docx", params$utFil)))




params <- list()
params$data <- "dGener.csv"
params$fra <- "AIRE"
params$til <- "ZNF532"
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "gen_result_fdr")
params$adj_p <- 0.05
params$adj_p_methods <-  "fdr"
params$utFil <- "gene_analyse_moderateT1_mot_control_H3.csv"
params$tidspunkt <- TRUE  # vil sammenlikne med 2 level av Status_tid!
params$outliers <- TRUE
params$subgr <- "H3_control"
params$nivaaer <-  c( "Moderate T2", "Control", "Severe T2","Severe T1", "Moderate T1")
params$nivaaerSex <-  c("Moderate T2 M",  "Moderate T2 F",  "Control M", "Control F", "Severe T2 M", "Severe T2 F", "Severe T1 M",  "Severe T1 F",  "Moderate T1 M", "Moderate T1 F")
rmarkdown::render(fs:::path(scriptPath, "gen_analyse.Rmd"), output_file = fs::path(params$utSti, gsub(".csv", ".docx", params$utFil)))





params <- list()
params$data <- "dGener.csv"
params$fra <- "AIRE"
params$til <- "ZNF532"
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "gen_result_fdr")
params$adj_p <- 0.05
params$adj_p_methods <-  "fdr"
params$utFil <- "gene_analyse_T1_mot_T2_H3.csv"
params$tidspunkt <- TRUE  # vil sammenlikne med 2 level av Status_tid!
params$outliers <- TRUE
params$subgr <- "H3_control"
rmarkdown::render(fs:::path(scriptPath, "gen_analyse_kunTid.Rmd"), output_file = fs::path(params$utSti, gsub(".csv", ".docx", params$utFil)))



params <- list()
params$data <- "dGener.csv"
params$fra <- "AIRE"
params$til <- "ZNF532"
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "gen_result_fdr")
params$adj_p <- 0.05
params$adj_p_methods <-  "fdr"
params$utFil <- "gene_analyse_T1_mot_Control_H3.csv"
params$tidspunkt <- TRUE  # vil sammenlikne med 2 level av Status_tid!
params$outliers <- TRUE
params$subgr <- "H3_control"
params$nivaaer <-  c("T1", "Control", "T2")
params$nivaaerSex <-  c( "T1 M", "T1 F", "Control M", "Control F", "T2 M",  "T2 F")
rmarkdown::render(fs:::path(scriptPath, "gen_analyse_kunTid.Rmd"), output_file = fs::path(params$utSti, gsub(".csv", ".docx", params$utFil)))


params <- list()
params$data <- "dGener.csv"
params$fra <- "AIRE"
params$til <- "ZNF532"
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "gen_result_fdr")
params$adj_p <- 0.05
params$adj_p_methods <-  "fdr"
params$utFil <- "gene_analyse_T2_mot_Control_H3.csv"
params$tidspunkt <- TRUE  # vil sammenlikne med 2 level av Status_tid!
params$outliers <- TRUE
params$subgr <- "H3_control"
params$nivaaer <-  c("T2", "Control", "T1")
params$nivaaerSex <-  c( "T2 M", "T2 F", "Control M", "Control F", "T1 M",  "T1 F")
rmarkdown::render(fs:::path(scriptPath, "gen_analyse_kunTid.Rmd"), output_file = fs::path(params$utSti, gsub(".csv", ".docx", params$utFil)))


rmarkdown::render(fs:::path(scriptPath, "vulcanoPlot_fdr.Rmd"), output_file = fs::path(scriptPath, "gen_result_fdr", paste0("Result_vulcanoplot.docx")))

