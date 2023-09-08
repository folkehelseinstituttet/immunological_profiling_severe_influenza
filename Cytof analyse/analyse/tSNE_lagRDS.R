posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel1_mars2022", "posNeg", "Data")
kort_fil_navn <- read.csv2(fs::path(posNeg_path, "kortFilNavn.csv"))
source("F:/Forskningsprosjekter/PDB 2794 - Immune responses aga_/Forskningsfiler/JOBO/CyTOF/Analyse i R OUS/CleanUpGatingMarch2022/Analyse/ploting_functions.R")


unsup_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "Cytof unsupperviced")
panel <- "Panel 1"
subData <- "All"
seed <- "seed 1345"
posNeg <- readRDS( fs::path(unsup_path, panel, subData, seed, paste0("posNeg_seed_", gsub("seed ", "", seed), ".RDS")))

ut_sti <- fs::path(unsup_path, panel, subData, seed)
dAlle <- readRDS( fs::path(unsup_path, panel, subData, seed, paste0("data_med_filnavn", gsub("seed ", "", seed), ".RDS")))
dAlle <- merge(dAlle, kort_fil_navn, by.x = "dataset", by.y = "filnavn")

set.seed(100) # to ensure same plot every time
temp <- random_events_vector(dAlle$dataset, 5000)

d <- dAlle[temp,]
posNegMat_mindre <- posNegMat[temp,]


d_signal <- d[,!grepl("kluster", colnames(d))]
d_signal <- d_signal[,!grepl("dataset", colnames(d_signal))]
d_signal <- d_signal[,!grepl("kort_filnavn", colnames(d_signal))]
d_signal <- d_signal[,!grepl("status", colnames(d_signal))]
d_signal <- d_signal[,!grepl("age", colnames(d_signal))]
d_signal <- d_signal[,!grepl("sex", colnames(d_signal))]

d <- merge(d, kort_fil_navn, by.x = "dataset", by.y = "filnavn")

d$statusTid <- d$status
set.seed(100) # to ensure same plot every time
umapAshinfac5 <- umap::umap(d_signal)


library(Rtsne)

set.seed(100) # to ensure same plot every time
res_tsne <- Rtsne(d_signal)

saveRDS(list(d = d,  posNegMat = posNegMat_mindre, umap = umapAshinfac5$layout,  res_tsne = res_tsne$Y), fs::path(ut_sti, "dataUMAP_tSNE_5000_nr2.RDS"))  #







##panel 2

posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel2_mars2022", "posNeg", "Data")
kort_fil_navn <- read.csv2(fs::path(posNeg_path, "kortFilNavn.csv"))
source("F:/Forskningsprosjekter/PDB 2794 - Immune responses aga_/Forskningsfiler/JOBO/CyTOF/Analyse i R OUS/CleanUpGatingMarch2022/Analyse/ploting_functions.R")


unsup_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "Cytof unsupperviced")
panel <- "Panel 2"
subData <- "All"
seed <- "seed 2234"

ut_sti <- fs::path(unsup_path, panel, subData, seed)
dAlle <- readRDS( fs::path(unsup_path, panel, subData, seed, paste0("data_med_filnavn", gsub("seed ", "", seed), ".RDS")))
dAlle <- merge(dAlle, kort_fil_navn, by.x = "dataset", by.y = "filnavn")
temp <- random_events_vector(dAlle$dataset, 5000)

d <- dAlle[temp,]


d_signal <- d[,!grepl("kluster", colnames(d))]
d_signal <- d_signal[,!grepl("dataset", colnames(d_signal))]
d_signal <- d_signal[,!grepl("kort_filnavn", colnames(d_signal))]
d_signal <- d_signal[,!grepl("status", colnames(d_signal))]
d_signal <- d_signal[,!grepl("age", colnames(d_signal))]
d_signal <- d_signal[,!grepl("sex", colnames(d_signal))]
d_signal <- d_signal[,!grepl("X", colnames(d_signal))]
d_signal <- d_signal[,!grepl("X.1", colnames(d_signal))]

d <- merge(d, kort_fil_navn, by.x = "dataset", by.y = "filnavn")

d$statusTid <- d$status
# set.seed(100) # to ensure same plot every time
# umapAshinfac5 <- umap::umap(d_signal)
# 

library(Rtsne)

set.seed(100) # to ensure same plot every time
res_tsne <- Rtsne(d_signal)

saveRDS(list(d = d,   res_tsne = res_tsne$Y), fs::path(ut_sti, "dataUMAP_tSNE_5000.RDS"))  #umap = umapAshinfac5$layout,






