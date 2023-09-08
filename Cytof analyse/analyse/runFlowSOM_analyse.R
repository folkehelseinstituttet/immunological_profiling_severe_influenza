rm(list = ls())

library(ComplexHeatmap)
library(gridExtra)
library(ggplot2)
library(grid)


orderPanel1 <- c("89Y_CD45", "116Cd_CD3", "145Nd_CD4", "113Cd_CD8", "111Cd_CD19", "152Sm_TCRgd", "166Er_TCRVa7.2", "149Sm_CD25",
                 "158Gd_CD27", "160Gd_CD28", "143Nd_CD127", "172Yb_CD38", "167Er_CCR7", "155Gd_CD45RA", "114Cd_HLADR", "146Nd_IgD", "159Tb_IgG",
                 "156Gd_CXCR3", "153Eu_CCR4", "141Pr_CCR6", "171Yb_CXCR5", "168Er_ICOS", "142Nd_KLRG1", "150Nd_CD134_OX40", "154Sm_TIGIT",
                 "161Dy_CD160", "164Dy_CD161", "162Dy_CD95", "163Dy_CRTH2", "165Ho_CD85j", "169Tm_NKG2A", "173Yb_CD141", "174Yb_CD279_PD-1",
                 "175Lu_CD14", "148Nd_CD16", "176Yb_CD56", "106Cd_CD57", "209Bi_CD11b", "147Sm_CD11c", "112Cd_CD5", "144Nd_CD15", "170Er_CD169", "151Eu_CD123")




scriptPath_UnsupAnalysis <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "Cytof unsupperviced")

source(fs::path(scriptPath_UnsupAnalysis, "readDataToAnalysePanel1TighterCleanUp.R"))
print("data read")


tamed <- as.character(dInfo$filnavn)
tamed <- tamed[!grepl("FHI005", tamed)] #denne oppfører seg helt rart og tas ut av analysen.
tamed <- tamed[!grepl("Ref", tamed)]

#seed 1234 panel 1 alle

params <- list()
params$panel <- "Panel 1"
params$seed <- 1234
params$ks <- c(10, 20, 30, 40, 50, 60)
params$n_per_file <- 25000
params$xdim <- 14
params$ydim <- 14
params$selectedEvents <- FALSE
params$markerSelectedEvents <- NA # have to be given if selectedEvents = TRUE

if(params$selectedEvents == T){
  set.seed(params$seed)
  random_events_data <- random_events_from_selected_events(posNeg = posNeg, marker = params$markerSelectedEvents, n = params$n_per_file)
  params$ext_name <- params$markerSelectedEvents
} else {
  set.seed(params$seed)
  number_of_events_data <-  number_of_events(data = fcs_data, file_names = file_names)
  random_events_data <- random_events(number_of_events_data, n = params$n_per_file)
  params$ext_name <- "All"
}
#

arcSindataMatrix <- list_to_matrix_selected_events(data = fcs_data,
                                                   kept_events = random_events_data,
                                                   file_names = file_names,
                                                   channels = kanalnavn,
                                                   archSin = T,
                                                   cofactor = 5,
                                                   scale = F
)
colnames(arcSindataMatrix)[1:length(kanaler)] <- kanaler

params$data = arcSindataMatrix[arcSindataMatrix$dataset %in% tamed, ]
params$kanaler = orderPanel1 #ta eventuelt ut de du ikke skal bruke
params$scaling = TRUE
params$column_cluster = FALSE #kluster ikke på kolonner.

params$utSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
if(!file.exists(params$utSti)){
  dir.create(params$utSti)
}
source(fs::path(scriptPath_UnsupAnalysis, "FlowSOM_analyse.R"))

write.csv2(tamed, fs::path(params$utSti, "radnavn.csv"))

#seed 1345 panel 1 alle

params <- list()
params$panel <- "Panel 1"
params$seed <- 1345
params$ks <- c(10, 20, 30, 40, 50, 60)
params$n_per_file <- 25000
params$xdim <- 14
params$ydim <- 14
params$selectedEvents <- FALSE
params$markerSelectedEvents <- NA # have to be given if selectedEvents = TRUE

if(params$selectedEvents == T){
  set.seed(params$seed)
  random_events_data <- random_events_from_selected_events(posNeg = posNeg, marker = params$markerSelectedEvents, n = params$n_per_file)
  params$ext_name <- params$markerSelectedEvents
} else {
  set.seed(params$seed)
  number_of_events_data <-  number_of_events(data = fcs_data, file_names = file_names)
  random_events_data <- random_events(number_of_events_data, n = params$n_per_file)
  params$ext_name <- "All"
}

arcSindataMatrix <- list_to_matrix_selected_events(data = fcs_data,
                                                   kept_events = random_events_data,
                                                   file_names = file_names,
                                                   channels = kanalnavn,
                                                   archSin = T,
                                                   cofactor = 5,
                                                   scale = F
)
colnames(arcSindataMatrix)[1:length(kanaler)] <- kanaler

params$data = arcSindataMatrix[arcSindataMatrix$dataset %in% tamed, ]
params$kanaler = orderPanel1 #ta eventuelt ut de du ikke skal bruke
params$scaling = TRUE
params$column_cluster = FALSE #kluster ikke på kolonner.

params$utSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
if(!file.exists(params$utSti)){
  dir.create(params$utSti)
}

source(fs::path(scriptPath_UnsupAnalysis, "FlowSOM_analyse.R"))


write.csv2(tamed, fs::path(params$utSti, "radnavn.csv"))



#seed 1456 panel 1 alle

params <- list()
params$panel <- "Panel 1"
params$seed <- 1456
params$ks <- c(10, 20, 30, 40, 50, 60)
params$n_per_file <- 25000
params$xdim <- 14
params$ydim <- 14
params$selectedEvents <- FALSE
params$markerSelectedEvents <- NA # have to be given if selectedEvents = TRUE

if(params$selectedEvents == T){
  set.seed(params$seed)
  random_events_data <- random_events_from_selected_events(posNeg = posNeg, marker = params$markerSelectedEvents, n = params$n_per_file)
  params$ext_name <- params$markerSelectedEvents
} else {
  set.seed(params$seed)
  number_of_events_data <-  number_of_events(data = fcs_data, file_names = file_names)
  random_events_data <- random_events(number_of_events_data, n = params$n_per_file)
  params$ext_name <- "All"
}

arcSindataMatrix <- list_to_matrix_selected_events(data = fcs_data,
                                                   kept_events = random_events_data,
                                                   file_names = file_names,
                                                   channels = kanalnavn,
                                                   archSin = T,
                                                   cofactor = 5,
                                                   scale = F
)
colnames(arcSindataMatrix)[1:length(kanaler)] <- kanaler

params$data = arcSindataMatrix[arcSindataMatrix$dataset %in% tamed, ]
params$kanaler = orderPanel1 #ta eventuelt ut de du ikke skal bruke
params$scaling = TRUE
params$column_cluster = FALSE #kluster ikke på kolonner.

params$utSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
if(!file.exists(params$utSti)){
  dir.create(params$utSti)
}

source(fs::path(scriptPath_UnsupAnalysis, "FlowSOM_analyse.R"))

write.csv2(tamed, fs::path(params$utSti, "radnavn.csv"))



#seed 1567 panel 1 alle

params <- list()
params$panel <- "Panel 1"
params$seed <- 1567
params$ks <- c(10, 20, 30, 40, 50, 60)
params$n_per_file <- 25000
params$xdim <- 14
params$ydim <- 14
params$selectedEvents <- FALSE
params$markerSelectedEvents <- NA # have to be given if selectedEvents = TRUE

if(params$selectedEvents == T){
  set.seed(params$seed)
  random_events_data <- random_events_from_selected_events(posNeg = posNeg, marker = params$markerSelectedEvents, n = params$n_per_file)
  params$ext_name <- params$markerSelectedEvents
} else {
  set.seed(params$seed)
  number_of_events_data <-  number_of_events(data = fcs_data, file_names = file_names)
  random_events_data <- random_events(number_of_events_data, n = params$n_per_file)
  params$ext_name <- "All"
}

arcSindataMatrix <- list_to_matrix_selected_events(data = fcs_data,
                                                   kept_events = random_events_data,
                                                   file_names = file_names,
                                                   channels = kanalnavn,
                                                   archSin = T,
                                                   cofactor = 5,
                                                   scale = F
)
colnames(arcSindataMatrix)[1:length(kanaler)] <- kanaler

params$data = arcSindataMatrix[arcSindataMatrix$dataset %in% tamed, ]
params$kanaler = orderPanel1 #ta eventuelt ut de du ikke skal bruke
params$scaling = TRUE
params$column_cluster = FALSE #kluster ikke på kolonner.

params$utSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
if(!file.exists(params$utSti)){
  dir.create(params$utSti)
}

source(fs::path(scriptPath_UnsupAnalysis, "FlowSOM_analyse.R"))

write.csv2(tamed, fs::path(params$utSti, "radnavn.csv"))

#panel 2

rm(list = ls())


orderPanel2 <- c("89Y_CD45", "116Cd_CD3", "113Cd_CD8", "145Nd_CD4",  "111Cd_CD19",  "153Eu_CXCR5-CD185",
                 "114Cd_HLADR", "175Lu_CD14", "209Bi_CD16",  "176Yb_CD56", "106Cd_CD57", "158Gd_CD27",
                 "160Gd_CD28", "169Tm_CD25", "172Yb_CD38", "112Cd_CD44", "143Nd_CD127-IL7Ra",  "167Er_CCR7-CD197",
                 "155Gd_CD45RA", "162Dy_FoxP3", "170Er_CTLA-4-CD152", "163Dy_CD33",  "154Sm_CD272-BTLA", "110Cd_CD107a",
                 "161Dy_IL-17A", "166Er_IL-10", "164Dy_Perforin", "171Yb_GranzymeB", "148Nd_CD274-PD-L1",
                 "173Yb_CD273-PD-L2", "159Tb_GM-CSF",   "168Er_CD154-CD40L", "174Yb_CD279-PD1", "141Pr_CD223-LAG3",
                 "147Sm_TIM-3",  "151Eu_CD137-4-1BB",  "165Ho_IFNg", "142Nd_IL-1b", "156Gd_IL-6", "166Er_IL-10",
                 "152Sm_TCRgd", "151Eu_CD137-4-1BB", "149Sm_IL-12p70", "150Nd_MIP-1b", "144Nd_IL-2",  "146Nd_TNFa")



scriptPath_UnsupAnalysis <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "Cytof unsupperviced")

source(fs::path(scriptPath_UnsupAnalysis, "readDataToAnalysepanel2TighterCleanUp.R"))

tamed <- as.character(dInfo$filnavn)
tamed <- tamed[!grepl("FHI005", tamed)] #denne oppfører seg helt rart og tas ut av analysen.
tamed <- tamed[!grepl("Ref", tamed)]


#seed 2234 Panel 2 alle

params <- list()
params$panel <- "Panel 2"
params$seed <- 2234
params$ks <- c(10, 20, 30, 40, 50, 60)
params$n_per_file <- 25000
params$xdim <- 14
params$ydim <- 14
params$selectedEvents <- FALSE
params$markerSelectedEvents <- NA # have to be given if selectedEvents = TRUE

if(params$selectedEvents == T){
  set.seed(params$seed)
  random_events_data <- random_events_from_selected_events(posNeg = posNeg, marker = params$markerSelectedEvents, n = params$n_per_file)
  params$ext_name <- params$markerSelectedEvents
} else {
  set.seed(params$seed)
  number_of_events_data <-  number_of_events(data = fcs_data, file_names = file_names)
  random_events_data <- random_events(number_of_events_data, n = params$n_per_file)
  params$ext_name <- "All"
}

arcSindataMatrix <- list_to_matrix_selected_events(data = fcs_data,
                                                   kept_events = random_events_data,
                                                   file_names = file_names,
                                                   channels = kanalnavn,
                                                   archSin = T,
                                                   cofactor = 5,
                                                   scale = F
)
colnames(arcSindataMatrix)[1:length(kanaler)] <- kanaler

params$data = arcSindataMatrix[arcSindataMatrix$dataset %in% tamed, ]
params$kanaler = orderPanel2 #ta eventuelt ut de du ikke skal bruke
params$scaling = TRUE
params$column_cluster = FALSE #kluster ikke på kolonner.

params$utSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
if(!file.exists(params$utSti)){
  dir.create(params$utSti)
}

source(fs::path(scriptPath_UnsupAnalysis, "FlowSOM_analyse.R"))
write.csv2(tamed, fs::path(params$utSti, "radnavn.csv"))


#seed 2345 Panel 2 alle

params <- list()
params$panel <- "Panel 2"
params$seed <- 2345
params$ks <- c(10, 20, 30, 40, 50, 60)
params$n_per_file <- 25000
params$xdim <- 14
params$ydim <- 14
params$selectedEvents <- FALSE
params$markerSelectedEvents <- NA # have to be given if selectedEvents = TRUE

if(params$selectedEvents == T){
  set.seed(params$seed)
  random_events_data <- random_events_from_selected_events(posNeg = posNeg, marker = params$markerSelectedEvents, n = params$n_per_file)
  params$ext_name <- params$markerSelectedEvents
} else {
  set.seed(params$seed)
  number_of_events_data <-  number_of_events(data = fcs_data, file_names = file_names)
  random_events_data <- random_events(number_of_events_data, n = params$n_per_file)
  params$ext_name <- "All"
}

arcSindataMatrix <- list_to_matrix_selected_events(data = fcs_data,
                                                   kept_events = random_events_data,
                                                   file_names = file_names,
                                                   channels = kanalnavn,
                                                   archSin = T,
                                                   cofactor = 5,
                                                   scale = F
)
colnames(arcSindataMatrix)[1:length(kanaler)] <- kanaler

params$data = arcSindataMatrix[arcSindataMatrix$dataset %in% tamed, ]
params$kanaler = orderPanel2 #ta eventuelt ut de du ikke skal bruke
params$scaling = TRUE
params$column_cluster = FALSE #kluster ikke på kolonner.

params$utSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
if(!file.exists(params$utSti)){
  dir.create(params$utSti)
}

source(fs::path(scriptPath_UnsupAnalysis, "FlowSOM_analyse.R"))
write.csv2(tamed, fs::path(params$utSti, "radnavn.csv"))




#seed 2456 Panel 2 alle

params <- list()
params$panel <- "Panel 2"
params$seed <- 2456
params$ks <- c(10, 20, 30, 40, 50, 60)
params$n_per_file <- 25000
params$xdim <- 14
params$ydim <- 14
params$selectedEvents <- FALSE
params$markerSelectedEvents <- NA # have to be given if selectedEvents = TRUE

if(params$selectedEvents == T){
  set.seed(params$seed)
  random_events_data <- random_events_from_selected_events(posNeg = posNeg, marker = params$markerSelectedEvents, n = params$n_per_file)
  params$ext_name <- params$markerSelectedEvents
} else {
  set.seed(params$seed)
  number_of_events_data <-  number_of_events(data = fcs_data, file_names = file_names)
  random_events_data <- random_events(number_of_events_data, n = params$n_per_file)
  params$ext_name <- "All"
}

arcSindataMatrix <- list_to_matrix_selected_events(data = fcs_data,
                                                   kept_events = random_events_data,
                                                   file_names = file_names,
                                                   channels = kanalnavn,
                                                   archSin = T,
                                                   cofactor = 5,
                                                   scale = F
)
colnames(arcSindataMatrix)[1:length(kanaler)] <- kanaler

params$data = arcSindataMatrix[arcSindataMatrix$dataset %in% tamed, ]
params$kanaler = orderPanel2 #ta eventuelt ut de du ikke skal bruke
params$scaling = TRUE
params$column_cluster = FALSE #kluster ikke på kolonner.

params$utSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
if(!file.exists(params$utSti)){
  dir.create(params$utSti)
}

source(fs::path(scriptPath_UnsupAnalysis, "FlowSOM_analyse.R"))

write.csv2(tamed, fs::path(params$utSti, "radnavn.csv"))



#seed 2567 Panel 2 alle

params <- list()
params$panel <- "Panel 2"
params$seed <- 2567
params$ks <- c(10, 20, 30, 40, 50, 60)
params$n_per_file <- 25000
params$xdim <- 14
params$ydim <- 14
params$selectedEvents <- FALSE
params$markerSelectedEvents <- NA # have to be given if selectedEvents = TRUE

if(params$selectedEvents == T){
  set.seed(params$seed)
  random_events_data <- random_events_from_selected_events(posNeg = posNeg, marker = params$markerSelectedEvents, n = params$n_per_file)
  params$ext_name <- params$markerSelectedEvents
} else {
  set.seed(params$seed)
  number_of_events_data <-  number_of_events(data = fcs_data, file_names = file_names)
  random_events_data <- random_events(number_of_events_data, n = params$n_per_file)
  params$ext_name <- "All"
}

arcSindataMatrix <- list_to_matrix_selected_events(data = fcs_data,
                                                   kept_events = random_events_data,
                                                   file_names = file_names,
                                                   channels = kanalnavn,
                                                   archSin = T,
                                                   cofactor = 5,
                                                   scale = F
)
colnames(arcSindataMatrix)[1:length(kanaler)] <- kanaler

params$data = arcSindataMatrix[arcSindataMatrix$dataset %in% tamed, ]
params$kanaler = orderPanel2 #ta eventuelt ut de du ikke skal bruke
params$scaling = TRUE
params$column_cluster = FALSE #kluster ikke på kolonner.

params$utSti <- fs::path(scriptPath_UnsupAnalysis, params$panel, params$ext_name, paste("seed", params$seed))
if(!file.exists(params$utSti)){
  dir.create(params$utSti)
}

source(fs::path(scriptPath_UnsupAnalysis, "FlowSOM_analyse.R"))
write.csv2(tamed, fs::path(params$utSti, "radnavn.csv"))

