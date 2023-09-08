# les inn stier

data_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel2_mars2022", "clean data")
scriptPath <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse")
posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel2_mars2022", "posNeg", "Data")

# les inn funksjoner 
source(fs::path(scriptPath,  "read_data_functions.R"))
source(fs::path(scriptPath,  "transformation_functions.R"))
source(fs::path(scriptPath,  "ploting_functions.R"))
source(fs::path(scriptPath,  "clustering_functions.R"))
source(fs::path(scriptPath,  "gating_functions.R"))
source(fs::path(scriptPath,  "analysis_functions.R"))



#lager dInfo fra filnavnene
posNegFilnavn <- as.character(readRDS(fs::path(posNeg_path,  "posNegFilnavn.rds")))
posNeg <- readRDS(fs::path(posNeg_path,  "posNeg.rds"))




# read all files in data_path into one dataset fcs_data

setwd(data_path)
fcs_data_with_info <- read_specific_data_from_folder(data_path = data_path, files_to_open = paste0(posNegFilnavn, ".fcs"))
fcs_data <- fcs_data_with_info$fcs_data
file_names <- fcs_data_with_info$file_names
rm(fcs_data_with_info)

params_fcs <- get_params_fcs_data(fcs_data[[1]])

kanaler <- params_fcs$desc[grepl("_",params_fcs$desc)]
kanaler <- kanaler[!grepl("_DNA", kanaler)]
kanaler <- kanaler[!grepl("_Cisp", kanaler)]
kanalnavn <- params_fcs$name[params_fcs$desc %in% kanaler]
print("antall markører med")
length(kanalnavn)


kortFilNavn <- read.csv2(fs::path(posNeg_path, "kortFilNavn.csv"))
rownames(kortFilNavn) <-  kortFilNavn$filnavn

kortFilNavn <- kortFilNavn[posNegFilnavn,]


dInfo <- kortFilNavn

dInfo$statusTid <- dInfo$status
dInfo$status <- gsub(" T1", "", dInfo$status)
dInfo$status <- gsub(" T2", "", dInfo$status)

dInfo$tid <- "Control"
dInfo$tid[grep("T1", dInfo$kort_filnavn)] <- "T1"
dInfo$tid[grep("T2", dInfo$kort_filnavn)] <- "T2"
dInfo$tid[grep("Ref", dInfo$kort_filnavn)] <- "Ref"

# 
# 
# dInfo$tid <- NA
# for(i in 1:nrow(tid)){
#   dInfo$tid[grep(as.character(tid$pasient[i]), dInfo$filenames)] <- as.character(tid$tid[i])
# }
# 
# dInfo$tid[grep(as.character("FHI81"), dInfo$filenames)] <- "28.06.2021"
# dInfo$tid[grep(as.character("FHI95"), dInfo$filenames)] <- "28.06.2021"
# 
