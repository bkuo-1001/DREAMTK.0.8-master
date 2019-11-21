#PRECONDITIONS, you have ran the SQL Script in the same folder and you have the model data in the R memory.

#now we load the files we created with the sql script
#file paths are hard coded oops.
#gotta make sure that "" is North America. (N/A values)
ac50_data = read.csv("C:/ProgramData/MySQL/MySQL Server 5.6/ac50_data.csv",header = TRUE, sep = ",", na.strings = "", stringsAsFactors = FALSE);
assay_data = read.csv("C:/ProgramData/MySQL/MySQL Server 5.6/assay_data.csv",header = TRUE, sep = ";", na.strings = "", stringsAsFactors = FALSE);
chemical_data = read.csv("C:/ProgramData/MySQL/MySQL Server 5.6/chemical_data.csv",header = TRUE, sep = ";", na.strings = "", stringsAsFactors = FALSE);
cpdat_data = read.csv("C:/ProgramData/MySQL/MySQL Server 5.6/cpdat_data.csv",header = TRUE, sep = "	", na.strings = "", stringsAsFactors = FALSE);
shedsinfo_data = read.csv("C:/ProgramData/MySQL/MySQL Server 5.6/shedsinfo_data.csv",header = TRUE, sep = "	", na.strings = "", stringsAsFactors = FALSE);

ac50_data$ac50 <- as.numeric(ac50_data$ac50);
ac50_data$aeid <- as.numeric(ac50_data$aeid);
ac50_data$ac_top <- as.numeric(ac50_data$ac_top);
ac50_data$ac_cutoff <- as.numeric(ac50_data$ac_cutoff);
ac50_data$casn <- as.character(ac50_data$casn);




chemical_data$casn <- as.character(chemical_data$casn);
chemical_data$name <- as.character(chemical_data$name);
chemical_data$cytotoxicity_um <- as.numeric(chemical_data$cytotoxicity_um);
chemical_data$cytotoxic <- as.character(chemical_data$cytotoxic);

#this will save the dataset onto the base R folder where you wil have to paste it where you need it .
save(ac50_data, assay_data, chemical_data, kow_data, model_3css_data, physiology_data, pka_data, shedsinfo_data, cpdat_data, compress = TRUE, compression_level = 9, file = "DreamTKv0.7.RData")


#}