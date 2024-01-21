library(DSSAT)
# setwd("D:/OneDrive - CGIAR/Personal/Personalresearch/Colleagues/Abellac/DSSAT_Spatialisation/PDSSAT/Nigeria/")
setwd('D:/OneDrive - CGIAR/Personal/Personalresearch/Colleagues/Abellac/DSSAT_Spatialisation/PDSSAT/Nigeria/filex/')
#FILEX <- "~/DATABANK/PDSSAT/Nigeria/filex/FILEX_DSSAT_NG312.xlsx"
# batch_file_path <- "All.v47"
batch_file_path1 <- "NigeriaBatch312.v47"

f1 <- list.files("FILEX/",pattern = ".SNX$")
# f2 <- list.files("ROI_files/",pattern = ".SNX$")
# path <- "D:/OneDrive - CGIAR/Personal/Personalresearch/Colleagues/Abellac/DSSAT_Spatialisation/PDSSAT/Nigeria/filex/FILEX"

path <- "C:\\DSSAT47\\Seasonal\\"
file <- paste0(path,"\\",f1)

batch_tbl <- data.frame()
for (i in file){
  df <-  data.frame(FILEX=i,TRTNO=1:99,RP=1,SQ=0,OP=0,CO=0)
  batch_tbl <- rbind(batch_tbl,df)
}

write_dssbatch(batch_tbl, file_name = batch_file_path1)

