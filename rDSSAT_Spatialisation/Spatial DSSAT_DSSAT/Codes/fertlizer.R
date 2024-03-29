# Sys.setenv(JAVA_HOME="C:/Program Files/Java/jdk1.8.0_181/")
# options(java.parameters = "-Xmx4g")
library(readxl)
library(gdata)
library(tidyverse)
library(miceadds)
library(data.table)
##setwd and path Major FILEX path
# start_y <- "3122"
# y1 <- as.Date(start_y,"%Y")
# y2 <- format(y1,"%y")
setwd("D:/OneDrive - CGIAR/Personal/Personalresearch/Colleagues/Abellac/DSSAT_Spatialisation/PDSSAT/Nigeria/filex/")
FILEX <- "D:/OneDrive - CGIAR/Personal/Personalresearch/Colleagues/Abellac/DSSAT_Spatialisation/PDSSAT/Nigeria/filex/FILEX_DSSAT_NG312.xlsx"

### Read write and treatments in DSSAT Part
sheet <- read_excel(FILEX,sheet = "Fertlizers_Inorganic")
## split sheet intoi multiple dat frames for creating mutple filex files
sheet_1 <- sheet[,c(1,2,4:(ncol(sheet)-3))]
#reshaping the second part for keeping all the depths in same columnn
sheet_2 <- melt(data = sheet_1,id.vars = c("Name","@F"))
## split name by _
n1 <- gsub("[[:digit:]]","",sheet_2$variable)
n2 <- gsub("[[:alpha:]]","",sheet_2$variable)
sheet_2$variable <- paste0(n1,"_",n2)
sheet_3 <- data.frame(do.call('rbind', strsplit(as.character(sheet_2$variable),'_',fixed=TRUE)))
sheet_2$variable <- sheet_3$X1
##REmove NA values rows
sheet_2[sheet_2=="NA"] <- NA
sheet_f <- na.omit(sheet_2)
names(sheet_f) <- c("Name","F" ,"variable","value")
## reshaping wide again to split at every station
sheet_f$samples <- rownames(sheet_f)
sheet_4 <- spread(sheet_f,variable,value)
#sheet_4 <- reshape(sheet_f,timevar = "variable",idvar = c("Name","F","samples"),direction = "wide")

### removing na values from one by one column 
FDATE <- na.omit(sheet_4[,c(1,2,10)])
FMCD <- na.omit(sheet_4[,c(1,2,13)])
FACD <- na.omit(sheet_4[,c(1,2,4)])
FDEP <- na.omit(sheet_4[,c(1,2,11)])
FAMN <- na.omit(sheet_4[,c(1,2,7)])
FAMP <- na.omit(sheet_4[,c(1,2,9)])
FAMK <- na.omit(sheet_4[,c(1,2,6)])
FAMC <- na.omit(sheet_4[,c(1,2,5)])
FAMO <- na.omit(sheet_4[,c(1,2,8)])
FOCD <- na.omit(sheet_4[,c(1,2,14)])
FERNAME <- na.omit(sheet_4[,c(1,2,12)])

##cbind all columns
sheet_5 <- cbind(FDATE[,1:3],FMCD[,3],FACD[,3],FDEP[,3],FAMN[,3],FAMP[,3],FAMK[,3],FAMC[,3],FAMO[,3],FOCD[,3],FERNAME[,3])
#sheet_5$FDATE <- paste0(sprintf("%03d",as.numeric(sheet_5$FDATE)))
names(sheet_5) <- c("Name","@F","FDATE","FMCD","FACD","FDEP","FAMN","FAMP","FAMK","FAMC","FAMO","FOCD","FERNAME")

## split sheet intoi multiple dat frames for creating mutple filex files
dfs <- split(sheet_5[2:ncol(sheet_5)],sheet_5$Name)

##creat directory for saving multple filex treat ments
dir.create(file.path(getwd(),"FR1"),showWarnings = FALSE)
dir.create(file.path(getwd(),"FR"),showWarnings = FALSE)
# splitting all filex into seperate csv files
for (i in names(dfs)){
  write.csv(dfs[[i]],paste0(getwd(),"/FR1/",i,".csv"),row.names = F)
}
# read all csv files and write in fixed widh text files and remove csv files 
FILES <- list.files(path = paste0(getwd(),"/FR1/"),pattern=".csv$",full.names = T)
for (i in 1:length(FILES)) {
  FILE=read.csv(file=FILES[i])
  colnames(FILE)[1] <- "F"
  FILE$F <- sprintf("%2d",as.integer(FILE$F))
  FILE$FDATE <- sprintf("%5d",as.integer(FILE$FDATE))
  FILE$FMCD <- sprintf("%5s",FILE$FMCD)
  FILE$FACD <- sprintf("%5s",FILE$FACD)
  FILE$FDEP <- sprintf("%5d",as.integer(FILE$FDEP))
  FILE$FAMN <- sprintf("%5d",as.integer(FILE$FAMN))
  FILE$FAMP <- sprintf("%5d",as.integer(FILE$FAMP))
  FILE$FAMK <- sprintf("%5d",as.integer(FILE$FAMK))
  FILE$FAMC <- sprintf("%5d",as.integer(FILE$FAMC))
  FILE$FAMO <- sprintf("%5d",as.integer(FILE$FAMO))
  FILE$FOCD <- sprintf("%5d",as.integer(FILE$FOCD))
  FILE$FERNAME <- sprintf("%-12s",FILE$FERNAME)
  colnames(FILE)[1] <- "@F"
  hdr <- "*FERTILIZERS (INORGANIC)"
  write.fwf(as.data.frame(hdr),"hdr.txt",quoteInfo = F,rownames = F,
            colnames = F)
  #w1 <- do.call(pmax, as.data.frame(nchar(as.matrix(FILE[,1:ncol(FILE)]))))
  write.table(FILE,file=paste0(sub(".csv","",FILES[i]),".txt"),row.names = F,col.names = T,quote = F)
}
#file.remove(FILES)

### Pate the header and dataframes 
genmodeldats <- function(inputfold,header,outdir){
  # Find common filenames in each folder
  filex <- tools::file_path_sans_ext(list.files(inputfold,pattern = "^.*\\.txt$"))
  # switchf <- tools::file_path_sans_ext(list.files(switchfold,pattern = "^.*\\.DAT$"))
  commonf <- Reduce(base::intersect,list(filex))
  # Now append all these files and write to outdir
  for(f in commonf){
    # New file name
    modeldat <- file.path(outdir,paste0(f,".txt"))
    #path for spliited filex
    filex1 <- file.path(inputfold,paste0(f,".txt"))
    cat("",file=modeldat,append=FALSE)
    for(i in c(header,filex1)){
      file.append(modeldat,i)
    }
    message("Filex written written to ",modeldat)
  }
}
genmodeldats(paste0(getwd(),"/FR1/"),header = "hdr.txt",paste0(getwd(),"/FR"))
## removing the old files 
file.remove("hdr.txt")
unlink("FR1",recursive = T)




####






















# # Sys.setenv(JAVA_HOME="C:/Program Files/Java/jdk1.8.0_181/")
# # options(java.parameters = "-Xmx4g")
# library(readxl)
# library(gdata)
# library(tidyverse)
# library(miceadds)
# library(data.table)
# ##setwd and path Major FILEX path
# # start_y <- "3122"
# # y1 <- as.Date(start_y,"%Y")
# # y2 <- format(y1,"%y")
# setwd("D:/OneDrive - CGIAR/Personal/Personalresearch/Colleagues/Abellac/DSSAT_Spatialisation/PDSSAT/Nigeria/filex/")
# FILEX <- "D:/OneDrive - CGIAR/Personal/Personalresearch/Colleagues/Abellac/DSSAT_Spatialisation/PDSSAT/Nigeria/filex/FILEX_DSSAT_NG312.xlsx"
# 
# ### Read write and treatments in DSSAT Part
# sheet <- read_excel(FILEX,sheet = "Fertlizers_Inorganic")
# ## split sheet intoi multiple dat frames for creating mutple filex files
# sheet_1 <- sheet[,c(1,2,4:(ncol(sheet)-3))]
# #reshaping the second part for keeping all the depths in same columnn
# sheet_2 <- melt(data = sheet_1,id.vars = c("Name","@F"))
# ## split name by _
# n1 <- gsub("[[:digit:]]","",sheet_2$variable)
# n2 <- gsub("[[:alpha:]]","",sheet_2$variable)
# sheet_2$variable <- paste0(n1,"_",n2)
# sheet_3 <- data.frame(do.call('rbind', strsplit(as.character(sheet_2$variable),'_',fixed=TRUE)))
# sheet_2$variable <- sheet_3$X1
# ##REmove NA values rows
# sheet_2[sheet_2=="NA"] <- NA
# sheet_f <- na.omit(sheet_2)
# names(sheet_f) <- c("Name","F" ,"variable","value")
# ## reshaping wide again to split at every station
# sheet_f$samples <- rownames(sheet_f)
# sheet_4 <- spread(sheet_f,variable,value)
# #sheet_4 <- reshape(sheet_f,timevar = "variable",idvar = c("Name","F","samples"),direction = "wide")
# 
# ### removing na values from one by one column 
# FDATE <- na.omit(sheet_4[,c(1,2,10)])
# FMCD <- na.omit(sheet_4[,c(1,2,13)])
# FACD <- na.omit(sheet_4[,c(1,2,4)])
# FDEP <- na.omit(sheet_4[,c(1,2,11)])
# FAMN <- na.omit(sheet_4[,c(1,2,7)])
# FAMP <- na.omit(sheet_4[,c(1,2,9)])
# FAMK <- na.omit(sheet_4[,c(1,2,6)])
# FAMC <- na.omit(sheet_4[,c(1,2,5)])
# FAMO <- na.omit(sheet_4[,c(1,2,8)])
# FOCD <- na.omit(sheet_4[,c(1,2,14)])
# FERNAME <- na.omit(sheet_4[,c(1,2,12)])
# 
# ##cbind all columns
# sheet_5 <- cbind(FDATE[,1:3],FMCD[,3],FACD[,3],FDEP[,3],FAMN[,3],FAMP[,3],FAMK[,3],FAMC[,3],FAMO[,3],FOCD[,3],FERNAME[,3])
# #sheet_5$FDATE <- paste0(sprintf("%03d",as.numeric(sheet_5$FDATE)))
# names(sheet_5) <- c("Name","@F","FDATE","FMCD","FACD","FDEP","FAMN","FAMP","FAMK","FAMC","FAMO","FOCD","FERNAME")
# 
# ## split sheet intoi multiple dat frames for creating mutple filex files
# dfs <- split(sheet_5[2:ncol(sheet_5)],sheet_5$Name)
# 
# ##creat directory for saving multple filex treat ments
# dir.create(file.path(getwd(),"FR1"),showWarnings = FALSE)
# dir.create(file.path(getwd(),"FR"),showWarnings = FALSE)
# # splitting all filex into seperate csv files
# for (i in names(dfs)){
#   write.csv(dfs[[i]],paste0(getwd(),"/FR1/",i,".csv"),row.names = F)
# }
# # read all csv files and write in fixed widh text files and remove csv files 
# FILES <- list.files(path = paste0(getwd(),"/FR1/"),pattern=".csv$",full.names = T)
# for (i in 1:length(FILES)) {
#   FILE=read.csv(file=FILES[i])
#   colnames(FILE)[1] <- "F"
#   FILE$F <- sprintf("%2d",as.integer(FILE$F))
#   FILE$FDATE <- sprintf("%5d",as.integer(FILE$FDATE))
#   FILE$FMCD <- sprintf("%5s",FILE$FMCD)
#   FILE$FACD <- sprintf("%5s",FILE$FACD)
#   FILE$FDEP <- sprintf("%5d",as.integer(FILE$FDEP))
#   FILE$FAMN <- sprintf("%5d",as.integer(FILE$FAMN))
#   FILE$FAMP <- sprintf("%5d",as.integer(FILE$FAMP))
#   FILE$FAMK <- sprintf("%5d",as.integer(FILE$FAMK))
#   FILE$FAMC <- sprintf("%5d",as.integer(FILE$FAMC))
#   FILE$FAMO <- sprintf("%5d",as.integer(FILE$FAMO))
#   FILE$FOCD <- sprintf("%5d",as.integer(FILE$FOCD))
#   FILE$FERNAME <- sprintf("%-12s",FILE$FERNAME)
#   colnames(FILE)[1] <- "@F"
#   hdr <- "*FERTILIZERS (INORGANIC)"
#   write.fwf(as.data.frame(hdr),"hdr.txt",quoteInfo = F,rownames = F,
#             colnames = F)
#   #w1 <- do.call(pmax, as.data.frame(nchar(as.matrix(FILE[,1:ncol(FILE)]))))
#   write.table(FILE,file=paste0(sub(".csv","",FILES[i]),".txt"),row.names = F,col.names = T,quote = F)
# }
# #file.remove(FILES)
# 
# ### Pate the header and dataframes 
# genmodeldats <- function(inputfold,header,outdir){
#   # Find common filenames in each folder
#   filex <- tools::file_path_sans_ext(list.files(inputfold,pattern = "^.*\\.txt$"))
#   # switchf <- tools::file_path_sans_ext(list.files(switchfold,pattern = "^.*\\.DAT$"))
#   commonf <- Reduce(base::intersect,list(filex))
#   # Now append all these files and write to outdir
#   for(f in commonf){
#     # New file name
#     modeldat <- file.path(outdir,paste0(f,".txt"))
#     #path for spliited filex
#     filex1 <- file.path(inputfold,paste0(f,".txt"))
#     cat("",file=modeldat,append=FALSE)
#     for(i in c(header,filex1)){
#       file.append(modeldat,i)
#     }
#     message("Filex written written to ",modeldat)
#   }
# }
# genmodeldats(paste0(getwd(),"/FR1/"),header = "hdr.txt",paste0(getwd(),"/FR"))
# ## removing the old files 
# file.remove("hdr.txt")
# unlink("FR1",recursive = T)
# 
# 
# 
# 
# ####