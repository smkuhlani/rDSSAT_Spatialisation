rm(list=ls()) # Cleaning the shit off
gc()

ptm <- proc.time() # start code time
zita1<-c('M312')
# zita2<-c('M171')
# zita3<-c('M181')

# f = list.dirs(path = "F:/DSSATwth/Uganda/CURR/W5E5", full.names = TRUE, recursive = TRUE)
f = list.dirs(path = "C:/DSSAT47/Weather/Nigeria8139", full.names = TRUE, recursive = TRUE)
f1 = f[2:length(f)]
for (i in f1){
  lf = list.files(i,pattern="*.WTH",full.names = T)
  file.copy(from=lf, to="C:/DSSAT47/Weather/",overwrite = TRUE, recursive = FALSE,copy.mode = TRUE)
  setwd("C:/DSSAT47/Seasonal/")
  # shell("C:/DSSAT47/DSCSM047.EXE CSCER047 A M1610001.SNX")
  shell("C:/DSSAT47/DSCSM047.EXE CSCER047 B NigeriaBatch312.v47")
  #names = substr(i, 20, nchar(i)) #number of letters in code before model name
  names = zita1 #number of letters in code before model name
      print(paste0("Simulating for ",i))
        print(paste0("model simulaion is completed for ",names))
        file.copy("C:/DSSAT47/Seasonal/Summary.OUT",
                  to=paste0("D:/~/rDSSAT_Spatialisation/processing/",names,".OUT"),
                  overwrite = TRUE, recursive = FALSE,copy.mode = TRUE)
}

time <- proc.time() - ptm # end code timer
time/60




