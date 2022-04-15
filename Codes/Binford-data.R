# Binford data


library(binford)

head(LRB) # 339 observations
LRB=as.data.frame(LRB)
LRBkey1=as.data.frame(LRBkey)

LRB$wc.area
# rm SOAM; ARC;

data=LRB[-which(LRB$wc.area=="SOAM"),]
data=data[-which(data$wc.area=="ARC"),]

LRB$wlocation
# rm Mexico and subartic (mostly canada)  
data=data[-which(data$wlocation=="Mexico"),]
data=data[-which(data$wlocation=="Alaska"),]
data=data[-which(data$wlocation=="Alberta"),]
data=data[-which(data$wlocation=="Northwest Territories"),]
data=data[-which(data$wlocation=="Quebec"),]
data=data[-which(data$wlocation=="Yukon"),]
data=data[-which(data$wlocation=="British Columbia"),]
data=data[-which(data$wlocation=="Saskatchewan"),]

data$wc.area
data$wlocation

# 233 HG people observations

# remove columns that will not be used
data=data[,-c(7:11,19:74,91:110,113:116,121:259,266:506)]
head(data)

var=names(data)
chave=LRBkey1[which(LRBkey1$X==var),]
head(chave)

popdata=data[,c(2:5,9,11:16)]
head(popdata)

setwd("/Users/christiellyborges/Library/Mobile Documents/com~apple~CloudDocs/Em-busca-do-PhD/Languages/Data/Data-ready")
write.table(popdata,"Binford-population.txt",row.names = F)

