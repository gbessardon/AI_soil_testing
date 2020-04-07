library(compositions)
library(caret)
library(earth)
library(readxl)
library(doParallel)
library(secr)
library(cluster)
library(raster)
library (rgeos)

texture_data <- read_xlsx('G:/TOPSOIL/LUCAS_TOPSOIL_v1.xlsx', sheet = "Sheet1", range = cell_cols("C:E"))
covariatesvar <- read_xlsx('G:/TOPSOIL/LUCAS_TOPSOIL_v1.xlsx', sheet = "Sheet1", range = cell_cols("F:M"))
lat_lon <- read_xlsx('G:/TOPSOIL//LUCAS_TOPSOIL_v1.xlsx', sheet = "Sheet1", range = cell_cols("P:Q"))
# perform additive log ratio tranformation of texture data (sand, silt, clay)
alr1 <- ilr(texture_data)

#fit models on the two ilr
ctrl <- safsControl(functions=caretSA, verbose=T, method="cv", allowParallel = TRUE, number=2)
training <- trainControl( allowParallel=T, method="cv", verbose=T, number=2, search = "random")
registerDoParallel(40)
model_earth_alr1 <- safs(x = texture_data , y =alr1[,1], selectSize = pickSizeBest,iters=50, safsControl = ctrl, trControl= training, metric = "RMSE", maximize=F, method="lm", tuneLength=20)
stopImplicitCluster()
closeAllConnections()

ctrl <- safsControl(functions=caretSA, verbose=T, method="cv", allowParallel = TRUE, number=2)
training <- trainControl( allowParallel=T, method="cv", verbose=T, number=2, search = "random")
registerDoParallel(40)
model_earth_alr2 <- safs(x = texture_data, y =alr1[,2], selectSize = pickSizeBest,iters=50, safsControl = ctrl, trControl= training, metric = "RMSE", maximize=F, method="lm", tuneLength=20)
stopImplicitCluster()
closeAllConnections()



#predict the two ilr rasters
#Load the model raster
fn<-"G:/TOPSOIL/Clay_eu23.tif"
Original<-raster(fn)
d<-dim(Original)
reso<-res(Original)
eo<-extent(Original)
projo<-projection(Original)
spg1<-na.omit(data.frame(lat_lon[,2],lat_lon[,1],alr1[,1]))
colnames(spg1) <- c('x', 'y', 'z')
coordinates(spg1) <- c("x", "y")
projection(spg1) <- CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
plot(spg1)
spg1p<-spTransform(spg1,projo)
dfspg1p<-data.frame(spg1p)
r<-raster(resolution=reso*10,eo,projo)
x<-rasterize(dfspg1p[,2:3],r,dfspg1p[,1],fun=max)
projection(x)<-CRS(projo)
plot(x)

beginCluster()
clusterR(x,sqrt,verbose=T)
endCluster()
beginCluster(20)
clusterR(x, predict, args=list(model=model_earth_alr1, progress="text"), filename="yourfilename_alr1.tif", overwrite=T)
endCluster()

spg2<-na.omit(data.frame(lat_lon,alr1[,2]))
colnames(spg2) <- c('x', 'y', 'z')
e <- extent(spg2[,1:2])
r<-raster(nrows=d[2],ncols=d[1],e,'+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs' )
x<-rasterize(spg2[,1:2],r,spg2[,3],fun=mean)
rp2<-projectRaster(x,Original)
beginCluster(20)
clusterR(rp2, predict, args=list(model=model_earth_alr2, progress="text"), filename="yourfilename_alr2.tif", overwrite=T)
endCluster()


#tranform back the ilr values into sand silt and clay

ilr.inv.ras <- function(x){
  if(anyNA(x)){
    return(rep(NA, 3))} else {
      return(ilrInv(x))}
}

alrs <- stack("yourfilename_alr1.tif", "yourfilename_alr2.tif")

calc(alrs, ilr.inv.ras, filename="yourfilename_texture.tif", overwrite=T)