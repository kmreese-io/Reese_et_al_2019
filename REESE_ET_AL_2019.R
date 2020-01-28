########################################################################################
## DYNAMIC COMMUNITIES ON THE MESA VERDE CUESTA
## KELSEY M. REESE, DONNA M. GLOWACKI, TIMOTHY A. KOHLER
## AMERICAN ANTIQUITY
## 84(4): 728--747
## WINTER 2019, DOI:10.1017/aaq.2019.74
########################################################################################
## AUTHOR: KELSEY M. REESE
########################################################################################

########################################################
# INSTALL REQUIRED PACKAGES AND CREATE PROJECT WORKSPACE
########################################################

r  <-  getOption('repos')
r['CRAN']  <-  'http://cran.us.r-project.org'
options(repos = r)

# Install and load all necessary packages.
install.packages(c('sp','Hmisc','rgeos','rgdal','raster','FedData','scales','splancs','graphics','apcluster','gdistance','grDevices','RColorBrewer','smoother','wordcloud'))
library('sp');library('Hmisc');library('rgeos');library('rgdal');library('raster');library('FedData');library('scales');library('splancs');library('graphics');library('apcluster');library('gdistance');library('grDevices');library('RColorBrewer');library('smoother');library('wordcloud')

# Create and set directory for the project.
dir.create('Users/USERNAME/Documents/PATH/TO/PROJECT/',recursive=T,showWarnings=F)
set.wd('Users/USERNAME/Documents/PATH/TO/PROJECT/')

# Create folders to store results of cost-distance and cluster analyses
dir.create('./ANALYSIS/SITES_cost-distance/',recursive=T,showWarnings=F)
dir.create('./ANALYSIS/SITES_phases/',recursive=T,showWarnings=F)
dir.create('./ANALYSIS/CLUSTERS_phases/',recursive=T,showWarnings=F)
dir.create('./ANALYSIS/CLUSTERS_apcluster/',recursive=T,showWarnings=F)
dir.create('./ANALYSIS/CLUSTERS_cost-distance/',recursive=T,showWarnings=F)
dir.create('./ANALYSIS/RESULTS/',recursive=T,showWarnings=F)

########################################################
# FUNCTIONS REQUIRED FOR ANALYSES
########################################################

# Cost-distance analysis, outputs cost-kilometers
costIterations <- function(x,...) {
  ph.sample <- coordinates(all.households)[sample(nrow(coordinates(all.households)),momentary.households),]  
  cost.km <- (costDistance(mv.cost,ph.sample,ph.sample)) * (5.036742/3.022045)
  diag(cost.km) <- NA
  write.table(cost.km,file=as.character(phase),append=T,row.names=F,col.names=F,sep=',')
  gc()
}

# Cluster analysis
communityClusters <- function(x,...) {
  ph.sample <- all.households[sample(nrow(all.households),momentary.households,replace=F),]
  write.table(ph.sample,file=as.character(sample.sites),append=T,sep=',',col.names=F)
  ph <- matrix(coordinates(ph.sample),ncol=2,byrow=F)
  mv.extent.sample <- sp::spsample(mv.dem.masked.coords,1000,type='random')
  mv.extent.sample <- as.matrix(coordinates(mv.extent.sample))
  getCellArray <- function(x,...) {
    cost.km <- (costDistance(mv.cost,x,mv.extent.sample)) * (5.036742/3.022045)
    return(cost.km)
  }
  cost.km <- apply(ph,1,FUN=function(x,...){ getCellArray(x,...) })
  cost.df <- as.data.frame(cost.km)
  cost.df[cost.df <= null.difference] <- 1
  cost.df[cost.df > null.difference] <- 0
  m.similarity <- negDistMat(t(cost.df),r=2,method='binary')
  n.clusters <- apcluster(m.similarity)
  df.clusters <- as.data.frame.APResult(n.clusters)
  write.table(df.clusters[order(df.clusters[,3],decreasing=F),],file=as.character(phase.clusters),append=T,sep=',',col.names=F)
  gc()
}

# Converts APResult to data frame, from flycircuit/R/clustering.R available on GitHub
as.data.frame.APResult <- function(x,row.names=NULL,optional=F,clusters,...) {
  if(missing(clusters))
    exemplars=names(x@exemplars)
  else
    exemplars=intersect(clusters,names(x@exemplars))
  clusterids=which(names(x@exemplars)%in%exemplars)
  clusters=x@clusters[clusterids]
  cls=sapply(clusters,length)
  ulc=unlist(clusters)
  df=data.frame(exemplar=factor(rep(exemplars,cls)),cluster=rep(clusterids,cls),idx=ulc,item=names(ulc),row.names=names(ulc),optional=optional,stringsAsFactors=FALSE,...=...)
  df
}

# Enclosing created clusters of households
communityEnclosures <- function(x,y) {
  i <- chull(x,y)
  return(areapl(cbind(x[i],y[i])))
}

########################################################
# SET PROJECTION SYSTEM, IMPORT STUDY AREA BOUNDARY, DOWNLOAD DEM, CREATE COST-RASTER, IMPORT SITE DATA
########################################################

# Define the projection systems to be used. 'master.projection' should be set to your projection system that uses UTMs, and 'longlat.projection' should be kept as is. 
master.projection <- sp::CRS('+proj=YOUR +datum=PROJ4 +zone=HERE')
longlat.projection <- sp::CRS('+proj=longlat +datum=WGS84 +ellps=WGS84')

# A shapefile of your study area. Place shapefile in the working directory before running.
study.area <- rgdal::readOGR('./',layer='study_area_shapefile')
projection(study.area) <- master.projection
longlat.study.area <- sp::spTransform(study.area,longlat.projection)

# Federal DEM download. Internet connection is required to run this chunk for the first time.
longlat.dem <- FedData::get_ned(longlat.study.area,label='study.area.dem',res=13,force.redo=F)
dem <- raster::projectRaster(longlat.dem,crs=master.projection)
study.area.dem <- raster::crop(dem,study.area)

# Create the cost-raster via 'gDistance' package. Number of directions can be set to 4, 8, or 16. The higher the number, the more computationally intensive.
n.directions <- 4
heightDiff <- function(x){x[2] - x[1]}
hd <- transition(study.area.dem,heightDiff,directions=n.directions,symm=F)
slope <- geoCorrection(hd)
adj <- adjacent(study.area.dem,cells=1:ncell(x),pairs=T,directions=n.directions)
speed <- slope
speed[adj] <- 6 * 1000 * exp(-3.5 * abs(slope[adj] + 0.05))
cost.raster <- geoCorrection(speed)

# Import your site information with locations.
site.information <- utils::read.csv('./IMPORT/YOUR/site_locations.csv')

# Change the 'UTMEast' and 'UTMNorth' to match the names of your columns with x- and y-coordinates.
COLUMN_EASTING <- base::grep('^UTMEast$',colnames(site.information))
COLUMN_NORTHING <- base::grep('^UTMNorth$',colnames(site.information))

site.coordinates <- base::matrix(NA,nrow=nrow(site.information),ncol=2)
site.coordinates[,1] <- site.coordinates[,COLUMN_EASTING]
site.coordinates[,2] <- site.coordinates[,COLUMN_NORTHING]
site.coordinates <- sp::SpatialPointsDataFrame(coords=site.coordinates,site.information,proj4string=master.projection)

# Clip site locations to your study area boundary, if necessary
site.coordinates <- site.coordinates[study.area,]

########################################################
# CREATE DATASET FOR NULL MODEL AND RUN NULL MODEL.
########################################################

study.area.coordinates <- raster::rasterToPoints(study.area.dem,spatial=T)
null.coordinates <- study.area.coordinates[study.area,]

phase <- base::as.character('./ANALYSIS/SITES_cost-distance/NULL.csv')
all.households <- sp::spsample(null.coordinates,1000,type='random')
momentary.households <- nrow(coordinates(all.households))
costIterations()

# Once the null values are calculated above, the results can be loaded below rather than re-running the null model.
# The 'null' MUST be loaded into the workspace before running the following analyses.
null <- base::as.character('./ANALYSIS/SITES_cost-distance/NULL.csv')

########################################################
# RUN THE COST DISTANCE AND CLUSTER ANALYSES.
########################################################

# Only one time period is shown here, but you can repeat this analysis for each time period in your study, just change the names of each file so they are saved separately.

# Name the output *.csv files by the time period you are calculating.
phase <- as.character('./ANALYSIS/SITES_cost-distance/PERIOD_15.csv')
sample.sites <- as.character('./ANALYSIS/SITES_phases/PERIOD_15.csv')
phase.clusters <- as.character('./ANALYSIS/CLUSTERS_phases/PERIOD_15.csv')
cluster.information <- as.character('./ANALYSIS/CLUSTERS_apcluster/PERIOD_15.csv')
cluster.cost <- as.character('./ANALYSIS/CLUSTERS_cost-distance/PERIOD_15.csv')

# Create an object with only sites occupied during the desired time period.
phase.households <- site.coordinates[which(site.coordinates$PERIOD_1 >= 1),]
all.households <- phase.households[rep(seq_len(dim(phase.households)[1]),phase.households$PERIOD_1),]

# Momentize the number of households, if necessary, by setting the estimated use-life of households and the total length of the study period for the corresponding object below.
# If momentizing households is not necessary, leave each object = 1
use.life.households <- 1
total.length.period <- 1
momentary.households <- round((use.life.households / total.length.period) * nrow(all.households))

# If you have momentized the households and would like to repeat the analysis multiple times, change the 'repeat.analysis' object to your desired number of iterations.
# If you do not want to run multiple iterations for your time period, leave the object = 1
repeat.analysis <- 1
replicate(repeat.analysis,costIterations())
d.null <- density(as.matrix(read.table(null,header=F,sep=',')),na.rm=T,from=0,to=momentary.households,n=512^2)
d.phase <- density(as.matrix(read.table(phase,header=F,sep=',')),na.rm=T,from=0,to=momentary.households,n=512^2)
d.bind <- cbind(as.matrix(d.null$x),as.matrix(d.phase$y-d.null$y))
null.difference <- mean(c(d.bind[which(d.bind[,2] < 0)[1] - 1],d.bind[which(d.bind[,2] < 0)[1] ]))
replicate(repeat.analysis,communityClusters())
write.csv(cbind(read.csv(phase.clusters,header=F),read.csv(sample.sites,header=F)),file=cluster.information)

########################################################
# CALCULATE THE SUMMARY STATISTICS FOR YOUR ANALYSIS
########################################################

## Phase information for output. Enter the time period number, beginning, middle, and end of the time period that you analyzed above.
modeling.phase <- 15
startpoint <- 1100
midpoint <- 1120
endpoint <- 1140

## Input final results of cluster analysis output for results calculations.
total.information <- read.csv(cluster.information,header=T,row.names=1)
household.count <- rep(1,nrow(total.information))
total.information <- cbind(total.information,household.count)
split.iterations <- split(total.information,rep(1:repeat.analysis,each=momentary.households))
split.clusters <- lapply(split.iterations,function(x) split(x,x[,3]))

## Number of clusters
n.clusters <- lapply(split.iterations,function(x,...) max(x[,3]))
clusters.matrix <- matrix(unlist(n.clusters,use.names=F))
mean.clusters <- round(mean(clusters.matrix))
median.clusters <- round(median(clusters.matrix))
sd.clusters <- round(sd(clusters.matrix))

## Households per cluster
n.households.list <- lapply(split.clusters,function(x) lapply(lapply(x,'[[',ncol(total.information)),sum) )
n.households.matrix <- matrix(unlist(n.households.list,use.names=F))
mean.households <- round(mean(n.households.matrix))
median.households <- round(median(n.households.matrix))
sd.households <- round(sd(n.households.matrix))

## Square-kilometers per cluster
km.cluster.list <- lapply(split.clusters,function(x) lapply(x,function(x) communityEnclosures(x[,COLUMN_EASTING+7],x[,COLUMN_NORTHING+7]) / 1000000 ))
km.cluster.matrix <- matrix(unlist(km.cluster.list,use.names=F))
mean.km.cluster <- mean(km.cluster.matrix)
median.km.cluster <- median(km.cluster.matrix)
sd.km.cluster <- sd(km.cluster.matrix)

## Square-kilometers per household
km.per.household <- km.cluster.matrix / n.households.matrix
mean.km.household <- mean(km.per.household)
median.km.household <- median(km.per.household)
sd.km.per.household <- sd(km.per.household)

## Mean cost-distance per cluster
coordinates.list <- lapply(split.clusters,function(x) lapply(x,function(x) cbind(x[,COLUMN_EASTING+7],x[,COLUMN_NORTHING+7]) ))
cost.distance.list <- rapply(coordinates.list,function(x) costDistance(mv.cost,x,x)*(5.036742/3.022045) )
write.csv(as.matrix(cost.distance.list),as.character(cluster.cost))
cost.distance.matrix <- matrix(unlist(cost.distance.list,use.names=F))
cost.distance.ordered <- as.matrix(cost.distance.matrix[order(cost.distance.matrix[,1],decreasing=F),],ncol=1)  
cost.distance.rm <- cost.distance.ordered[-(1:nrow(total.information)),]
mean.distance <- mean(cost.distance.rm)
median.distance <- median(cost.distance.rm)
sd.distance <- sd(cost.distance.rm)

## Mean maximum cost-distance per cluster
cluster.costs <- read.csv(as.character(cluster.cost))
n.digits <- nchar(sub('^0+','',sub('\\.','',cluster.costs[,1])))
truncated <- as.numeric(substr(cluster.costs[,1],1,n.digits-1))
combined <- cbind(truncated,cluster.costs[,2])
splitting <- split(combined[,2],combined[,1])
cluster.maximums <- lapply(splitting,function(x) max(x))
cluster.maximums.unlist <- matrix(unlist(cluster.maximums,use.names=F))
max.mean.distance <- mean(cluster.maximums.unlist)
median.max.distance <- median(cluster.maximums.unlist)
sd.max.distance <- sd(cluster.maximums.unlist)

results <- cbind(modeling.phase,startpoint,midpoint,endpoint,momentary.households,null.difference,mean.clusters,median.clusters,sd.clusters,mean.households,median.households,sd.households,mean.km.cluster,median.km.cluster,sd.km.cluster,mean.km.household,median.km.household,sd.km.per.household,mean.distance,median.distance,sd.distance,max.mean.distance,median.max.distance,sd.max.distance)
write.table(results,file='./ANALYSIS/RESULTS/RESULTS.csv',append=T,row.names=F,col.names=!file.exists('./BASE/RESULTS/RESULTS.csv'),sep=',')
