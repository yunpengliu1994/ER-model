##############################################################
####raster file and worldclim data ################################
##############################################################
## raster of ENA at 10arcmin res ##
library(sf);library(dplyr);library(raster)
rst.10m <- raster("data analysis/wc2.1_10m_bio/wc2.1_10m_elev.tif")
ENA.name=c("ME","NH","VT","NY","MA","RI","CT","NJ","PA","DE","MD","MI","OH","IN","IL","WI","WV","VA","NC","TN","KY","SC","GA","AL","MS","FL")
ENA <- read_sf('data analysis/cb_2018_us_state_20m/cb_2018_us_state_20m.shp')%>%filter(STUSPS%in%ENA.name)%>%as(., "Spatial")
biolist <- list.files("data analysis/wc2.1_10m_bio/", pattern = "*.tif", full.names = TRUE)
biolist <- gtools::mixedsort(sort(biolist))
biostack <- raster::stack(biolist)
#aggregate from 10arcmin resolution to 15arcmin (factor = 1.5)
ENA.rst=ENA%>%raster::mask(biostack,.)%>%raster::crop(extent(ENA))#%>%aggregate(fact=1.5)

domain.ena=raster::as.data.frame(ENA.rst,xy=TRUE)#
colnames(domain.ena)=c("x","y",paste("bio",c(1:19),sep=""),"ele")
domain.ena$Cell.id=seq(1:nrow(domain.ena))

## generate a cell adjacency list that is used for identifying neighbouring cells when simulating spread
# find the row an column of each cell
ncol.grid=length(unique(domain.ena$x))
domain.ena$row <- ((domain.ena$Cell.id - 1) %/% ncol.grid) + 1
domain.ena$col <- ((domain.ena$Cell.id -1 ) %% ncol.grid) + 1

# get the neighbour cells 
domain.ena$below.cell <- domain.ena$Cell.id - ncol.grid
domain.ena$above.cell <- domain.ena$Cell.id + ncol.grid
domain.ena$left.cell <- ifelse(domain.ena$col == 1, domain.ena$Cell.id + (ncol.grid-1) , domain.ena$Cell.id - 1)
domain.ena$right.cell <- ifelse(domain.ena$col == ncol.grid, domain.ena$Cell.id - (ncol.grid-1) , domain.ena$Cell.id + 1)
#remove ocean cells
domain.ena=na.omit(domain.ena)

# add domain ID - numbering used in the range.cells list to be used as indexing
domain.ena$domain.ID <- seq(along=domain.ena$Cell.id)    
# find adjacent cells
# get a rook cell adjacency structure from the domain 
# turn the neighbours from the domain into a matrix  and then a list of cells
domain.NB.mat <- as.matrix(subset(domain.ena, select=c(below.cell, above.cell, left.cell , right.cell)))
dimnames(domain.NB.mat) <- NULL
adjCells.ena <- lapply(seq(along=domain.ena$Cell.id), function(x) domain.NB.mat[x,])
names(adjCells.ena) <- domain.ena$Cell.ID
# now filter out cells which aren't don't appear in the vector of Cell.id
adjCells.ena <- lapply(adjCells.ena, function(X) X[X %in% domain.ena$Cell.id])
adjCells.ena <- lapply(adjCells.ena, function(X) with(domain.ena, domain.ID[match(X, Cell.id)]))
names(adjCells.ena) <- domain.ena$domain.ID

save(domain.ena,adjCells.ena,file="data analysis/Domain.ena.rda")

# make graph object to identify connected sets of cells 
library(graph)
domainGraph.ena <-  new("graphNEL", nodes=as.character(domain.ena$domain.ID), edgeL=adjCells.ena,edgemode="directed")
g1cc <- connComp(domainGraph.ena)
save(g1cc,file="data analysis/connectedGraph.ena.rda")

##############################################################
####get US plat distributions ################################
##############################################################
## step 1.1 --- data downloading of native plants
library(dplyr)
dis <- unique(rWCVPdata::wcvp_distributions$region)
region.us=dis[grep("U.S.A.",dis)]
splist=rWCVPdata::wcvp_distributions
splist.us=subset(splist,splist$region%in%region.us)
spnames <- rWCVPdata::wcvp_names%>%
	filter(accepted_plant_name_id%in%unique(splist.us$plant_name_id)&taxon_status%in%"Accepted"&taxon_rank%in%c("Species","Subspecies","Variety","Form","Subform"))%>%
	select(taxon_name)%>%na.omit()
write.table(spnames,"spnames.us.txt")

# dis <- rWCVPdata::wcvp_distributions[,c("continent","region")]%>%distinct()%>%filter(continent%in%c("ASIA-TROPICAL","ASIA-TEMPERATE","EUROPE","AFRICA")&(!region%in%c("Southern Africa","Macaronesia")))
# splist=rWCVPdata::wcvp_distributions
# splist.native=subset(splist,splist$region%in%dis$region)
# spnames <- rWCVPdata::wcvp_names%>%
	# filter(accepted_plant_name_id%in%unique(splist.native$plant_name_id)&taxon_status%in%"Accepted"&taxon_rank%in%c("Species","Subspecies","Variety","Form","Subform"))%>%
	# select(taxon_name)%>%na.omit()
# write.table(spnames,"spnames.native.txt")

spnames=read.table("data analysis/occ_data.plants/spnames.us.txt")

seql=c(seq(1,nrow(spnames),100),26355)
seql2=lapply(2:length(seql),function(i,seql){if (i<length(seql)){c(seql[i-1],(seql[i]-1))
	} else {
	c(seql[i-1],seql[i])
	}},seql)

# al=do.call(c,lapply(seql2,function(i)i[[2]]))
# done=list.files("data analysis/occ_data.plants", pattern = "*.Rdata", full.names = TRUE)%>%gsub("data analysis/occ_data.plants/","",.)%>%gsub(".Rdata","",.)%>%gsub("occ_data.plants","",.)
# remain=al[!al%in%done]
# seql3=lapply(seql2,function(i,done) if(i[[2]]%in%done) return(i),done)
# seql3=seql3[-which(sapply(seql3, is.null))]
mf=function(i,spnames){
	#occ_data=rgbif::occ_data(scientificName = sort(spnames$taxon_name)[i[1]:i[2]],country="US",hasCoordinate = TRUE, decimalLatitude = "22, 50",decimalLongitude = "-130, -60",limit=100000)
	occ_data <- spocc::occ(query = sort(spnames$taxon_name)[i[1]:i[2]],geometry = c(-130,22,-60,50),from=c("gbif","idigbio","inat"),limit=100000,has_coords=TRUE)	
	save(occ_data,file=paste("data analysis/occ_data.plants/",i[2],".Rdata",sep=""))
	rm(occ_data)
}

library(parallel)
no_cores <- 59#detectCores() - 1
mycl <- makePSOCKcluster(no_cores);
parSapply(cl=mycl,X=seql3,mf,spnames)
stopCluster(mycl)	

## step 1.2 --- data downloading of FL potiental invaders
library(dplyr)
#splist=read.csv("data analysis/horizon_scan_100_plant_list.csv")%>%filter(!Dis%in%"US native")
#Global Register of Introduced and Invasive Species - United States (Contiguous) (ver.2.0, 2022)
splist=read.csv("data analysis/US-RIIS.csv")%>%.[grep("invasive",.[,"degreeOfEstablishment"]),]%>%filter(kingdom%in%c("Animalia","Plantae"))
seql=c(seq(1,nrow(splist),100),nrow(splist))
seql2=lapply(2:length(seql),function(i,seql){if (i<length(seql)){c(seql[i-1],(seql[i]-1))
	} else {
	c(seql[i-1],seql[i])
	}},seql)
	
mf=function(i,splist){
	occ_data <- spocc::occ(query = sort(splist$scientificName)[i[1]:i[2]],geometry = c(-93,24.48,-66,48.3),from=c("gbif","idigbio","inat"),limit=100000,has_coords=TRUE)	
	save(occ_data,file=paste("data analysis/occ_data.invasive/",i[2],".Rdata",sep=""))
	rm(occ_data);gc()
}

library(parallel)
no_cores <- 40#detectCores() - 1
mycl <- makePSOCKcluster(no_cores);
parSapply(cl=mycl,X=seql2,mf,splist)
stopCluster(mycl)	

## step 2 --- data cleaning
library(dplyr);library(raster);library(sf)
ele.us <- raster("data analysis/wc2.1_10m_bio/wc2.1_10m_elev.tif")
rasterResolution <- max(res(ele.us))
#occlist <- list.files("data analysis/occ_data.plants", pattern = "*.Rdata", full.names = TRUE)
occlist <- list.files("data analysis/occ_data.invasive", pattern = "*.Rdata", full.names = TRUE)
clean.occdata=function(filename,rasterResolution,us=TRUE){
	require(dplyr)
	occ=get(load(filename))
	occ=subset(occ,names(occ)%in%c("gbif","idigbio","inat","bison"))
	dat.all=c()
	for (dat in 1:length(occ)){
		dat.t=vctrs::list_drop_empty(occ[[dat]][[2]])
		if (length(dat.t)>0) {
			dat.t2=do.call(rbind,lapply(1:length(dat.t),function(i,dat.t) data.frame(Sci_name=gsub("_"," ",names(dat.t)[i]),dat.t[[i]][,c("longitude","latitude")]),dat.t))
			dat.t2$latitude <- round(as.numeric(dat.t2$latitude), digits = 2);dat.t2$longitude <- round(as.numeric(dat.t2$longitude), digits = 2)
			dat.t2=unique(dat.t2)
			#dat.t2$scour=names(occ)[dat]
		}
		dat.all=rbind(dat.all,dat.t2)
	}
	dat.all=unique(dat.all)%>%na.omit()
	#if (us) dat.all=dat.all%>%filter(longitude>=-124.7258&(longitude<=-66.94989&latitude>=24.49813 & latitude<=49.38436))%>%na.omit()
	
	splist=unique(dat.all$Sci_name)
	clean.occ2=function(j,dat.all,rasterResolution){	
		edd.sp=dat.all%>%filter(Sci_name%in%j)
		# remove one point in a pair that occurs within one pixel
		while(min(spatstat.geom::nndist(edd.sp[,c("latitude","longitude")])) < rasterResolution){
		  nnD <- spatstat.geom::nndist(edd.sp[,c("latitude","longitude")])
		  edd.sp <- edd.sp[-(which(min(nnD) == nnD) [1]), ]
		}
		return(edd.sp)
	}
	# occ.clean.t=c()
	# for(i in splist){
	 # tmp=clean.occ2(i,dat.all,rasterResolution)
	# occ.clean.t=rbind(occ.clean.t,tmp)
	# print(i)
	# }
	
	occ.clean.t=do.call(rbind,lapply(splist,clean.occ2,dat.all,rasterResolution))
	return(occ.clean.t)
}

library(parallel)
no_cores <- 57#detectCores() - 1
mycl <- makePSOCKcluster(no_cores);
dis=parLapply(cl=mycl,X=occlist,clean.occdata,rasterResolution,us=FALSE)
dis2=do.call(rbind,dis)%>%distinct()
stopCluster(mycl)	
#save(dis2,file="data analysis/occ_data.plants/occ.us.RData")
save(dis2,file="data analysis/occ_data.invasive/occ.ena.invasive.RData")

# dis=list()
# for(i in 1:length(occlist)){
	# dis[[i]]=clean.occdata(occlist[i],rasterResolution,us=FALSE)
	# print(i)
# }

## step 3 --- get occurenc data from BIEN and do the data cleanning
library(BIEN)
spdis=BIEN_occurrence_country(country = "United States",cultivated =TRUE,native.status=TRUE)
save(spdis,file="data analysis/spdis.us.bien.rawdata.RData")#including introduced records
library(dplyr);
spdis2=spdis%>%filter(longitude>=-124.7258&longitude<=-66.94989&latitude>=24.49813 & latitude<=49.38436&(!is.na(scrubbed_species_binomial)))
spdis2$latitude <- round(as.numeric(spdis2$latitude), digits = 2);spdis2$longitude <- round(as.numeric(spdis2$longitude), digits = 2)
spdis2=spdis2%>%distinct()	
save(spdis2,file="data analysis/spdis.us.bien.cleangeo.RData")#including introduced records

spdis=BIEN_occurrence_country(country = "United States")#Exclude detected introduced species
spdis2=spdis%>%filter(longitude>=-124.7258&longitude<=-66.94989&latitude>=24.49813 & latitude<=49.38436&(!is.na(scrubbed_species_binomial)))
spdis2$latitude <- round(as.numeric(spdis2$latitude), digits = 2);spdis2$longitude <- round(as.numeric(spdis2$longitude), digits = 2)
spdis2=spdis2%>%select(scrubbed_species_binomial,latitude,longitude)%>%distinct()	
spdis=spdis2
save(spdis,file="data analysis/spdis.us.bien.cleangeo.native.RData")#Exclude detected introduced species

## step 3.5 --- correct spnames using POWO
load("data analysis/spdis.us.bien.cleangeo.native.RData")
scrubbed_species_binomial=unique(spdis$scrubbed_species_binomial)
spnames <- rWCVPdata::wcvp_names
splist.mat = as.data.frame(scrubbed_species_binomial) %>% left_join(spnames,by=c("scrubbed_species_binomial"="taxon_name"),relationship ="many-to-many") %>% filter(taxon_status%in%c("Accepted","Synonym")) %>% 
	dplyr::select(scrubbed_species_binomial,powo_id,plant_name_id,accepted_plant_name_id,taxon_rank,taxon_status,primary_author)
sp.unmat=as.data.frame(scrubbed_species_binomial) %>% filter(!scrubbed_species_binomial%in%unique(splist.mat$scrubbed_species_binomial)) 

library(TNRS)	
res=TNRS(taxonomic_names = sp.unmat$scrubbed_species_binomial,sources = "wcvp")
res.unmat=res%>%filter(Overall_score<0.9|Taxonomic_status=="No opinion")

res2=TNRS(taxonomic_names = res.unmat$Name_submitted)
res.all=rbind(res[!res$Name_submitted%in%res.unmat$Name_submitted,],res2)
#write.csv(res.all,"res.all.csv")#mannually check
sp.mat=splist.mat[order(splist.mat$taxon_status),]%>% .[!duplicated(.[,'scrubbed_species_binomial']),]%>%dplyr::select('scrubbed_species_binomial','taxon_status','accepted_plant_name_id')%>%
	left_join(spnames[,c('plant_name_id','taxon_name')],by=c('accepted_plant_name_id'='plant_name_id'))
res.all=res.all[,c("Name_submitted",'Taxonomic_status',"Accepted_name_id","Accepted_name")]
colnames(sp.mat)=colnames(res.all)
spname.cor=rbind(sp.mat,res.all)
#save(spname.cor,file="data analysis/BIENnames.POWO.Rdata")
load("data analysis/BIENnames.POWO.Rdata")
spdis.powo=spdis%>%left_join(spname.cor,by=c("scrubbed_species_binomial"="Name_submitted"))

## step 4 --- combine BIEN and GBIF data
load("data analysis/occ_data.plants/occ.us.RData")
spdis.powo=spdis.powo[,c("Accepted_name","longitude", "latitude")]%>%distinct()
colnames(spdis.powo)=colnames(dis2)
dis.all=rbind(dis2,spdis.powo)%>%distinct()

#data cleaning
library(raster)
ele.us <- raster("data analysis/wc2.1_10m_bio/wc2.1_10m_elev.tif")
US <- sf::read_sf('data analysis/cb_2018_us_state_20m/cb_2018_us_state_20m.shp')%>%filter(!NAME%in%c("Hawaii","Alaska","Puerto Rico"))%>%as(., "Spatial")
US=US%>%mask(ele.us,.)%>%crop(extent(US))
rasterResolution <- max(res(US))
splist=unique(dis.all$Sci_name)
	clean.occ2=function(j,dat.all,rasterResolution){	
		require(dplyr)
		edd.sp=dat.all%>%filter(Sci_name%in%j)
		# remove one point in a pair that occurs within one pixel
		while(min(spatstat.geom::nndist(edd.sp[,c("latitude","longitude")])) < rasterResolution){
		  nnD <- spatstat.geom::nndist(edd.sp[,c("latitude","longitude")])
		  edd.sp <- edd.sp[-(which(min(nnD) == nnD) [1]), ]
		}
		return(edd.sp)
	}
rm(dis2,spdis,spdis.powo,spname.cor,US,ele.us)	
library(parallel)
no_cores <- 59#detectCores() - 1
mycl <- makePSOCKcluster(no_cores);
occ.clean0=parLapply(cl=mycl,X=splist,fun=clean.occ2,dis.all,rasterResolution)
occ.clean=do.call(rbind,occ.clean0)
stopCluster(mycl)
save(occ.clean,file="data analysis/occ_data.plants/occ.us.plusBIEN.RData")

	
#combine BIEN data
library(dplyr)
load("data analysis/sp.try.inv.Rdata")
load("data analysis/spdis.us.bien.rawdata.RData")
spdis.inv=spdis%>%filter(scrubbed_species_binomial%in%sp.try.inv$species)

#data cleaning
splist=unique(dat.all$Sci_name)
	clean.occ2=function(j,dat.all,rasterResolution){	
		edd.sp=dat.all%>%filter(Sci_name%in%j)
		# remove one point in a pair that occurs within one pixel
		while(min(spatstat.geom::nndist(edd.sp[,c("latitude","longitude")])) < rasterResolution){
		  nnD <- spatstat.geom::nndist(edd.sp[,c("latitude","longitude")])
		  edd.sp <- edd.sp[-(which(min(nnD) == nnD) [1]), ]
		}
		return(edd.sp)
	}
	occ.clean.t=do.call(rbind,lapply(splist,clean.occ2,dat.all,0.5))

##############################################################
## create convex hull for each species and then intersect with raster layer
##############################################################
library(raster);library(dplyr)
#load("data analysis/Domain.us.rda")
load("data analysis/Domain.ena.rda")
ele.us <- raster("data analysis/wc2.1_10m_bio/wc2.1_10m_elev.tif")
#US <- sf::read_sf('data analysis/cb_2018_us_state_20m/cb_2018_us_state_20m.shp')%>%filter(!NAME%in%c("Hawaii","Alaska","Puerto Rico"))%>%as(., "Spatial")
#US=US%>%mask(ele.us,.)%>%crop(extent(US))
ENA.name=c("ME","NH","VT","NY","MA","RI","CT","NJ","PA","DE","MD","MI","OH","IN","IL","WI","WV","VA","NC","TN","KY","SC","GA","AL","MS","FL")
ENA <- sf::read_sf('data analysis/cb_2018_us_state_20m/cb_2018_us_state_20m.shp')%>%filter(STUSPS%in%ENA.name)%>%as(., "Spatial")
ENA.rst=ENA%>%raster::mask(ele.us,.)%>%raster::crop(extent(ENA))#%>%aggregate(fact=3)

# occur_points=get(load("data analysis/occ_data.plants/occ.us.plusBIEN.RData"))%>%filter(longitude>=-92.83333&longitude<=-66.83333&latitude>=24.16667 & latitude<=48.16667&(!Sci_name%in%""))
# splist=sort(unique(occur_points$Sci_name))

occur_points.inv=get(load("data analysis/occ_data.invasive/occ.all.RData"))#%>%filter(longitude>=-92.83333&longitude<=-66.83333&latitude>=24.16667 & latitude<=48.16667&(!Sci_name%in%""))
splist.inv=sort(unique(occur_points.inv$Sci_name))

get.frag=function(species,occur_points,US,domain.us){
	librariesToload=c("sp","dplyr","rangeBuilder","geosphere","raster","sf","graph")
	sapply(librariesToload,require,character=TRUE)
	#species=splist[i]
	print(species)	
	edd.geo=occur_points%>%filter(Sci_name%in%species)%>%as.data.frame();rownames(edd.geo)=NULL
	alldfsp <- SpatialPointsDataFrame(edd.geo[,c("longitude","latitude")], data = edd.geo, proj4string=CRS("EPSG:4326"))
	## STEP 2-2.1 cluster occurence points based on the shortest distance between two points on an ellipsoid (default is WGS84 ellipsoid) -----
	if(length(alldfsp)>2&length(unique(alldfsp$longitude))>1&length(unique(alldfsp$latitude))>1){
		mdist <- distm(alldfsp) # generate a geodesic distance matrix in meters
		# recl <- recluster::recluster.region(mdist, method="average", mincl=4, maxcl=8)
		# nclust=min(recl[[1]][recl[[1]][,"ex.diss"]>=0.8,"k"])
		# cluster all points using a hierarchical clustering approach
		hc <- hclust(as.dist(mdist), method="average")
		# define clusters based on a tree "height" cutoff "distance thres" and add them to the SpDataFrame
		alldfsp$clust <- cutree(hc, h=300000) #distance thres =300km (ca.18 cells)
		#plot(alldfsp, col=rainbow(length(unique(alldfsp$clust)))[factor(alldfsp$clust)])

		## STEP 2-2.2 Find the convex hull of each cluster of points being plotted, and take the 50th quantile of the max distance between points as buffer size
		## STEP 2-2.3 mask with grid cells 
		get.buffer=function(cluster,alldfsp,US){
			dfsp=subset(alldfsp,alldfsp$clust==cluster)
			#Create alpha hull
			if(length(unique(dfsp$longitude))==1|length(unique(dfsp$latitude))==1) {
					re=buffer(dfsp,20000)%>%mask(US,.)%>%raster::as.data.frame(xy=TRUE)
			}else{
				result <- try({hull <- getDynamicAlphaHull(x = dfsp@coords,clipToCoast = "no")}, silent = TRUE)
				if (class(result)== "try-error"){
					re=buffer(dfsp,20000)%>%mask(US,.)%>%raster::as.data.frame(xy=TRUE)
				}else{				
					hull <- getDynamicAlphaHull(x = dfsp@coords,clipToCoast = "no")
					hullTrans <- st_transform(hull[[1]], "+proj=cea +lat_ts=0 +lon_0=0")
					alldfspTrans <- spTransform(dfsp, "+proj=cea +lat_ts=0 +lon_0")
					#buffDist <- quantile(x = (apply(spDists(alldfspTrans), 2, FUN = function(x) sort(x)[2])), probs = 0.80, na.rm = TRUE) #Here we take the 80th quantile of the max distance between points as buffer size
					buffDist=20000 #ca. 1 cells
					buffer <- st_buffer(hullTrans,buffDist,endCapStyle="SQUARE")%>%st_transform(CRS("EPSG:4326"))
					#plot(dfsp); plot(hull[[1]],add=T);plot(buffer,add=T)		
					re=as(buffer, "Spatial")%>%mask(US,.)%>%raster::as.data.frame(xy=TRUE)
				}
			}				
			colnames(re)[3]="spdis"
			re=re%>%left_join(domain.us,by=c("x","y"))%>%filter(!is.na(spdis))
			return(re$domain.ID)
		}	
		frags = lapply(unique(alldfsp$clust),get.buffer,alldfsp,US)
		# for(cluster in unique(alldfsp$clust)){
			# tmp=get.buffer(cluster,alldfsp,US)
			# print(cluster)
		# }
	}else{
		re=buffer(alldfsp,20000)%>%mask(US,.)%>%raster::as.data.frame(xy=TRUE)
		colnames(re)[3]="spdis"
		re=re%>%left_join(domain.us,by=c("x","y"))%>%filter(!is.na(spdis))
		frags=as.list(re$domain.ID)
	}
	## STEP 2-2.4 Merging Listed Vectors that share Elementsand then make identify connected sets of cells (i.e., fragments)
	input <- frags
	repeat {
	  # Get A Count Table
	  tbl <- table(unlist(input))
	  # No Duplicated Items Then break Out
	  if (length(tbl[tbl > 1]) == 0) { break }
	  # Take A First Duplicated Item And Get Index Where The Item Is Located
	  idx <- which(sapply(seq_len(length(input)), function(i) {
		any(input[[i]] == names(tbl[tbl > 1])[1])
	  }))
	  # Create New vector By Union
	  newvec <- sort(unique(unlist(input[idx])))
	  # Append newvec To cbnl And Remove Original vectors
	  input <- c(input, list(newvec))[-idx]
	}
	names(input)=rep(species,length(input))
	return(input)	
}

library(parallel)
no_cores <- 57
mycl <- makePSOCKcluster(no_cores);
#RangeFragments=parLapply(cl=mycl,X=splist,fun=get.frag,occur_points,US,domain.us)
#RangeFragments=parLapply(cl=mycl,X=splist,fun=get.frag,occur_points,ENA.rst,domain.ena)
RangeFragments=parLapply(cl=mycl,X=splist.inv,fun=get.frag,occur_points.inv,ENA.rst,domain.ena)#potienal invaders, 81 species
stopCluster(mycl)

#debug
RangeFragments=list()
for(i in 1:length(splist.inv)){
 tmp=get.frag(splist.inv[i],occur_points.inv,ENA.rst,domain.ena)
	#check distribution
	# length(unlist(tmp));length(unique(unlist(tmp)))
	# summary(tmp)
	plot(domain.ena$x,domain.ena$y,cex=0.5)
	id=unlist(tmp)
	#id=tmp[[1]]
	focalDom<-domain.ena%>%filter(domain.ID%in%id)
	points(focalDom$x,focalDom$y,col="green",pch=16,cex=0.5) 
	edd.geo=occur_points.inv%>%filter(Sci_name%in%splist.inv[i])
	points(edd.geo$longitude,edd.geo$latitude,col="red",pch=16,cex=0.5)
 RangeFragments=c(RangeFragments,tmp)
}

#names(RangeFragments)=splist;
names(RangeFragments)=splist.inv;
get.fin=function(RangeFragments){
	RangeFragments2=vctrs::list_drop_empty(RangeFragments)#some points are in the ocean and  result in no match with grids and thus may induce frag of 0 length
	name=do.call(c,lapply(names(RangeFragments2),function(i,RangeFragments2) rep(i,length(RangeFragments2[[i]])),RangeFragments2))
	RangeFragments.fin=do.call(c,RangeFragments2);names(RangeFragments.fin)=name
	#RangeFragments.us=vctrs::list_drop_empty(RangeFragments.fin)
	#save(RangeFragments.us,file="data analysis/NativeRangeFragments.plants.US.RData")
	RangeFragments.ena=vctrs::list_drop_empty(RangeFragments.fin)
	return(RangeFragments.ena)
}
#RangeFragments.ena=get.fin(RangeFragments)
RangeFragments.ena.inv=get.fin(RangeFragments)
# save(RangeFragments.ena,file="data analysis/NativeRangeFragments.plants.ENA.RData")
# save(RangeFragments.ena.inv,file="data analysis/AlienRangeFragments.plants.ENA.RData")
save(RangeFragments.ena.inv,file="data analysis/AlienRangeFragments.alltaxa.ENA.RData")
##############################################################
#########remove introduced records based on powo #############
##############################################################
## POWO
library(rWCVP);library(dplyr)
library(sp);library(sf);library(raster);library(data.table)
#load("data analysis/Domain.us.rda")
load("data analysis/Domain.ena.rda")
dis <- rWCVPdata::wcvp_distributions
dis.us=dis[grep("U.S.A.",dis$region),]%>% filter(!is.na(area_code_l3))
# ele.us <- raster("data analysis/wc2.1_10m_bio/wc2.1_10m_elev.tif")
# grids=rasterFromXYZ(domain.us[,c("x","y","domain.ID")],res=res(ele.us),crs=crs(ele.us))%>%stars::st_as_stars()%>%st_as_sf()#%>% mutate(area.grids = st_area(.) %>% as.numeric())
ele.us <- raster("data analysis/wc2.1_10m_bio/wc2.1_10m_elev.tif")#%>%aggregate(fact=3)
grids=rasterFromXYZ(domain.ena[,c("x","y","domain.ID")],res=res(ele.us),crs=crs(ele.us))%>%stars::st_as_stars()%>%st_as_sf()#%>% mutate(area.grids = st_area(.) %>% as.numeric())

powo.us <- sf::read_sf("data analysis/WGSRPD-level3/level3.shp")%>% filter(LEVEL3_COD%in%dis.us$area_code_l3)%>%
		st_transform(crs=st_crs(grids))%>%st_make_valid()%>%  mutate(area.powo = st_area(.) %>% as.numeric())
int <- st_intersection(grids, powo.us) 
powo2grid=int%>%as.data.frame()%>%dplyr::select(domain.ID,LEVEL3_COD)%>%distinct()

#RangeFragments<-get(load("data analysis/NativeRangeFragments.plants.US.RData"))
RangeFragments<-get(load("data analysis/NativeRangeFragments.plants.ENA.RData"))
splist=unique(names(RangeFragments))
spnames <- rWCVPdata::wcvp_names%>%filter(taxon_name%in%splist)
spnames=spnames[order(spnames$taxon_status),]%>% .[!duplicated(.[,'taxon_name']),]%>%dplyr::select(plant_name_id,taxon_name)%>%distinct()
distributions <- rWCVPdata::wcvp_distributions%>%filter(plant_name_id%in%spnames$plant_name_id)%>%dplyr::select(plant_name_id,introduced,area_code_l3)%>%distinct()
disPOWO=distributions%>%left_join(spnames,by="plant_name_id",relationship =  "many-to-many")%>%filter(introduced==1&(!is.na(taxon_name)))%>%
	left_join(powo2grid,by=c("area_code_l3"="LEVEL3_COD"),relationship ="many-to-many")%>%
	filter(!is.na(domain.ID))%>%dplyr::select(taxon_name,domain.ID,introduced)%>%distinct()

## glonaf
region.glonaf=read.csv("data analysis/Region_gloNAF.csv")
taxa.glonaf=as.data.frame(fread("data analysis/Taxon_gloNAF.csv"))
glonaf=taxa.glonaf%>%left_join(region.glonaf,by="region_id")%>%filter(country%in%"United States of America (the)")%>%
	dplyr::select(standardized_name,tdwg3)%>%filter(!tdwg3%in%c("ASK","HAW"))%>%distinct()%>%
	left_join(powo2grid,by=c("tdwg3"="LEVEL3_COD"),relationship ="many-to-many")%>%
	filter(!is.na(domain.ID))	
#using POWO to correct names
trenames=glonaf%>%dplyr::select(standardized_name)%>%distinct()
spnames=rWCVPdata::wcvp_names
splist.mat =trenames %>% left_join(spnames,by=c("standardized_name"="taxon_name"),relationship = "many-to-many") %>% filter(taxon_status%in%c("Accepted","Synonym")) %>% 
	dplyr::select(standardized_name,accepted_plant_name_id,taxon_rank,taxon_status)
sp.unmat=trenames %>% filter(!standardized_name%in%splist.mat$standardized_name) 
library(TNRS)	
res=TNRS(taxonomic_names = sp.unmat$standardized_name,sources = "wcvp")
res.unmat=res%>%filter(Overall_score<0.9|Taxonomic_status=="No opinion"|Name_matched_rank%in%"genus")
res2=TNRS(taxonomic_names = res.unmat$Name_submitted)
res.all=rbind(res[!res$Name_submitted%in%res.unmat$Name_submitted,],res2)%>%filter((!Taxonomic_status%in%"No opinion")&(!Name_matched_rank%in%"genus")&Genus_score==1&Overall_score>=0.9)

sp.mat=splist.mat[order(splist.mat$taxon_status),] %>% .[!duplicated(.[,'standardized_name']),]%>%dplyr::select('standardized_name','taxon_status','accepted_plant_name_id')%>%
	left_join(spnames[,c('plant_name_id','taxon_name')],by=c('accepted_plant_name_id'='plant_name_id'))
res.fin=res.all%>%dplyr::select(Name_submitted,Taxonomic_status,Accepted_name_id,Accepted_name)
colnames(sp.mat)=colnames(res.fin)
spname.cor=rbind(sp.mat,res.fin)
glonaf=glonaf%>%left_join(spname.cor,by=c("standardized_name"="Name_submitted"))%>%dplyr::select(Accepted_name,domain.ID)%>%distinct()
glonaf$introduced=1;colnames(glonaf)[1]="taxon_name"

introduce.us=rbind(glonaf,disPOWO)%>% distinct()
save(introduce.us,file="data analysis/introducedData.US.Rdata")
save(introduce.us,file="data analysis/introducedData.ENA.Rdata")

load("data analysis/introducedData.ENA.Rdata")
RangeFragments<-get(load("data analysis/NativeRangeFragments.plants.ENA.RData"))
#RangeFragments<-get(load("data analysis/NativeRangeFragments.plants.US.RData"))
mf=function(i,RangeFragments) {return(data.frame(Tip_Label=names(RangeFragments)[i],domain.ID=RangeFragments[[i]],fragNum=i))}
rangeData=do.call(rbind,lapply(1:length(RangeFragments),mf,RangeFragments))%>%
	left_join(introduce.us,by=c("Tip_Label"="taxon_name","domain.ID"))%>%filter(is.na(introduced))
#correct fragNum
frag=rangeData%>%dplyr::select(Tip_Label)%>%distinct()
frag$fragNum=1:nrow(frag)
rangeData=rangeData[,-c(3:4)]%>%left_join(frag,by="Tip_Label")
#save(rangeData,file="data analysis/NativeRangeFragments.plants.US.clean.RData")
save(rangeData,file="data analysis/NativeRangeFragments.plants.ENA.clean.RData")

#################################################
############ get trait data from BIEN ###############
#################################################
# db_version db_release_date
#       4.2.6      2022-10-19
library(BIEN);library(dplyr)
bien.sp=BIEN_list_country(country = "United States")
#traitname=BIEN_trait_list()
trait.selct=c("leaf area per leaf dry mass","leaf nitrogen content per leaf dry mass","maximum whole plant height","stem wood density","whole plant height")
traits=list()
for (i in trait.selct) traits[[i]]=BIEN_trait_mean(species=bien.sp$scrubbed_species_binomial,trait=i)
save(traits,file="data analysis/traits.us.bien.RData")

## data cleaning and correct names
## correct spnames using POWO
spnames <- rWCVPdata::wcvp_names
splist.mat = bien.sp %>% left_join(spnames,by=c("scrubbed_species_binomial"="taxon_name"),relationship ="many-to-many") %>% filter(taxon_status%in%c("Accepted","Synonym")) %>% 
	dplyr::select(scrubbed_species_binomial,powo_id,plant_name_id,accepted_plant_name_id,taxon_rank,taxon_status,primary_author)
sp.unmat=bien.sp %>% filter(!scrubbed_species_binomial%in%unique(splist.mat$scrubbed_species_binomial)) 

library(TNRS)	
res=TNRS(taxonomic_names = sp.unmat$scrubbed_species_binomial,sources = "wcvp")
res.unmat=res%>%filter(Overall_score<0.9|Genus_score<1|Taxonomic_status=="No opinion"|Name_matched_rank%in%c("","genus"))

res2=TNRS(taxonomic_names = res.unmat$Name_submitted)
res.all=rbind(res[!res$Name_submitted%in%res.unmat$Name_submitted,],res2)
#write.csv(res.all,"res.all.csv")#mannually check
sp.mat=splist.mat[order(splist.mat$taxon_status),]%>% .[!duplicated(.[,'scrubbed_species_binomial']),]%>%dplyr::select('scrubbed_species_binomial','taxon_status','accepted_plant_name_id')%>%
	left_join(spnames[,c('plant_name_id','taxon_name')],by=c('accepted_plant_name_id'='plant_name_id'))
res.fin=res.all%>%filter((!Name_matched_rank%in%c("","genus"))&(!Taxonomic_status%in%"No opinion"))%>%select(Name_submitted,Taxonomic_status,Accepted_name_id,Accepted_name)
colnames(sp.mat)=colnames(res.fin)
spname.cor=rbind(sp.mat,res.fin)
save(spname.cor,file="data analysis/BIENtrait.names.POWO.Rdata")

load("data analysis/BIENtrait.names.POWO.Rdata")
load("data analysis/traits.us.bien.RData")
mf=function(trt,spname.cor) trt.powo=trt%>%left_join(spname.cor,by=c("species"="Name_submitted"))
traits.powo=lapply(traits,mf,spname.cor)
save(traits.powo,file="data analysis/traits.us.bien.powo.RData")

#trait data for potienal invaders
bien.sp=BIEN_list_all()
## correct spnames using POWO
spnames <- rWCVPdata::wcvp_names
splist.mat = bien.sp %>% left_join(spnames,by=c("species"="taxon_name"),relationship ="many-to-many") %>% filter(taxon_status%in%c("Accepted","Synonym")) %>% 
	dplyr::select(species,accepted_plant_name_id,taxon_rank,taxon_status)
sp.unmat=bien.sp %>% filter(!species%in%unique(splist.mat$species)) 

library(TNRS)	
res=TNRS(taxonomic_names = sp.unmat$species,sources = "wcvp")
res.unmat=res%>%filter(Overall_score<0.9|Genus_score<1|Taxonomic_status=="No opinion"|Name_matched_rank%in%c("","genus"))

res2=TNRS(taxonomic_names = res.unmat$Name_submitted)
res.all=rbind(res[!res$Name_submitted%in%res.unmat$Name_submitted,],res2)
a=res.all%>%filter((Overall_score<0.9&Genus_score==1&(!Name_matched_rank%in%c("","genus"))))
write.csv(a,"a.csv")#mannually check
b=res.all%>%filter((Overall_score>=0.9&Genus_score==1&(!Name_matched_rank%in%c("","genus"))))
res.all=rbind(read.csv("a.csv"),b)

sp.mat=splist.mat[order(splist.mat$taxon_status),]%>% .[!duplicated(.[,'species']),]%>%dplyr::select('species','taxon_status','accepted_plant_name_id')%>%
	left_join(spnames[,c('plant_name_id','taxon_name')],by=c('accepted_plant_name_id'='plant_name_id'))
res.fin=res.all%>%select(Name_submitted,Taxonomic_status,Accepted_name_id,Accepted_name)
colnames(res.fin)=colnames(sp.mat)
spname.cor=rbind(sp.mat,res.fin)
splist=read.csv("data analysis/horizon_scan_100_plant_list.csv")%>%filter(!Dis%in%"US native")
sp.try.inv=spname.cor%>%filter(taxon_name%in%splist$acc_name)
save(sp.try.inv,file="sp.try.inv.Rdata")
trait.selct=c("leaf area per leaf dry mass","leaf nitrogen content per leaf dry mass","maximum whole plant height","stem wood density","whole plant height")
traits=list()
for (i in trait.selct) traits[[i]]=BIEN_trait_mean(species=sp.try.inv$species,trait=i)
save(traits,file="data analysis/traits.invasive.bien.RData")

#combine trait data of native and potienal invaders
traits.nat=get(load("data analysis/traits.us.bien.powo.RData"))
traits.inv=get(load("data analysis/traits.invasive.bien.RData"))
load("data analysis/BIENtrait.names.POWO.Rdata")
load("data analysis/sp.try.inv.Rdata")
colnames(spname.cor)=colnames(sp.try.inv)
spnames=rbind(sp.try.inv,spname.cor)%>%distinct()
mf=function(i,spnames) {re=i%>%left_join(spnames,by="species",relationship="many-to-many")%>%select(taxon_name,mean_value,trait,unit)%>%distinct();return(re)}
biendata=rbind(do.call(rbind,lapply(traits.inv,mf,spnames)),do.call(rbind,lapply(traits.nat,mf,spnames)))%>%distinct()
rownames(biendata)=NULL
save(biendata,file="data analysis/traits.bien.powo.RData")

#################################################
############ get trait data from TRY ###############
#################################################
try.sp=read.delim("data analysis/TryAccSpecies.txt",header=T)
## correct spnames using POWO
spnames <- rWCVPdata::wcvp_names
splist.mat = try.sp %>% left_join(spnames,by=c("AccSpeciesName"="taxon_name"),relationship ="many-to-many") %>% filter(taxon_status%in%c("Accepted","Synonym")) %>% 
	dplyr::select(AccSpeciesName,AccSpeciesID,plant_name_id,accepted_plant_name_id,taxon_rank,taxon_status)
sp.unmat=try.sp %>% filter(!AccSpeciesName%in%unique(splist.mat$AccSpeciesName)) 

library(TNRS)	
res=TNRS(taxonomic_names = sp.unmat$AccSpeciesName,sources = "wcvp")
res.unmat=res%>%filter(Overall_score<0.9|Genus_score<1|Taxonomic_status=="No opinion"|Name_matched_rank%in%c("","genus"))

res2=TNRS(taxonomic_names = res.unmat$Name_submitted)
res.all=rbind(res[!res$Name_submitted%in%res.unmat$Name_submitted,],res2)
a=res.all%>%filter((Overall_score<0.9&Genus_score==1&(!Name_matched_rank%in%c("","genus"))))
write.csv(a,"a.csv")#mannually check
b=res.all%>%filter((Overall_score>=0.9&Genus_score==1&(!Name_matched_rank%in%c("","genus"))))
res.all=rbind(read.csv("a.csv"),b)

sp.mat=splist.mat[order(splist.mat$taxon_status),]%>% .[!duplicated(.[,'AccSpeciesName']),]%>%dplyr::select('AccSpeciesName','taxon_status','accepted_plant_name_id')%>%
	left_join(spnames[,c('plant_name_id','taxon_name')],by=c('accepted_plant_name_id'='plant_name_id'))
res.fin=res.all%>%select(Name_submitted,Taxonomic_status,Accepted_name_id,Accepted_name)
colnames(res.fin)=colnames(sp.mat)
spname.cor=rbind(sp.mat,res.fin)
trp.sp.powo=try.sp%>%left_join(spname.cor,by="AccSpeciesName")
save(trp.sp.powo,file="data analysis/TryAccSpecies.POWO.Rdata")

load("data analysis/TryAccSpecies.POWO.Rdata")
occur_points=get(load("data analysis/occ_data.plants/occ.us.plusBIEN.RData"))%>%filter(longitude>=-92.83333&longitude<=-66.83333&latitude>=24.16667 & latitude<=48.16667&(!Sci_name%in%""))
splist1=sort(unique(occur_points$Sci_name))
splist2=read.csv("data analysis/horizon_scan_100_plant_list.csv")%>%filter(!Dis%in%"US native")
splist=unique(c(splist2$acc_name,splist1))
sp.in.ena=trp.sp.powo%>%filter(AccSpeciesName%in%splist)
seql=c(seq(1,nrow(sp.in.ena),9990),nrow(sp.in.ena))
seql2=lapply(2:length(seql),function(i,seql){if (i<length(seql)){c(seql[i-1],(seql[i]-1))
	} else {
	c(seql[i-1],seql[i])
	}},seql)
for(i in 1:length(seql2)){
	spname=paste(paste(sp.in.ena$AccSpeciesID[seql2[[1]][1]:seql2[[1]][2]],collapse=","),",",sep="")
	write.table(spname,paste("data analysis/try.query.name",i,"txt",sep="."),quote=FALSE,row.names = FALSE,col.names = FALSE)
}

#ref:Kattge, J, Bönisch, G, Díaz, S, et al. TRY plant trait database – enhanced coverage and open access. Glob Change Biol. 2020; 26: 119– 188. https://doi.org/10.1111/gcb.14904
require(data.table);library(dplyr)
TRYdata <- fread("data analysis/TRY data/30740.txt", header = T, sep = "\t", dec = ".", quote = "", data.table = T)
TRYdata2=TRYdata%>%filter(AccSpeciesID%in%sp.in.ena$AccSpeciesID&(!is.na(TraitID)))%>%
	select(AccSpeciesID,TraitID,TraitName,DataName,OriglName,OrigValueStr,OrigUnitStr,StdValue,UnitName,Reference)%>%
	distinct()%>%left_join(sp.in.ena[,c("AccSpeciesID","taxon_name")],by="AccSpeciesID")
TRYdata <- fread("data analysis/TRY data/30421.txt", header = T, sep = "\t", dec = ".", quote = "", data.table = T)
TRYdata.clean=TRYdata%>%filter(AccSpeciesID%in%sp.in.ena$AccSpeciesID&(!is.na(TraitID)))%>%
	select(AccSpeciesID,TraitID,TraitName,DataName,OriglName,OrigValueStr,OrigUnitStr,StdValue,UnitName,Reference)%>%
	distinct()%>%left_join(sp.in.ena[,c("AccSpeciesID","taxon_name")],by="AccSpeciesID")
TRYdata2=rbind(TRYdata2,TRYdata.clean)%>%distinct()
save(TRYdata2,file="data analysis/traits.try.powo.RData")

#################################################
############ merge all trait data ###############
#################################################
##combine try and bien data
library(dplyr)
load("data analysis/traits.bien.powo.RData")
load("data analysis/traits.try.powo.RData")
bien.trat=unique(biendata[,c("trait", "unit")])
biendata[biendata$trait=="maximum whole plant height","trait"]="whole plant height"
try.trat=unique(TRYdata2[,c("TraitID","TraitName","UnitName")])
sla=try.trat%>%filter(UnitName=="mm2 mg-1");sla$trait="leaf area per leaf dry mass"
ln=try.trat%>%filter(UnitName=="mg/g");ln$trait="leaf nitrogen content per leaf dry mass"
height=try.trat%>%filter(UnitName=="m"&(!TraitName%in%"Root rooting depth"));height$trait="whole plant height"
rootdep=try.trat%>%filter(TraitName%in%"Root rooting depth");rootdep$trait="root depth"
wd=try.trat%>%filter(UnitName=="g/cm3");wd$trait="stem wood density"
traitnames=rbind(sla,ln,height,rootdep,wd)
TRYdata=TRYdata2%>%left_join(traitnames,by="TraitID")%>%select(taxon_name,StdValue,trait)%>%distinct()%>%na.omit()
colnames(biendata)[2]="StdValue"
alldata=rbind(TRYdata,biendata[-4])%>%distinct()

## other trait data
## max root depth from Tumber-Dávila, S.J., et. al. (2022). New Phytol.
#Update WD using FIA’s REF_SPECIES table (the variable is called “WOOD_SPGR_GREENVOL_DRYWT”)
Fia.sp=read.csv("data analysis/FIAsp_clean.csv")
wd=read.csv("data analysis/REF_SPECIES.csv")%>% select(SPCD,WOOD_SPGR_GREENVOL_DRYWT)%>%left_join(Fia.sp,by="SPCD")%>%na.omit()%>%select(Spname,WOOD_SPGR_GREENVOL_DRYWT)
root.dep=read.csv("data analysis/nph18031-sup-0001-datasets1.csv")%>% select(Species,Dr)%>%filter(!(is.na(Dr)))
wd$trait="stem wood density";root.dep$trait="root depth"
colnames(root.dep)=colnames(wd)=c("Spname","StdValue","trait")
alldata2=rbind(wd,root.dep)

#correct names
sp=data.frame(AccSpeciesName=unique(alldata2$Spname))
spnames <- rWCVPdata::wcvp_names
splist.mat = sp %>% left_join(spnames,by=c("AccSpeciesName"="taxon_name"),relationship ="many-to-many") %>% filter(taxon_status%in%c("Accepted","Synonym")) %>% 
	dplyr::select(AccSpeciesName,plant_name_id,accepted_plant_name_id,taxon_rank,taxon_status)
sp.unmat=sp %>% filter(!AccSpeciesName%in%unique(splist.mat$AccSpeciesName)) 
write.csv(sp.unmat$AccSpeciesName,"sp.unmat.csv")

res=read.csv("tnrs_result.csv")
sp.mat=splist.mat[order(splist.mat$taxon_status),]%>% .[!duplicated(.[,'AccSpeciesName']),]%>%dplyr::select('AccSpeciesName','taxon_status','accepted_plant_name_id')%>%
	left_join(spnames[,c('plant_name_id','taxon_name')],by=c('accepted_plant_name_id'='plant_name_id'))
res.fin=res%>%select(Name_submitted,Taxonomic_status,Accepted_name_id,Accepted_name)
colnames(res.fin)=colnames(sp.mat)
spname.cor=rbind(sp.mat,res.fin)
alldata2=alldata2%>%left_join(spname.cor,by=c("Spname"="AccSpeciesName"))
traitdata=rbind(alldata2[,c("taxon_name","StdValue","trait")],alldata)%>%distinct()%>%na.omit()
traitdata$StdValue=as.numeric(traitdata$StdValue)
traitdata=traitdata%>%na.omit()%>%group_by(taxon_name,trait)%>%summarize(StdValue=mean(StdValue))
save(traitdata,file="data analysis/traits.all.RData")

#imputation missing values and filter ENA species
trait_mat0=traitdata %>% tidyr::pivot_wider(names_from = taxon_name, values_from = StdValue) %>% tibble::column_to_rownames("trait") %>%t()
rangeData=get(load("data analysis/NativeRangeFragments.plants.ENA.clean.RData"))%>%filter(Tip_Label%in%rownames(trait_mat0))
trait_mat0=subset(trait_mat0,rownames(trait_mat0)%in%unique(rangeData$Tip_Label))
#imputation missing values
library(Rphylopars);library(ape);
tree0=get(load("data analysis/Phylo.us.Rdata"))
trait_mat=subset(trait_mat0,rownames(trait_mat0)%in%tree0$tip)
trait_mat=trait_mat[,-5]#too many NAs (root dep: 6068 out of 8044)
tree=keep.tip(tree0,tip=rownames(trait_mat))

trait_data=as.data.frame(trait_mat)
trait_data$species=rownames(trait_mat)
trait_data=trait_data[,c("species",colnames(trait_mat))]
trait_mat1=phylopars(trait_data ,tree,
    pheno_error = TRUE,phylo_correlated = TRUE,pheno_correlated = TRUE)
trait_mat2=	trait_mat1$anc_recon[rownames(trait_mat),]
trait_mat.nophy=na.omit(subset(trait_mat0,!rownames(trait_mat0)%in%tree0$tip)[,-5])
trait_mat3=scale(rbind(trait_mat2,trait_mat.nophy))
colnames(trait_mat3)=c("SLA","LN","WD","Hmax")
save(trait_mat3,file="data analysis/traits.ENA.RData")

##using deep learning method to predict the porb of invasion given the similarity of a pair of grid Cell
## ref: https://saturncloud.io/blog/deep-learning-with-r/
## for each species,
## first, we randomly pick up a starting cell and caculated the likelyhood of invasion to its neighbour cell based on similarity
## then we repeated the above process for all the invaded cells, and overlay their likelyhood maps, the invasion likehood of each cell is the product multipled by all maps
## third, we link the likelyhood 
## made projections for other uninvaded cells. and then we repeated this
## how to intergrete dispersal process into this process?

################################################################################################################
# use native ranges to generate a community similarity matrix for calculating environmental resistance
################################################################################################################
## ER1: species assemblages ----
# generate a list where each element is a vector of the species in each grid cell
get.sim.sp=function(rangeData){
	speciesByCell<-list()
	speciesByCellsubset<-split(rangeData$Tip_Label,rangeData$domain.ID)
	domain.ID<-sort(unique(rangeData$domain.ID))
	speciesByCell[domain.ID]<-speciesByCellsubset

	# generate a matrix to store community similarity between each pair of cells
	nCells<-length(speciesByCell)
	CommSimilarity<-matrix(ncol=nCells,nrow=nCells)
	cellIN<-sapply(speciesByCell,length)  # number of species per cell
	for(i in 1:nCells){ # loop over each cell
	  focalN<-cellIN[i] # richness of cell i
	  if(focalN>0){
		cellIspecies<-speciesByCell[[i]] # species in cell i 
		for(j in 1:nCells){
		  if(cellIN[j]>0){
			cellJspecies<-speciesByCell[[j]] # species in  cell j 
			# proportion of species in cell i that are also in cell j 
			CommSimilarity[i,j]<-length(intersect(cellIspecies,cellJspecies))/focalN  
		  }else{
			CommSimilarity[i,j]<-0 # if the richness of cell j is zero, set simmilarity to zero
		  }
		}
	  }else{
		CommSimilarity[i,]<-0 # if the richness of cell i is zero, set simmilarity to zero
	  }
	  CommSimilarity[i,i]<-NA # set comparisons between cell i and i to NA
		print(i)
	}
	return(CommSimilarity)
}

rangeData=get(load("data analysis/NativeRangeFragments.plants.ENA.clean.RData"))
CommSimilarity.sp=get.sim.sp(rangeData)
save(CommSimilarity.sp,file=paste0("data analysis/CommSimilarity.ENA.sp.rda"))

#species with phylo infromation
library(ape);library(dplyr)	
#load("data analysis/Phylo.tipnames.POWO.Rdata")
tree0=get(load("data analysis/Phylo.us.Rdata"))
rangeData.phy=rangeData%>%filter(Tip_Label%in%tree0$tip)
CommSimilarity.sp2=get.sim.sp(rangeData.phy)
save(CommSimilarity.sp2,file=paste0("data analysis/CommSimilarity.ENA.sp2.rda"))

################################################################################################################
# use native ranges to generate a phylogenetic similarity matrix for calculating environmental resistance
################################################################################################################
## ER2: phylogenetic similarity ----
##using POWO to correct names
library(ape);library(data.table);library(dplyr)
tre0=read.tree("data analysis/FromAo_Smith_100treesV2/ALLMB.tre")		#obtained:Nov 7,2023
tip=data.frame(tip=sort(tre0$tip))
tip$Taxon_name=gsub("_", " ",tip$tip)
write.table(tip,"data analysis/FromAo_Smith_100treesV2/tipnames.txt")

# tre0=read.tree("data analysis/NorthAmerica.ExaML.treePL.best.tre")#Park DS,2020. Replication code and data, Zenodo. doi: 10.5281/zenodo.3755913.
# tip=data.frame(tip=sort(tre0$tip))
# tip$Taxon_name=gsub("_", " ",tip$tip)
# write.table(tip,"data analysis/tipnames.park.txt")

trenames=fread("data analysis/FromAo_Smith_100treesV2/tipnames.csv")
spnames <- rWCVPdata::wcvp_names
splist.mat = trenames %>% left_join(spnames,by=c("Taxon_name"="taxon_name"),relationship = "many-to-many") %>% filter(taxon_status%in%c("Accepted","Synonym")) %>% 
	select(tip,Taxon_name,accepted_plant_name_id,taxon_rank,taxon_status)
acc=spnames%>% filter(taxon_status%in%"Accepted")
trenames.powo=splist.mat%>% left_join(acc[,c("accepted_plant_name_id","taxon_name")],by="accepted_plant_name_id",relationship =  "many-to-many")%>%distinct()
sp.unmat=trenames %>% filter(!Taxon_name%in%unique(splist.mat$Taxon_name)) 

library(TNRS)	
res=TNRS(taxonomic_names = sp.unmat$Taxon_name,sources = "wcvp")
res.unmat=res%>%filter(Overall_score<0.9|Taxonomic_status=="No opinion")
res2=TNRS(taxonomic_names = res.unmat$Name_submitted)
res.all=rbind(res[!res$Name_submitted%in%res.unmat$Name_submitted,],res2)%>%filter((!Taxonomic_status%in%"No opinion")&(!Name_matched_rank%in%"genus"))
a=res.all%>%filter(Overall_score<0.9|Genus_score<1)
write.csv(a,"a.csv")#mannually check
a=read.csv("a.csv")
b=res.all%>%filter(Overall_score>=0.9&Genus_score==1)%>%rbind(read.csv("a.csv"))

sp.mat=trenames.powo[order(trenames.powo$taxon_status),]%>% filter(!is.na(accepted_plant_name_id)) %>% .[!duplicated(.[,'Taxon_name']),]%>%dplyr::select('tip','Taxon_name','taxon_status','accepted_plant_name_id')%>%
	left_join(spnames[,c('plant_name_id','taxon_name')],by=c('accepted_plant_name_id'='plant_name_id'))
res.fin=b%>%left_join(trenames,by=c("Name_submitted"="Taxon_name"))%>%select(tip,Name_submitted,Taxonomic_status,Accepted_name_id,Accepted_name)
colnames(sp.mat)=colnames(res.fin)
spname.cor=rbind(sp.mat,res.fin)
save(spname.cor,file="data analysis/Phylo.tipnames.POWO.Rdata")

RangeFragments<-get(load("data analysis/NativeRangeFragments.plants.US.RData"))
splist=unique(names(RangeFragments))
spname.cor2=spname.cor%>%.[order(.[,"Taxonomic_status"]),]%>%.[!duplicated(.[,'Accepted_name']),]
tipname=spname.cor2%>%filter(Accepted_name%in%splist)
tre=keep.tip(tre0,tip=tre0$tip[tre0$tip%in%tipname$tip])
tre$tip.label=tipname[match(tre$tip.label,tipname$tip),"Accepted_name"]
save(tre,file="data analysis/Phylo.us.Rdata")
save(tre,file="data analysis/Phylo.us.park.Rdata")
## caculate phylo similarity
## phylogenetic beta diversity function
source("data analysis/read_Phylogenetic_Tree.R")
phyBeta <- function(comm, tree, ancestors=NULL, method=c("sim", "sim2")) {
	if (length(method) > 1) method <- method[1]

	comm1 <- comm[[1]]
	comm2 <- comm[[2]]

	if (is.null(ancestors)) ancs <- ancestor(tree, tip=tree$tip.label, TRUE)
	else ancs <- ancestors

	a1 <- unique(unlist(ancs[tree$tip.label %in% comm1]))
	a2 <- unique(unlist(ancs[tree$tip.label %in% comm2]))
	shared <- a1[a1 %in% a2]
	uniq.1 <- a1[!a1 %in% a2]
	uniq.2 <- a2[!a2 %in% a1]
	
	aa <- sum(tree$edge.length[which(tree$edge[,2] %in% shared)])
	bb <- sum(tree$edge.length[which(tree$edge[,2] %in% uniq.1)])
	cc <- sum(tree$edge.length[which(tree$edge[,2] %in% uniq.2)])

	if (method == "sim") result <- 1 - aa/(aa + min(bb,cc))
	if (method == "sim2") result <- list(i=1 - aa/(aa + bb),j=1 - aa/(aa + cc))
	return(result)
	}

library(ape);library(dplyr)	
#load("data analysis/Phylo.tipnames.POWO.Rdata")
tree0=get(load("data analysis/Phylo.us.Rdata"))
rangeData=get(load("data analysis/NativeRangeFragments.plants.ENA.clean.RData"))%>%filter(Tip_Label%in%tree0$tip)
tree=keep.tip(tree0,tip=unique(rangeData$Tip_Label))

d=rangeData[,1:2];d$val=1;rownames(d)=NULL
d=d %>% tidyr::pivot_wider(names_from = Tip_Label, values_from = val,values_fill = 0) %>% tibble::column_to_rownames("domain.ID")
#phylobeta=betapart::phylo.beta.pair(d, tree, index.family="sorensen")
## distribution data
for (i in 1:ncol(d)) d[,i] <- as.logical(d[,i])
d <- as.matrix(d)
load("data analysis/CommSimilarity.ENA.phylo2.rda")
# pBeta.mat <- matrix(NA, nrow=nrow(d), ncol=nrow(d))
# rownames(pBeta.mat) <- colnames(pBeta.mat) <- rownames(d)
ancs <- ancestor(tree=tree, tip=tree$tip.label, include.tip = TRUE)

ii <- 1;
for (i in 1:(nrow(d)-1)) {
	for (j in (i+1):nrow(d)) {
		comm1 <- colnames(d)[which(as.logical(d[i,]))]
		comm2 <- colnames(d)[which(as.logical(d[j,]))]
		comm <- list(comm1, comm2);names(comm)=paste(i,j,sep="_")
		pBeta.mat[i,j] <- phyBeta(comm,tree,ancs,method="sim")#simpson beta
		ii <- ii + 1
		}
	print(i)
	}
pBeta.mat[is.na(pBeta.mat)]=0
psim.t=pBeta.mat+t(pBeta.mat)
diag(psim.t)=NA
save(psim.t, file="data analysis/CommSimilarity.ENA.phylo.rda")

ii <- 1;
for (i in 1:(nrow(d)-1)) {
	for (j in (i+1):nrow(d)) {
		comm1 <- colnames(d)[which(as.logical(d[i,]))]
		comm2 <- colnames(d)[which(as.logical(d[j,]))]
		comm <- list(comm1, comm2);names(comm)=paste(i,j,sep="_")
		all.p <- phyBeta(comm,tree,ancs,method="sim2")# modified simpson like taxonmic beta
		pBeta.mat[i,j]=all.p[[1]]
		pBeta.mat[j,i]=all.p[[2]]
		ii <- ii + 1
		}
	print(i)
	}
save(pBeta.mat, file="data analysis/CommSimilarity.ENA.phylo2.rda")

## functional similarity
## Measure functional beta
library(dplyr)	
load("data analysis/traits.ENA.RData")
source("data analysis/MVNH_functions.R") ## ref:A unifying framework for quantifying and comparing ndimensional hypervolumes https://github.com/lvmuyang/MVNH
#A simple R package for quantifying and comparing multivariate normal hypervolumes and the partitioned components. 
# MVNH_det() to calculate niche volue and MVNH_dissimilarity() to calculate niche overlap. The dissimilarity of two hypervolumes with the Bhattacharyya distance ranges from 0 to unlimited
rangeData=get(load("data analysis/NativeRangeFragments.plants.ENA.clean.RData"))%>%filter(Tip_Label%in%rownames(trait_mat3))
d=rangeData[,1:2];d$val=1;rownames(d)=NULL
d=d %>% tidyr::pivot_wider(names_from = Tip_Label, values_from = val,values_fill = 0) %>% tibble::column_to_rownames("domain.ID")
#phylobeta=betapart::phylo.beta.pair(d, tree, index.family="sorensen")
## distribution data
for (i in 1:ncol(d)) d[,i] <- as.logical(d[,i])
d <- as.matrix(d)
Beta.turnover <- matrix(NA, nrow=nrow(d), ncol=nrow(d));rownames(Beta.turnover) <- colnames(Beta.turnover) <- rownames(d)
Beta.nest <- matrix(NA, nrow=nrow(d), ncol=nrow(d));rownames(Beta.nest) <- colnames(Beta.nest) <- rownames(d)

ii <- 1;
for (i in 1:(nrow(d)-1)) {
	for (j in (i+1):nrow(d)) {
		comm1 <- colnames(d)[which(as.logical(d[i,]))]
		comm2 <- colnames(d)[which(as.logical(d[j,]))]
		db1=subset(trait_mat3,rownames(trait_mat3)%in%comm1)
		db2=subset(trait_mat3,rownames(trait_mat3)%in%comm2)
		re=MVNH_dissimilarity(db1,db2,var.names = colnames(trait_mat3))
		Beta.turnover[i,j] <- re$M[1]#Mahalanobis distance between the two niche centroids
		Beta.nest[i,j] <- re$D[1]#determinant ratio component measuring the difference of the niche volumes		
		ii <- ii + 1
		}
	print(i)
	}

mf=function(pBeta.mat){
	pBeta.mat[is.na(pBeta.mat)]=0
	psim.t=pBeta.mat+t(pBeta.mat)
	diag(psim.t)=NA
	return(psim.t)
}
Beta.turnover=mf(Beta.turnover);Beta.nest=mf(Beta.nest)
Beta.traits=Beta.turnover+Beta.nest
save(Beta.turnover, file="data analysis/CommSimilarity.ENA.trait.turnover.rda")
save(Beta.traits, file="data analysis/CommSimilarity.ENA.trait.rda")

## Measure climate, population, urban and soil distance
## step 1 --- intergrate all the envir data -----
library(raster);library(dplyr);library(sf);
load("data analysis/Domain.ena.rda")
ENA.name=c("ME","NH","VT","NY","MA","RI","CT","NJ","PA","DE","MD","MI","OH","IN","IL","WI","WV","VA","NC","TN","KY","SC","GA","AL","MS","FL")
ENA <- read_sf('data analysis/cb_2018_us_state_20m/cb_2018_us_state_20m.shp')%>%filter(STUSPS%in%ENA.name)%>%as(., "Spatial")
## population ---
toload=list.files("data analysis/population/")%>%.[grep(".tif",.)]
pops <- raster::stack(paste0("data analysis/population/",toload))
pops.ena=ENA%>%raster::mask(pops,.)%>%aggregate(fact=4)%>%raster::crop(extent(ENA))
save(pops.ena,file="data analysis/population/pops.ena.rda")
load("data analysis/population/pops.ena.rda")
domain=domain.ena%>%select(x,y,domain.ID,row,col)
domain$x=round(domain$x,3);domain$y=round(domain$y,3)
pop.domain=raster::as.data.frame(pops.ena,xy=TRUE)%>%na.omit()
pop.domain$x=round(pop.domain$x,3);pop.domain$y=round(pop.domain$y,3)
colnames(pop.domain)[-c(1:2)]=paste("pop",seq(2000,2020,5),sep="")
pop.domain=pop.domain%>%left_join(domain,by=c("x","y"))%>%na.omit()
pop.domain[pop.domain<0.01]=NA
#imputation on grids with 0 population (74 out of 8490) using neighbour cells 
miss=pop.domain%>%filter(is.na(pop2000))
imp=function(i,pop.domain,vars){
	cells.row=(miss[i,]$row-5):(miss[i,]$row+5)
	cells.col=(miss[i,]$col-5):(miss[i,]$col+5)
	impute=pop.domain%>%filter(row%in%cells.row&col%in%cells.col)%>%dplyr::select(contains(vars))%>%as.matrix()%>%robustbase::colMedians(na.rm =TRUE)
	return(impute)
}
miss.imp=c()
for (i in paste("pop",seq(2000,2020,5),sep="")){
miss.imp.t=do.call(rbind,lapply(1:nrow(miss),imp,pop.domain,i))
miss.imp=cbind(miss.imp,miss.imp.t)
}
miss.imp=cbind(miss[,c("x","y")],miss.imp,miss[,c("domain.ID","row","col")])
pop.domain=pop.domain%>%filter(!is.na(pop2000))%>%rbind(miss.imp)%>%.[order(.[,"domain.ID"]),]%>%select(domain.ID,starts_with("pop"))

## urban ---
toload=list.files("data analysis/Urbanization/")%>%.[grep(".tif",.)]
urban <- raster::stack(paste0("data analysis/Urbanization/",toload))
urban.ena=ENA%>%raster::mask(urban,.)%>%raster::crop(extent(ENA))
load("data analysis/Urbanization/urban.ena.rda")
urban.domain=raster::as.data.frame(urban.ena,xy=TRUE)%>%na.omit()
colnames(urban.domain)[-c(1:2)]=c(paste("built",c(1975,1990,2000,2014),sep=""),paste("smod",c(1975,1990,2000,2015),sep=""))
dfr <- rasterFromXYZ(domain.ena[,c("x","y","domain.ID")],res=res(urban.ena),crs=crs(urban.ena))%>%stars::st_as_stars()%>%st_as_sf()
urban.domain$urban.ID=1:nrow(urban.domain)
urban.dfr=rasterFromXYZ(urban.domain[,c("x","y","urban.ID")],res=res(urban.ena),crs=crs(urban.ena))%>%stars::st_as_stars()%>%st_as_sf()
#intersec= st_intersection(dfr,urban.dfr)
intersec=st_intersects(dfr,urban.dfr)
mf=function(i,urban.domain) urban.domain%>%filter(urban.ID%in%i)%>%colMeans()
urban.domain2=do.call(rbind,lapply(intersec,mf,urban.domain))
urban.domain2=cbind(domain.ena[,c("x","y","domain.ID")],urban.domain2[,c(paste("built",c(1975,1990,2000,2014),sep=""),paste("smod",c(1975,1990,2000,2015),sep=""))])
#p=rasterFromXYZ(urban.domain2[,c("x","y","smod2015")],res=res(pops.ena),crs=crs(pops.ena))
urban.pop=urban.domain2%>%left_join(pop.domain,by="domain.ID")

## road density ---
ele.us <- raster("data analysis/wc2.1_10m_bio/wc2.1_10m_elev.tif")
road=raster("data analysis/road/grip4_total_dens_m_km2.asc")
road.ena=ENA%>%raster::mask(road,.)%>%raster::crop(extent(ENA))%>%stars::st_as_stars()%>%st_as_sf()
dfr <- rasterFromXYZ(domain.ena[,c("x","y","domain.ID")],res=res(ele.us),crs=crs(ele.us))%>%stars::st_as_stars()%>%st_as_sf()
intersec=st_intersects(dfr,road.ena)
road.domain=data.frame(road.ID=1:nrow(road.ena),road.density=as.data.frame(road.ena)[,1])
mf=function(i,road.domain) if (length(i)>0) road.domain%>%filter(road.ID%in%i)%>%as.matrix()%>%colMeans(na.rm =TRUE) else NA
road.domain2=do.call(rbind,lapply(intersec,mf,road.domain))
road.domain3=data.frame(domain.ena[,c("x","y","domain.ID","row","col")],road.density=road.domain2[,2])
miss=road.domain3%>%filter(is.na(road.density))
imp=function(i,pop.domain,vars){
	cells.row=(miss[i,]$row-1):(miss[i,]$row+1)
	cells.col=(miss[i,]$col-1):(miss[i,]$col+1)
	impute=pop.domain%>%filter(row%in%cells.row&col%in%cells.col)%>%dplyr::select(contains(vars))%>%as.matrix()%>%robustbase::colMedians(na.rm =TRUE)
	return(impute)
}
miss.imp=do.call(rbind,lapply(1:nrow(miss),imp,road.domain3,"road.density"))
miss.imp=cbind(miss[,-6],miss.imp)
road.domain3=road.domain3%>%filter(!is.na(road.density))%>%rbind(miss.imp)%>%.[order(.[,"domain.ID"]),]%>%select(domain.ID,"road.density")

## airport and trainstation ---
airport=read.csv("data analysis/road/airport.us.csv")
spdf <- SpatialPoints(coords = airport[,c("lon","lat")], proj4string = crs(ele.us))%>%raster::crop(extent(ENA))%>%st_as_sf()
intersec=st_intersects(dfr,spdf)
trans.cout=data.frame(domain.ena[,c("x","y","domain.ID")],airport=do.call(rbind,lapply(intersec,length)))%>%left_join(road.domain3,by="domain.ID")

#p=rasterFromXYZ(trans.cout[,c("x","y","road.density")],res=res(ele.us),crs=crs(ele.us));plot(p)

## soil ---
ele.us <- raster("data analysis/wc2.1_10m_bio/wc2.1_10m_elev.tif")
toload=list.files("data analysis/Soil conditions/")
soil0=raster::stack(paste0("data analysis/Soil conditions/",toload))
soil.ena=projectRaster(soil0,crs = crs(ele.us))%>%raster::crop(extent(ENA))
save(soil.ena,file="data analysis/Soil conditions/soil.ena.rda")
load("data analysis/Soil conditions/soil.ena.rda")
soil.domain=raster::as.data.frame(soil.ena,xy=TRUE)%>%na.omit()
colnames(soil.domain)[-c(1:2)]=gsub("cm_mean_5000","",colnames(soil.domain)[-c(1:2)])

dfr <- rasterFromXYZ(domain.ena[,c("x","y","domain.ID")],res=res(ele.us),crs=crs(ele.us))%>%stars::st_as_stars()%>%st_as_sf()
soil.domain$soil.ID=1:nrow(soil.domain)
soil.dfr=rasterFromXYZ(soil.domain[,c("x","y","soil.ID")],res=res(soil.ena),crs=crs(soil.ena))%>%stars::st_as_stars()%>%st_as_sf()
intersec=st_intersects(dfr,soil.dfr)
mf=function(i,soil.domain) soil.domain%>%filter(soil.ID%in%i)%>%colMeans(na.rm=T)
soil.domain2=do.call(rbind,lapply(intersec,mf,soil.domain))
soil.domain2=cbind(domain.ena[,c("x","y","domain.ID","row","col")],soil.domain2[,-c(1,2,82)])

#imputation on grids with NA (43 out of 8490) using neighbour cells 
miss=soil.domain2%>%filter(is.na(bdod_0.5))
miss.imp=c()
for (i in colnames(soil.domain2)[-c(1:5)]){
miss.imp.t=do.call(rbind,lapply(1:nrow(miss),imp,soil.domain2,i))
miss.imp=cbind(miss.imp,miss.imp.t)
}
miss.imp=cbind(miss[,c("x","y")],miss.imp,miss[,c("domain.ID","row","col")])
soil.domain3=soil.domain2%>%filter(!is.na(bdod_0.5))%>%rbind(miss.imp)%>%.[order(.[,"domain.ID"]),]%>%select(-row,-col)
#p=rasterFromXYZ(soil.domain3[,c("x","y","bdod_0.5")],res=res(ele.us),crs=crs(ele.us))
envir.ena=list(urban.pop,soil.domain3,trans.cout,domain.ena) %>% purrr::reduce(left_join,by=c("domain.ID","x","y"))
save(envir.ena,file="data analysis/envir.ena.rda")

## step 2 --- caculate envir dist -----
load("data analysis/envir.ena.rda")
library(dplyr)
clim=envir.ena%>%dplyr::select(ele,starts_with("bio"))%>%ade4::dudi.pca(center = TRUE,scale = T,scannf = FALSE,nf = 5)
human=data.frame(pop=rowMeans(envir.ena%>%dplyr::select(starts_with("pop"))),
	built=rowMeans(envir.ena%>%dplyr::select(starts_with("built"))),
	urban=rowMeans(envir.ena%>%dplyr::select(starts_with("smod"))),
	envir.ena%>%dplyr::select(airport,road.density))
#human=envir.ena%>%dplyr::select(starts_with(c("pop","smod","built")),airport,road.density)
human[human==0]=0.000001
human=human%>%#log()%>%
	ade4::dudi.pca(center = TRUE,scale = T,scannf = FALSE,nf = 4)
soil=envir.ena%>%dplyr::select(contains("_"))%>%ade4::dudi.pca(center = TRUE,scale = T,scannf = FALSE,nf = 5)
clim.soil=envir.ena%>%dplyr::select(ele,starts_with("bio"),contains("_"))%>%ade4::dudi.pca(center = TRUE,scale = T,scannf = FALSE,nf = 5)
	
clim.d=dist(clim$li)%>%as.matrix()%>%scales::rescale();diag(clim.d)=NA
soil.d=dist(soil$li)%>%as.matrix()%>%scales::rescale();diag(soil.d)=NA
clim.soil.d=dist(clim.soil$li)%>%as.matrix()%>%scales::rescale();diag(clim.soil.d)=NA
human.d=dist(human$li)%>%as.matrix()%>%scales::rescale();diag(human.d)=NA
clim.soil.pop.d=dist(cbind(clim.soil$li,human$li))%>%as.matrix()%>%scales::rescale();diag(clim.soil.pop.d)=NA
clim.human.d=dist(cbind(clim$li,human$li))%>%as.matrix()%>%scales::rescale();diag(clim.human.d)=NA
#p=rasterFromXYZ(cbind(envir.ena[,c("x","y")],clim.soil.d[,1]),res=res(ele.us),crs=crs(ele.us));plot(p)
#envir.d=list(clim.d,soil.d,clim.soil.d,human.d,clim.soil.pop.d,clim.human.d);names(envir.d)=c("clim","soil","clim.soil","human","clim.soil.pop","clim.human")
envir.d=list(clim.soil.pop.d,clim.human.d);names(envir.d)=c("clim.soil.pop","clim.human")
rm(envir.ena,clim,human,soil,clim.d,human.d,soil.d,clim.soil,clim.soil.d,clim.soil.pop.d)
gc()
#################################################### 
############## Simulations   ####################### 
#################################################### 
#ref:https://doi.org/10.6084/m9.figshare.13270406
library(data.table);library(dplyr)
load("data analysis/Domain.ena.rda")
#SimFiles=list.files("data analysis/")
#toload<-SimFiles[grep("CommSimilarity.ENA.",SimFiles)]
#RangeFragments<-get(load("data analysis/AlienRangeFragments.plants.ENA.RData"))
load("data analysis/envir.ena.rda")

RangeFragments<-get(load("data analysis/AlienRangeFragments.alltaxa.ENA.RData"))
fragmentTable1=data.frame(fragBinomial=unique(names(RangeFragments)),fragSpecies=1:length(unique(names(RangeFragments))))
fragmentTable2=data.frame(fragNum=1:length(RangeFragments),fragSize=unlist(lapply(RangeFragments,length)),fragBinomial=names(RangeFragments))
fragmentTable=fragmentTable2%>%left_join(fragmentTable1,by="fragBinomial")
rm(fragmentTable1,fragmentTable2)

mf=function(i,RangeFragments) {return(data.frame(Tip_Label=names(RangeFragments)[i],domain.ID=RangeFragments[[i]],fragNum=i))}
IntroductionsFiltered=do.call(rbind,lapply(1:length(RangeFragments),mf,RangeFragments))

get.sim=function(REP,fragmentTable,RangeFragments,EnvironmentData,IntroductionsFiltered,adjCells.ena,CommSimilarity,mu,mat.name,inputProbSurface=TRUE){
	source("data analysis/Functions.R")
	#source("Code/2.Functions/Functions.R")
	simResults<-SimulateSpread(Range="Alien",fragmentTable=fragmentTable,RangeFragments=RangeFragments,
									 EnvironmentData=EnvironmentData,IntroductionsFiltered=IntroductionsFiltered,adjCells=adjCells.ena,
									 SpeciesOptimum=NULL,CommSimilarity=CommSimilarity,plotIt=FALSE,envVar=c(),mu=mu,
									 seed="Random",niche="ClimateMatching",inputProbSurface=inputProbSurface)
	names(simResults)=names(RangeFragments)		
    outfile<-paste0("/blue/matthewthomas1/yunpeng.liu/data analysis/simulation.ena2/",mat.name,"/ER_",mat.name,"_",mu,"_",REP,".rda")
	save(simResults,file=outfile)
}

library(parallel)
no_cores <- 40#detectCores() - 1
### simulate environmental resistance model with different mu values
muval<-sort(c(0.5,1.5,3,5,10,30,100),decreasing =TRUE)
#muval=c(10,5,0.5,30)
#mat.list=c("human","clim.human","clim.soil.pop","null","clim","sp","phylo2","soil","clim.soil","trait","trait.turnover")	
mat.list=c("clim.human","clim.soil.pop","null")

gc() 
for (mat in 1:length(mat.list)){
	mat.name=mat.list[mat]
	print(mat.name)
	if (mat.name=="null"){
		mycl <- makePSOCKcluster(no_cores);
		parSapply(cl=mycl,X=101:200,get.sim,fragmentTable,RangeFragments,EnvironmentData=domain.ena,IntroductionsFiltered,adjCells.ena,
			CommSimilarity=NULL,mu=NA,mat.name="null",inputProbSurface=FALSE)
		stopCluster(mycl)	
	}else{
		if (mat.name%in%names(envir.d)){		
			simmat=envir.d[[mat.name]]
		}else{
			simmat=get(load(paste0("data analysis/CommSimilarity.ENA.",mat.name,".rda")))
		}
		if (grepl("sp",mat.name)+grepl("sp2",mat.name)==0) simmat=1-simmat
		if (max(simmat,na.rm=T)>1) simmat=scales::rescale(simmat)
		# nstarts=ifelse(mat.name=="clim",7,1)	
		# if (nstarts==6) {toload=list.files("data analysis/simulation.ena/clim")%>%.[grep("ER_clim_1.5_",.)]%>%gsub("ER_clim_1.5_","",.)%>%gsub(".rda","",.)%>%as.numeric();toload=c(1:100)[-toload]
		# }else{toload=1:100}
	for(i in 1:length(muval)){
		mycl <- makePSOCKcluster(no_cores);		
		parSapply(cl=mycl,X=101:200,get.sim,fragmentTable,RangeFragments,EnvironmentData=domain.ena,IntroductionsFiltered,adjCells.ena,
			CommSimilarity=simmat,mu=muval[i],mat.name)
	  # simResults<- SimulateSpread(Range=Range,fragmentTable=fragmentTable,RangeFragments=RangeFragments,
									# EnvironmentData=EnvironmentData,IntroductionsFiltered=IntroductionsFiltered,adjCells=adjCells.fl,
									# SpeciesOptimum=SpeciesOptimum,CommSimilarity=CommSimilarity,plotIt=FALSE,envVar=c(),mu=muval[i],
									# seed=seed,niche="ClimateMatching",inputProbSurface=TRUE)
	   print(muval[i])
	   stopCluster(mycl)
	   gc()  
	  }	
	}	
}
		
# plot observed range
	windows()
	library(scales)
	tmp=get(load("data analysis/simulation/Alien_Random_ClimateResistance_ele_30_17.rda"))
	currspeciesname=unique(names(RangeFragments))[45]
	plot(EnvironmentData$x,EnvironmentData$y,cex=0.5)
	focalDom<-EnvironmentData%>%filter(domain.ID%in%unlist(subset(RangeFragments,names(RangeFragments) %in% currspeciesname)))
    points(focalDom$x,focalDom$y,col="green",pch=16,cex=0.5)  
	focalDom2<-EnvironmentData%>%filter(domain.ID%in%unlist(subset(tmp,names(tmp) %in% currspeciesname)))
    points(focalDom2$x,focalDom2$y,col=alpha("red", 0.5),pch=16,cex=0.5)	
########################################################################################
# for each model and alien range fragment calculate the probability that each cell is invaded across the 100 replicate runs
########################################################################################
library(dplyr)
#mat.list=c("clim","clim.pca","phylo","phylo2","sp","sp2","trait","trait.turnover","null")	
mat.list=c("null","clim","soil","clim.soil","human","clim.soil.pop","sp","phylo2","trait","trait.turnover")	
RangeFragments<-get(load("data analysis/AlienRangeFragments.alltaxa.ENA.RData"))
fragmentTableA=data.frame(fragBinomial=unique(names(RangeFragments)),fragSpecies=1:length(unique(names(RangeFragments))))
fragmentTableB=data.frame(fragNum=1:length(RangeFragments),fragSize=unlist(lapply(RangeFragments,length)),fragBinomial=names(RangeFragments))
AlienfragmentTable=fragmentTableB%>%left_join(fragmentTableA,by="fragBinomial")
get.prob=function(i,toload,mat){
		 get.simdis=function(j,i){
			#re=get(load(paste0("/blue/matthewthomas1/yunpeng.liu/data analysis/simulation.ena2/",mat,"/",j)))[[i]]  
			re=get(load(paste0("data analysis/simulation.ena2/",mat,"/",j)))[[i]]  
		   return(re)
		  }
		  simOutput.t<-do.call(c,lapply(toload,get.simdis,i))
		  simOutput.freq=table(simOutput.t)/length(toload)
		  return(simOutput.freq)
		  rm(re,simOutput.t);gc()
	  }
intersectRanges<-function(x,y){return(length(intersect(x,y)))}	  
library(parallel)
no_cores <- detectCores() - 1
muval<-sort(c(0.5,1.5,3,5,10,30,100),decreasing =TRUE)
nFrag=nrow(AlienfragmentTable);nReps=100;
overlapArrayMuVals=list();
for (mat in mat.list){
	if (mat=="null"){
		overlapArrayMuVals[[mat]]=array(dim=c(nFrag,1,nReps))
		toload<-list.files(paste0("/blue/matthewthomas1/yunpeng.liu/data analysis/simulation.ena/",mat,"/"))
		toload=gtools::mixedsort(sort(toload))
		#tmp=get(load(paste0("data analysis/simulation.ena/",mat,"/",toload[1])))
		#AlienfragmentTable=AlienfragmentTable.all%>%filter(fragBinomial%in%names(tmp))
		print(c(mat,length(toload)))
		
		mycl <- makePSOCKcluster(no_cores);
		InvasionProb<-parLapply(cl=mycl,X=1:nFrag,get.prob,toload,mat)  
		names(InvasionProb)<-AlienfragmentTable$fragNum	
		outfile<-paste("/blue/matthewthomas1/yunpeng.liu/data analysis/InvasionProbabilites.ena/Alien_Random",mat,"rda",sep=".")
		save(InvasionProb,file=outfile)	 
		stopCluster(mycl);
		rm(InvasionProb);gc()
		
		#Fragments=subset(RangeFragments,names(RangeFragments)%in%names(tmp))	
		for(i in 1:length(toload)){
			load(paste0("/blue/matthewthomas1/yunpeng.liu/data analysis/simulation.ena/",mat,"/",toload[i]))
			overlap<-mapply(intersectRanges,RangeFragments,simResults)/AlienfragmentTable$fragSize
			overlapArrayMuVals[[mat]][,1,i]<-overlap
		}				
	}else{
		overlapArrayMuVals[[mat]]=array(dim=c(nFrag,length(muval),nReps))		
		for(x in 1:length(muval)){
			#SimFiles<-list.files(paste0("/blue/matthewthomas1/yunpeng.liu/data analysis/simulation.ena2/",mat,"/"))	
			SimFiles<-list.files(paste0("data analysis/simulation.ena2/",mat,"/"))	
			toload<-SimFiles[grep(paste0("ER_",mat,"_",muval[x],"_"),SimFiles)]
			toload=gtools::mixedsort(sort(toload))
			print(c(mat,muval[x],length(toload)))
			
			mycl <- makePSOCKcluster(no_cores);		
			InvasionProb<-parLapply(cl=mycl,X=1:nFrag,get.prob,toload,mat)  
			names(InvasionProb)<-AlienfragmentTable$fragNum	
			#outfile<-paste("/blue/matthewthomas1/yunpeng.liu/data analysis/InvasionProbabilites.ena/Alien_Random",mat,muval[x],"rda",sep=".")
			outfile<-paste("data analysis/InvasionProbabilites.ena/Alien_Random",mat,muval[x],"rda",sep=".")
			save(InvasionProb,file=outfile)	  	
			stopCluster(mycl);
			rm(InvasionProb);gc()
			
			for(i in 1:length(toload)){
			  #load(paste0("/blue/matthewthomas1/yunpeng.liu/data analysis/simulation.ena/",mat,"/",toload[i]))
			  load(paste0("data analysis/simulation.ena2/",mat,"/",toload[i]))
			  overlap<-mapply(intersectRanges,RangeFragments,simResults)/AlienfragmentTable$fragSize
			  overlapArrayMuVals[[mat]][,x,i]<-overlap		 
			}	
		  }
		}
	  #rm(overlapArrayMuVals)
}
save(overlapArrayMuVals,file=paste0("/blue/matthewthomas1/yunpeng.liu/data analysis/InvasionProbabilites.ena/ModelOverlapScores.rda"))
save(overlapArrayMuVals,file=paste0("data analysis/InvasionProbabilites.ena/ModelOverlapScores.human.rda"))	  
########################################################################################
# calculate species level accuracy (i.e. overlap) scores
########################################################################################
library(dplyr);library(matrixStats)
load("/blue/matthewthomas1/yunpeng.liu/data analysis/InvasionProbabilites.ena/ModelOverlapScores.rda")
load("data analysis/AlienRangeFragments.alltaxa.ENA.RData")

# calculate species accuracy scores for each model
source("data analysis/Functions.R")
muval<-sort(c(0.5,1.5,3,5,10,30,100),decreasing =TRUE)
mat.list=c("null","clim","sp","soil","human","phylo2","trait","clim.soil","clim.soil.pop","trait.turnover")	
se=function(x) 1.96*sd(x)/sqrt(length(x))
fragmentTable1=data.frame(fragBinomial=unique(names(RangeFragments.ena.inv)),fragSpecies=1:length(unique(names(RangeFragments.ena.inv))))
fragmentTable2=data.frame(fragNum=1:length(RangeFragments.ena.inv),fragSize=unlist(lapply(RangeFragments.ena.inv,length)),fragBinomial=names(RangeFragments.ena.inv))
AlienfragmentTable=fragmentTable2%>%left_join(fragmentTable1,by="fragBinomial")
	
speciesAlien=unique(AlienfragmentTable$fragBinomial);nSpeciesAlien=length(speciesAlien)
OScoresMuvalsAlienMedian<-OScoresMuvalsAlienSE<-OScoresMuvalsAlienMean<-list()
AlienfragmentTable$FocalModel<-NA
	
for (mat in mat.list){
	print(mat)
	if (mat=="null") {
		OScoresMuvalsAlienMean[[mat]]<-matrix(ncol=1,nrow=nSpeciesAlien)
	}else {
		OScoresMuvalsAlienMean[[mat]]<-matrix(ncol=length(muval),nrow=nSpeciesAlien);colnames(OScoresMuvalsAlienMean[[mat]])=muval
	}
	rownames(OScoresMuvalsAlienMean[[mat]])<-speciesAlien;
	OScoresMuvalsAlienMedian[[mat]]<-OScoresMuvalsAlienSE[[mat]]<-OScoresMuvalsAlienMean[[mat]]
	if (mat=="null") {
		AlienfragmentTable$FocalModel<-1
		scoreMatrix<-calculateSpeciesOverlapScores(FragTable=AlienfragmentTable,overlapArray=overlapArrayMuVals[[mat]])
		OScoresMuvalsAlienMean[[mat]][,1]<-rowMeans(scoreMatrix)
		OScoresMuvalsAlienMedian[[mat]][,1]<-rowQuantiles(scoreMatrix,probs=c(0.5))
		OScoresMuvalsAlienSE[[mat]][,1]<-apply(scoreMatrix,1,se)	  
	}else {		
		for(x in 1:length(muval)){
			AlienfragmentTable$FocalModel<-x
			scoreMatrix<-calculateSpeciesOverlapScores(FragTable=AlienfragmentTable,overlapArray=overlapArrayMuVals[[mat]])
			OScoresMuvalsAlienMean[[mat]][,x]<-rowMeans(scoreMatrix)
			OScoresMuvalsAlienMedian[[mat]][,x]<-rowQuantiles(scoreMatrix,probs=c(0.5))
			OScoresMuvalsAlienSE[[mat]][,x]<-apply(scoreMatrix,1,se)
			  #print(x)
		}
	}	
}
save(OScoresMuvalsAlienMean,OScoresMuvalsAlienMedian,OScoresMuvalsAlienSE,
	file=paste0("/blue/matthewthomas1/yunpeng.liu/data analysis/InvasionProbabilites.ena/SpeciesModelOverlapScores.rda"))

########################################################################################
## model comparison
########################################################################################
library(ggpubr);library(dplyr);library(tidyr)
load("data analysis/InvasionProbabilites.ena/SpeciesModelOverlapScores.rda")
load("data analysis/AlienRangeFragments.alltaxa.ENA.RData")
#mat.list=c("null","clim","clim.pca","phylo","phylo2","sp","sp2","trait","trait.turnover")	
mat.list=c("null","clim","soil","human","clim.soil","clim.soil.pop","sp","phylo2","trait","trait.turnover")	
muval<-sort(c(0.5,1.5,3,5,10,30,100),decreasing =TRUE)
splist=read.csv("data analysis/US-RIIS.csv")%>%select(kingdom,phylum,class,scientificName,habitat)%>%filter(scientificName%in%names(RangeFragments.ena.inv))%>%
	filter(!phylum%in%c("Rhodophyta","Chlorophyta","Charophyta"))
	#filter(!class%in%c("Teleostei","Ascidiacea","Malacostraca","Maxillopoda","Ostracoda","Branchiopoda"))
# aquatic=rbind(splist[grep("Aquatic",splist$habitat),],splist[grep("Marine",splist$habitat),])
# splist=splist%>%filter(!scientificName%in%aquatic$scientificName)
taxa=c("Chordata","Arthropoda","Other animal","Tracheophyta")#other animal incl. Mollusca,Annelida,Nematoda, etc.
splist%>%group_by(kingdom,phylum)%>%tally()

##step 1: statistic
re=c()
for (taxa.i in taxa){
	sp=splist[splist$phylum==taxa.i,"scientificName"]
	if (taxa.i%in%"Other animal") sp=splist[splist$kingdom=="Animalia"&(!splist$phylum%in%c("Chordata","Arthropoda")),"scientificName"]
	
	null.mean=OScoresMuvalsAlienMean[["null"]]%>%subset(rownames(.)%in%sp)
	null.se=OScoresMuvalsAlienSE[["null"]]%>%subset(rownames(.)%in%sp)
	for (mat in mat.list[-1]){
		for (mu in 1:length(muval)){
		er.mean=as.matrix(OScoresMuvalsAlienMean[[mat]][,mu])%>%subset(rownames(.)%in%sp);colnames(er.mean)=muval[mu]
		er.se=as.matrix(OScoresMuvalsAlienSE[[mat]][,mu])%>%subset(rownames(.)%in%sp);colnames(er.se)=muval[mu]
		
		#how many species does ER outperform random dispersal?
		diffScore<-er.mean-null.mean
		a=length(which(diffScore>0))/length(diffScore)

		# how many species can we reject the null model?
		diffScore<-er.mean-(null.mean+null.se)
		b=length(diffScore[diffScore>0])/length(diffScore)

		sta=t.test(er.mean, null.mean,paired=T)
		tmp=data.frame(taxa=taxa.i,model=mat,mu=muval[mu],a,b,mean.s=mean(er.mean),diff.with.null=sta$est,p=sta$p.value)
		re=rbind(re,tmp)
		}
	}
}
re2=re%>%group_by(taxa,model)%>%summarize(a=max(a),b=max(b),mean.s=max(mean.s),diff.with.null=max(diff.with.null),p=p[which(diff.with.null==max(diff.with.null))])%>%as.data.frame()
re.mu=re%>%group_by(model,taxa)%>%reframe(a=mu[which(a==max(a))],b=mu[which(b==max(b))],mean.s=mu[which(mean.s==max(mean.s))],diff.with.null=mu[which(diff.with.null==max(diff.with.null))])%>%as.data.frame()

## step2: compare different muvals and models
range.size=data.frame(fragSize=unlist(lapply(RangeFragments.ena.inv,length)),Species=names(RangeFragments.ena.inv))%>%group_by(Species)%>%summarize(rangesize=sum(fragSize))
score=c()
for (mat in mat.list){
	if (mat!="null"){
		score.mean=OScoresMuvalsAlienMean[[mat]]%>%as.data.frame()%>%mutate(Species=rownames(.),model=mat)%>%pivot_longer(!(Species|model),names_to = "mu", values_to = "score.mean")
		score.median=OScoresMuvalsAlienMedian[[mat]]%>%as.data.frame()%>%mutate(Species=rownames(.))%>%pivot_longer(!Species,names_to = "mu", values_to = "score.median")
		score.se=OScoresMuvalsAlienSE[[mat]]%>%as.data.frame()%>%mutate(Species=rownames(.))%>%pivot_longer(!Species,names_to = "mu", values_to = "score.se")	
		score.t=list(score.mean,score.se,score.median) %>% purrr::reduce(left_join,by=c('Species','mu'))
	}else{
		score.t=cbind(OScoresMuvalsAlienMean[[mat]],OScoresMuvalsAlienSE[[mat]],OScoresMuvalsAlienMedian[[mat]]) %>% as.data.frame() 
		colnames(score.t)=c("score.mean","score.se","score.median")
		score.t=score.t%>%mutate(Species=rownames(.),model=mat,mu=NA)%>%select(Species:mu,score.mean:score.median)
	}	
	score.t=score.t%>%left_join(range.size,by="Species")
	score=rbind(score,score.t)
}

dat=c()
for (taxa.i in taxa){
	sp=splist[splist$phylum==taxa.i,"scientificName"]
	if (taxa.i%in%"Other animal") sp=splist[splist$kingdom=="Animalia"&(!splist$phylum%in%c("Chordata","Arthropoda")),"scientificName"]
	score2=c()
	for (mat in mat.list){
		if (mat=="null") tmp=score%>%filter(Species%in%sp&is.na(mu)) else{
			bestmu=re.mu[re.mu$model==mat&re.mu$taxa==taxa.i,"a"]
			tmp=score%>%filter(Species%in%sp&model==mat&mu==bestmu)}
		score2=rbind(score2,tmp)
	}
	score3=score2
	score3$score.mean=round(score3$score.mean,2)
	score3=score3%>%group_by(Species)%>%reframe(bestmodel=model[which(score.mean==max(score.mean))])
	# if multiple models are tied for first place for a species, downweight these so each species contributes equally
	nModelsTied<-score3%>%group_by(Species)%>%tally()
	Weights<-1/nModelsTied$n
	WeightsVec<-rep(Weights,nModelsTied$n)
	BestModelsTable<-cbind(score3,WeightsVec)%>%group_by(bestmodel)%>%summarize(spnum=sum(WeightsVec))
	BestModelsTable$bestmodel=factor(BestModelsTable$bestmodel,levels=mat.list)
	dat=rbind(dat,data.frame(taxa=taxa.i,BestModelsTable))
}

theme=theme_bw()+theme(axis.text = element_text(size=12,color='black'),
			axis.text.x = element_text(angle=30),
			axis.title = element_text(size=18,color='black'),
			legend.background=element_rect(fill='transparent'),
			legend.text=element_text(face="bold",size=15),
			legend.title=element_text(face="bold",size=18),
			strip.text=element_text(size=18,color='black',face="bold"),
		strip.background = element_rect(fill = 'white', colour = 'black', size = rel(2), linetype = 2))
dat$model=factor(dat$bestmodel,levels=mat.list)	
ggplot(dat, aes(x=bestmodel, y=spnum)) + 
		geom_bar(stat = "identity")+
		#annotate("text",x=7,y=0.95,size=8,label=taxa.i,color="red")+
		labs(x="Type of models",y="Number of invasive species")+   
		facet_wrap(~ taxa,scales="free_y")+theme
		
# #compare different models
	# p[[taxa.i]]=ggplot(score2,aes(x = rangesize,y=score.mean,fill=model))+
		# geom_ribbon(aes(x = rangesize,ymin=score.mean-score.se, ymax=score.mean+score.se,fill=model),alpha=0.3,show.legend=FALSE)+
		# geom_line(aes(x=rangesize, y=score.mean,colour=model),size=0.8,linetype=1,alpha=0.8)+
		# # scale_color_discrete(labels = c("Random dispersal", "Species similarity","Phylogenetic similarity","Trait distance","Trait turnover",
			# # "Climatic similarity","Soil similarity","Human activity","Cliamte and soil", "Climate, soil and human activity"))+
		# # geom_errorbar(aes(ymin=score.mean-score.se, ymax=score.mean+score.se,color=model),linewidth=1,width=150,show.legend=FALSE,alpha=0.6)+
		# # geom_point(size=3,shape=21,color="black",alpha=0.6) + 
		# #geom_smooth(aes(x = rangesize,y=score.mean,color=model),data=score2,method = "loess",size=1.5,show.legend=FALSE,se =TRUE,span=0.8,linetype=1,alpha=0.5)+
		# annotate("text",x=2000,y=0.95,size=8,label=taxa.i,color="red")+
		# labs(x="Range sizes of species",y="Model accuracy scores",colour="Model") +theme+theme(legend.position =c(0.4,0.8))		

############################################################
######### map predicted ranges for different ER model #########
############################################################
library(dplyr);
load("data analysis/domain.ena.rda")
mat.list=c("null","clim","soil","clim.soil","human","clim.soil.pop","sp","phylo2","trait","trait.turnover")
inv.prob.files=list.files("data analysis/InvasionProbabilites.ena/")%>%.[grep("Alien_Random",.)]
load("data analysis/AlienRangeFragments.alltaxa.ENA.RData")
AlienfragmentTable=data.frame(fragNum=1:length(RangeFragments.ena.inv),fragSize=unlist(lapply(RangeFragments.ena.inv,length)),fragBinomial=names(RangeFragments.ena.inv))

get.rich=function(mat,AlienfragmentTable){
	require(matrixStats);require(dplyr)
	mf=function(x)data.frame(domain.ID=names(x),Freq=as.numeric(x))
	if(mat%in%c("clim","trait")) mu=10
	if(mat%in%c("soil","clim.soil","clim.soil.pop")) mu=5
	if(mat%in%c("phylo2","sp","human")) mu=0.5
	if(mat%in%"trait.turnover") mu=30
	if (mat%in%"null") ER.model=get(load(paste("data analysis/InvasionProbabilites.ena/Alien_Random",mat,"rda",sep="."))) else {
		ER.model=get(load(paste("data analysis/InvasionProbabilites.ena/Alien_Random",mat,mu,"rda",sep=".")))}
	names(ER.model)=AlienfragmentTable$fragBinomial
	rich.ER=do.call(rbind,lapply(ER.model,mf))
	rich.ER[,1]=as.numeric(rich.ER[,1])
	rich.ER=rich.ER%>%group_by(domain.ID)%>%summarize(ER=sum(Freq))
	colnames(rich.ER)[2]=mat
	return(rich.ER)
}

rich.ER.all=list()
for (mat in mat.list) rich.ER.all[[mat]]=get.rich(mat,AlienfragmentTable);print(mat)

library(parallel)
no_cores=detectCores()-1
mycl <- makePSOCKcluster(no_cores);	
rich.ER.all=parLapply(cl=mycl,mat.list,get.rich,AlienfragmentTable)
stopCluster(mycl)
# for (taxa.i in taxa){
	# if(taxa.i=="Plantae"){
		# RangeFragments=subset(RangeFragments.ena.inv,names(RangeFragments.ena.inv)%in%splist[splist$kingdom=="Plantae","scientificName"])
	# }else{
		# if (taxa.i%in%c("Chordata","Arthropoda")) RangeFragments=subset(RangeFragments.ena.inv,names(RangeFragments.ena.inv)%in%splist[splist$phylum==taxa.i,"scientificName"])
		# if (taxa.i%in%"Other animal") RangeFragments=subset(RangeFragments.ena.inv,names(RangeFragments.ena.inv)%in%splist[splist$kingdom=="Animalia"&(!splist$phylum%in%c("Chordata","Arthropoda")),"scientificName"])
	# }
	# fragmentTable1=data.frame(fragBinomial=unique(names(RangeFragments)),fragSpecies=1:length(unique(names(RangeFragments))))
	# fragmentTable2=data.frame(fragNum=1:length(RangeFragments),fragSize=unlist(lapply(RangeFragments,length)),fragBinomial=names(RangeFragments))
	# AlienfragmentTable=fragmentTable2%>%left_join(fragmentTable1,by="fragBinomial")
	# AlienRichness<-table(unlist(RangeFragments))

	# inv.prob.files=list.files(paste0("data analysis/InvasionProbabilites.ena/",taxa.i,"/"))
	# rich.ER.all=lapply(mat.list,get.rich,AlienfragmentTable,inv.prob.files)
	# names(rich.ER.all)=mat.list

	
# }
AlienRichness<-table(unlist(RangeFragments))
domain=domain.ena[,c("domain.ID","x","y")]
domain$CurrentRichness<-NA
domain$CurrentRichness[match(as.numeric(names(AlienRichness)),domain$domain.ID)]<-as.numeric(AlienRichness)
domain=c(list(domain),rich.ER.all) %>% purrr::reduce(left_join,by="domain.ID")
ma=ceiling(max(apply(domain[,-c(1:3)],2,max,na.rm=T)))

#plot richness
library(ggpubr)
the=theme_bw()+theme(legend.background=element_rect(fill='transparent'),
		legend.text=element_text(size=12),
		legend.title=element_text(face="bold",size=10))
		#legend.position=c(0.75,0.2))
p1=ggplot() +  geom_tile(data = domain, aes(x = x, y = y, fill = CurrentRichness)) +
  coord_quickmap() +  labs(x="",y="",fill="Richness")+
  scale_fill_gradientn(colours = viridis::viridis(99),na.value = "lightgray",limits=c(0,ma))+the+
  annotate("text",x=-72,y=30,size=4,label="Alien \nspecies \nRichness",color="red")
  
p2=ggplot() +  geom_tile(data = domain, aes(x = x, y = y, fill = null)) +
  coord_quickmap() +  labs(x="",y="",fill="Richness")+
  scale_fill_gradientn(colours = viridis::viridis(99),na.value = "lightgray",limits=c(0,ma))+the+
  annotate("text",x=-72,y=25,size=6,label=paste0("r = ",round(cor(na.omit(domain[,c("CurrentRichness","null")])),2)[1,2]))+
  annotate("text",x=-72,y=30,size=4,label="ER Projected \nnull model",color="red")
  
p3=ggplot() +  geom_tile(data = domain, aes(x = x, y = y, fill = clim)) +
  coord_quickmap() +  labs(x="",y="",fill="Richness")+
  scale_fill_gradientn(colours = viridis::viridis(99),na.value = "lightgray",limits=c(0,ma))+the+
  annotate("text",x=-72,y=25,size=6,label=paste0("r = ",round(cor(na.omit(domain[,c("CurrentRichness","clim")])),2)[1,2]))+
  annotate("text",x=-72,y=30,size=4,label="ER Projected\n(Climate \ndistance)",color="red")
  
p4=ggplot() +  geom_tile(data = domain, aes(x = x, y = y, fill = soil)) +
  coord_quickmap() +  labs(x="",y="",fill="Richness")+
  scale_fill_gradientn(colours = viridis::viridis(99),na.value = "lightgray",limits=c(0,ma))+the+
  annotate("text",x=-72,y=25,size=6,label=paste0("r = ",round(cor(na.omit(domain[,c("CurrentRichness","soil")])),2)[1,2]))+
  annotate("text",x=-72,y=30,size=4,label="ER Projected\n(soil \ndistance)",color="red")

p5=ggplot() +  geom_tile(data = domain, aes(x = x, y = y, fill = clim.soil)) +
  coord_quickmap() +  labs(x="",y="",fill="Richness")+
  scale_fill_gradientn(colours = viridis::viridis(99),na.value = "lightgray",limits=c(0,ma))+the+
  annotate("text",x=-72,y=25,size=6,label=paste0("r = ",round(cor(na.omit(domain[,c("CurrentRichness","clim.soil")])),2)[1,2]))+
  annotate("text",x=-72,y=30,size=4,label="ER Projected\n(Climate and soil \ndistance)",color="red")
  
p6=ggplot() +  geom_tile(data = domain, aes(x = x, y = y, fill = human)) +
  coord_quickmap() +  labs(x="",y="",fill="Richness")+
  scale_fill_gradientn(colours = viridis::viridis(99),na.value = "lightgray",limits=c(0,ma))+the+
  annotate("text",x=-72,y=25,size=6,label=paste0("r = ",round(cor(na.omit(domain[,c("CurrentRichness","human")])),2)[1,2]))+
  annotate("text",x=-72,y=30,size=4,label="ER Projected\n(human disturbance)",color="red")
  
p7=ggplot() +  geom_tile(data = domain, aes(x = x, y = y, fill = clim.soil.pop)) +
  coord_quickmap() +  labs(x="",y="",fill="Richness")+
  scale_fill_gradientn(colours = viridis::viridis(99),na.value = "lightgray",limits=c(0,ma))+the+
  annotate("text",x=-72,y=25,size=6,label=paste0("r = ",round(cor(na.omit(domain[,c("CurrentRichness","clim.soil.pop")])),2)[1,2]))+
  annotate("text",x=-72,y=30,size=4,label="ER Projected\n(climate, soil and\nhuman disturbance)",color="red")

p8=ggplot() +  geom_tile(data = domain, aes(x = x, y = y, fill = sp)) +
  coord_quickmap() +  labs(x="",y="",fill="Richness")+
  scale_fill_gradientn(colours = viridis::viridis(99),na.value = "lightgray",limits=c(0,ma))+the+
  annotate("text",x=-72,y=25,size=6,label=paste0("r = ",round(cor(na.omit(domain[,c("CurrentRichness","sp")])),2)[1,2]))+
  annotate("text",x=-72,y=30,size=4,label="ER Projected\n(species assemblages)",color="red")
  
p9=ggplot() +  geom_tile(data = domain, aes(x = x, y = y, fill = phylo2)) +
  coord_quickmap() +  labs(x="",y="",fill="Richness")+
  scale_fill_gradientn(colours = viridis::viridis(99),na.value = "lightgray",limits=c(0,ma))+the+
  annotate("text",x=-72,y=25,size=6,label=paste0("r = ",round(cor(na.omit(domain[,c("CurrentRichness","phylo2")])),2)[1,2]))+
  annotate("text",x=-72,y=30,size=4,label="ER Projected\n(phylogenetic assemblages)",color="red")
  
p10=ggplot() +  geom_tile(data = domain, aes(x = x, y = y, fill = trait)) +
  coord_quickmap() +  labs(x="",y="",fill="Richness")+
  scale_fill_gradientn(colours = viridis::viridis(99),na.value = "lightgray",limits=c(0,ma))+the+
  annotate("text",x=-72,y=25,size=6,label=paste0("r = ",round(cor(na.omit(domain[,c("CurrentRichness","trait")])),2)[1,2]))+
  annotate("text",x=-72,y=30,size=4,label="ER Projected\n(trait \ndistance)",color="red")

p11=ggplot() +  geom_tile(data = domain, aes(x = x, y = y, fill = trait.turnover)) +
  coord_quickmap() +  labs(x="",y="",fill="Richness")+
  scale_fill_gradientn(colours = viridis::viridis(99),na.value = "lightgray",limits=c(0,ma))+the+
  annotate("text",x=-72,y=25,size=6,label=paste0("r = ",round(cor(na.omit(domain[,c("CurrentRichness","trait.turnover")])),2)[1,2]))+
  annotate("text",x=-72,y=30,size=4,label="ER Projected\n(trait turnover)",color="red")
    
ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,nrow=2,ncol=5,labels="auto",font.label = list(size = 20),common.legend=TRUE,legend="right")	

#####################################################################
### calculate expected alien range size for different ER thresholds
#####################################################################
library(dplyr);
load("data analysis/connectedGraph.ena.rda")
#inputER<-c(0.1,0.2,0.3,seq(0.4,1,0.01))
inputER<-c(0.1,0.3,0.4,0.5,seq(0.55,0.98,0.01))
#inputER<-seq(0.55,0.95,0.03)
#RangeFragments<-get(load("data analysis/AlienRangeFragments.plants.ENA.RData"))
load("data analysis/AlienRangeFragments.alltaxa.ENA.RData")

get.simrange=function(x,RangeFragments,simmat,g1cc,mat.name){
	require(dplyr)
	fragmentTable1=data.frame(fragBinomial=unique(names(RangeFragments)),fragSpecies=1:length(unique(names(RangeFragments))))
	fragmentTable2=data.frame(fragNum=1:length(RangeFragments),fragSize=unlist(lapply(RangeFragments,length)),fragBinomial=names(RangeFragments))
	FragmentTable=fragmentTable2%>%left_join(fragmentTable1,by="fragBinomial")
	#for(x in 1:inputER){  
	source("data analysis/Functions.R")
	output<-runConnectedComponentsModel(threshold=x,nRandomStarts=100,fragmentTable=FragmentTable,RangeFragments=RangeFragments,CommSimilarity=simmat,g1cc=g1cc)
	#save the simulated ranges
	fragmentTable<-output[[1]]
	simRanges<-output[[2]]
	save(simRanges,fragmentTable,file=paste0("data analysis/simulation.ena.predict/PredictedRangeSize_",x,".",mat.name,".rda")) 
	rm(fragmentTable1,fragmentTable2,output,RangeFragments,fragmentTable,simRanges);gc()
	}
	
library(parallel)
no_cores <- detectCores() - 1	
mat.list=c("sp","clim","soil","phylo2","human","trait","clim.soil","clim.soil.pop","trait.turnover")	
	
for (mat in mat.list){
	if (mat%in%names(envir.d)){
		simmat=envir.d[[mat]]
	}else{
		simmat=get(load(paste0("data analysis/CommSimilarity.ENA.",mat,".rda")))
	}
	if (mat%in%"sp") simmat=1-simmat
	if (max(simmat,na.rm=T)>1) simmat=scales::rescale(simmat)
	print(mat)
	mycl <- makePSOCKcluster(no_cores);
	parSapply(cl=mycl,X=inputER,get.simrange,RangeFragments.ena.inv,simmat,g1cc,mat)
	stopCluster(mycl)
	gc()	
}

# fit quantile regressions between observed and predicted range sizes for different ER thresholds and, for quantile egressions, tau values 
splist=read.csv("data analysis/US-RIIS.csv")%>%select(kingdom,phylum,class,scientificName)%>%filter(scientificName%in%names(RangeFragments.ena.inv))
taxa=c("Chordata","Arthropoda","Other animal","Plantae")

get.ERQ=function(mat,taxa,inputER){
	require(Metrics);require(quantreg);
	# different ER thresholds and tau values to explore
	input<-expand.grid(inputER=inputER,
						#tau=c(0.5,0.9,0.95,0.99),
						   tau=0.99,
						   #rmse=NA,,r2=NA,
						   rmseQ=NA,
						   rmseQspecies=NA)
	input.taxa[[taxa.i]]=input
	# # store predicted range sizes and regression stats that we will need for plotting
	# rsMeanMat<-matrix(ncol=nrow(input),nrow=frags)
	for(x in 1:nrow(input)){
	  load(paste0("data analysis/simulation.ena.predict/PredictedRangeSize_",input$inputER[x],".",mat,".rda"))	  
	  fragmentTable$fragSizeLOG<-log10(fragmentTable$fragSize)
	  fragmentTable$predfragSizeLOG<-log10(fragmentTable$PredictedFragSize)
	  #rsMeanMat[,x]<-fragmentTable$PredictedFragSize
	  predfragSizeLOG<-fragmentTable$predfragSizeLOG  
	  for (taxa.i in taxa){
	  fragmentTable.t=fragmentTable%>%filter(fragBinomial%in%)
	  
	  
	  # quantile regression predicting maximum range size: 在传统回归中，我们构建回归模型由自变量求出因变量的条件期望；而在分位数回归中，我们构建回归模型由自变量求出因变量的条件分位数
	  rqfit <- rq(fragSizeLOG ~ predfragSizeLOG, data = fragmentTable,tau=input$tau[x],method="fn")
	  predicted<-predict(rqfit,newdata =as.data.frame(predfragSizeLOG))
	  input$rmseQ[x]<-rmse(predfragSizeLOG, predicted)
	  
	  # quantile regression predicting maximum range size - at species level
	  predfragSizeLOG<-log10(sapply(split(fragmentTable$PredictedFragSize,fragmentTable$fragSpecies),sum))
	  fragSizeLOG<-log10(sapply(split(fragmentTable$fragSize,fragmentTable$fragSpecies),sum))
	  rqfit <- rq(fragSizeLOG ~ predfragSizeLOG, tau=input$tau[x],method="fn")
	  predicted<-predict(rqfit,newdata =as.data.frame(predfragSizeLOG))
	  input$rmseQspecies[x]<-rmse(predfragSizeLOG, predicted)
	  }
	  # # linear model predicting mean range size
	  # mod<-lm(fragSizeLOG ~ predfragSizeLOG, data = fragmentTable)
	  # input$r2[x]<-summary(mod)$r.squared
	  # predicted<-predict(mod,newdata =as.data.frame(predfragSizeLOG))
	  # input$rmse[x]<-rmse(fragmentTable$fragSizeLOG, fragmentTable$predfragSizeLOG)
	  
	    
	}
	save(input,file=paste0("data analysis/simulation.ena.predict/",taxa.i,"/BestERvalues_",mat,".rda"))
}

mycl <- makePSOCKcluster(no_cores);
for (taxa.i in taxa) parSapply(cl=mycl,mat.list,get.ERQ,taxa.i,inputER)
stopCluster(mycl)

## plot root mean square error of different ERs under different climate scenerio
library(ggpubr);library(dplyr)
mat.list=c("clim","soil","clim.soil","human","clim.soil.pop","sp","phylo2","trait","trait.turnover")	
taxa=c("Chordata","Arthropoda","Other animal","Plantae")

rmse=c();
for (taxa.i in taxa){
for (mat in mat.list){
	load(paste0("data analysis/simulation.ena.predict/",taxa.i,"/BestERvalues_",mat,".rda"))
	rmse=rbind(rmse,data.frame(taxa=taxa.i,model=mat,input[,c("inputER","rmseQ","rmseQspecies")]))
}
}
stat.ERQ=rmse%>%group_by(taxa,model)%>%summarize(bestERQ=inputER[which(rmseQ==min(rmseQ))],bestERQspecies=inputER[which(rmseQspecies==min(rmseQspecies))])#%>%
	#group_by(taxa)%>%summarize(bestERQ=round(mean(bestERQ),2),bestERQspecies=round(mean(bestERQspecies),2))
# taxa         model bestERQ bestERQspecies
  # <chr>        <chr>   <dbl>          <dbl>
# 1 Arthropoda   sp        0.3           0.81
# 2 Chordata     sp        0.3           0.81
# 3 Other animal sp        0.5           0.3 
# 4 Plantae      sp        0.1           0.82

theme=theme_bw()+theme(axis.text = element_text(size=15,color='black'),
			axis.title = element_text(size=18,color='black'),
			legend.background=element_rect(fill='transparent'),
			legend.text=element_text(face="bold",size=15),
			legend.title=element_text(face="bold",size=18))
rmse$model=factor(rmse$model,levels=mat.list)

p=ggplot(rmse,aes(x = inputER,y=rmseQspecies,colour=model))+
	geom_line(size=0.8,linetype=1,alpha=0.8)+
	#scale_color_discrete(labels = c("Null model", "Climatic distance","Species composition","Phylogenetic composition","Trait distance","Trait turnover"))+
	labs(x="Environmental resistance values",y="Root mean square error",colour="Model") +theme+theme(legend.position =c(0.1,0.85))

plt=list()
for (taxa.i in taxa){
p=ggplot(subset(rmse,rmse$taxa==taxa.i),aes(x = inputER,y=rmseQspecies,colour=model))+
	geom_line(size=0.8,linetype=1,alpha=0.8)+
	#scale_color_discrete(labels = c("Null model", "Climatic distance","Species composition","Phylogenetic composition","Trait distance","Trait turnover"))+
	labs(x="Environmental resistance values",y="Root mean square error",colour="Model") +theme+theme(legend.position =c(0.3,0.75))
for (mat in 1:length(mat.list)) {
	bestERQ<-stat.ERQ$bestERQspecies[which(stat.ERQ$model==mat.list[mat]&stat.ERQ$taxa==taxa.i)]
	p=p+geom_vline(xintercept=bestERQ,col=scales::hue_pal()(length(mat.list))[mat],size=1,linetype="longdash")+
		annotate("text", x=bestERQ-0.03 , y=2.5+mat*0.1 ,size=5,label=bestERQ,color=scales::hue_pal()(length(mat.list))[mat])
}
plt[[taxa.i]]=p
}
ggarrange(plt[[1]],plt[[2]],plt[[3]],plt[[4]],labels="auto",common.legend=TRUE)

########################################################################################  #####
# generate matrix of invasion probabilities for each species in each cell using best ER model of each species
########################################################################################  
library(dplyr)
load("data analysis/Domain.ena.rda")
load("data analysis/AlienRangeFragments.alltaxa.ENA.RData")
splist=read.csv("data analysis/US-RIIS.csv")%>%select(kingdom,phylum,class,scientificName)%>%filter(scientificName%in%names(RangeFragments.ena.inv))
taxa=c("Chordata","Arthropoda","Other animal","Plantae")
nRandomStarts<-100
AlienfragmentTable=data.frame(fragNum=1:length(RangeFragments.ena.inv),fragSize=unlist(lapply(RangeFragments.ena.inv,length)),fragBinomial=names(RangeFragments.ena.inv))

if(taxa.i=="Plantae"){
	sp=splist[splist$kingdom=="Plantae","scientificName"]
}else{
	sp=splist[splist$phylum==taxa.i,"scientificName"]
	if (taxa.i%in%"Other animal") sp=splist[splist$kingdom=="Animalia"&(!splist$phylum%in%c("Chordata","Arthropoda")),"scientificName"]
}
nSpecies<-length(sp)

load("data analysis/speciesANDbestmodel.rda")
score.count=score3%>%group_by(Species)%>%tally()%>%filter(n>1)
score3[score3$Species%in%score.count$Species,"bestmodel"]="sp"
score=unique(score3)

# bestERQ=0.67;mat="sp"
# bestERQ=0.88;mat="clim.pca"
# load(paste0("data analysis/simulation.ena.predict/PredictedRangeSizeConnectedComponentsERThreshold_",bestERQ,".",mat,".rda"))

InvasionProbabilityMatrix<-matrix(0,ncol=nrow(domain.ena),nrow=nSpecies)
for(i in 1:nSpecies){ # loop over each species  
  focRowAlien<-which(AlienfragmentTable$fragBinomial==sp[i])
  nFrags<-length(focRowAlien)
  
  # bestmodel=score%>%filter(Species==unique(AlienfragmentTable[focRowAlien,"fragBinomial"]))
  # mat=bestmodel$bestmodel
  # if (mat=="clim.pca") bestERQ=0.88
  # if (mat=="phylo2") bestERQ=0.84
  # if (mat=="sp") bestERQ=0.67
  # if (mat%in%c("trait","trait.turnover")) bestERQ=0.99  
  # load(paste0("data analysis/simulation.ena.predict/PredictedRangeSizeConnectedComponentsERThreshold_",bestERQ,".",mat,".rda"))
  
  # for each alien range fragment calculate the probability of invading each cell
  InvasionProbs<-matrix(0,ncol=nrow(domain.ena),nrow=nFrags)
  for(x in 1:nFrags){
    if(fragmentTable$RandomStart[focRowAlien[x]]==1){
      fragSize<-AlienfragmentTable$fragSize[focRowAlien[x]]
      InvasionProbs[x,as.numeric(names(simRanges[[focRowAlien[x]]]))]<-as.numeric(simRanges[[focRowAlien[x]]])/min(nRandomStarts,fragSize)
    }else{
      InvasionProbs[x,simRanges[[focRowAlien[x]]]]<-1
    }
  }
  
  # the overall probability of an alien species invading a cell (a cell can be invaded from multiple fragments)
  InvasionProbVec<-1-matrixStats::colProds(1-InvasionProbs)
  
  # # set invasion probabilities of cells that are in the species' native range to zero
  # focRowNative<-which(NativefragmentTable$fragSpecies==i)
  # focNativeRange<-unique(unlist(NativeRangeFragments[focRowNative]))
  # InvasionProbVec[focNativeRange]<-0
  InvasionProbabilityMatrix[i,]<-InvasionProbVec
  
  print(i)
}

#devided by range size
rangesize=AlienfragmentTable%>%group_by(fragBinomial)%>%summarize(ranges=sum(fragSize))%>%as.data.frame()
cat1=subset(rangesize,rangesize[,"ranges"]>=quantile(rangesize[,"ranges"],0.75))[,"fragBinomial"]
cat2=subset(rangesize,rangesize[,"ranges"]<=quantile(rangesize[,"ranges"],0.25))[,"fragBinomial"]

domain=domain.ena
# Projected alien richness
domain$InvasionProb<-colSums(InvasionProbabilityMatrix)
# Current alien richness
RangeFragments=subset(RangeFragments.ena.inv,names(RangeFragments.ena.inv)%in%sp)
AlienRichness<-table(unlist(RangeFragments))
domain$CurrentRichness<-NA
domain$CurrentRichness[match(as.numeric(names(AlienRichness)),domain$domain.ID)]<-as.numeric(AlienRichness)


#plot richness
library(ggpubr)
the=theme_bw()+theme(legend.background=element_rect(fill='transparent'),
		legend.text=element_text(size=10),
		legend.title=element_text(face="bold",size=10))
		
p1=ggplot() +  geom_tile(data = domain, aes(x = x, y = y, fill = InvasionProb)) +
  coord_quickmap() +  labs(x="",y="",fill="InvasionProb")+
  scale_fill_gradientn(colours = viridis::viridis(99),na.value = "lightgray")+the+theme(legend.position=c(0.8,0.25))	
p2=ggplot() +  geom_tile(data = domain, aes(x = x, y = y, fill = CurrentRichness)) +
  coord_quickmap() +  labs(x="",y="",fill="Richness")+
  scale_fill_gradientn(colours = viridis::viridis(99),na.value = "lightgray")+the+theme(legend.position=c(0.8,0.25))

ggarrange(p1,p2,labels="auto",font.label = list(size = 20))	

############################################################################################################
####get future climate data from worldclim and project future distributions and conduct SDM ################################
############################################################################################################
#ref for select the best CMIP6 models: Projected Changes in Temperature and Precipitation Over the United States, Central America, and the Caribbean in CMIP6 GCMs (Table 1)
# ACCESS-CM2;EC-Earth3-Veg;MPI-ESM1-2-HR;UKESM1-0-LL are models with both low terrestrial precipitation and temperature bias
library(sf);library(dplyr);library(raster)
clim.list=list.files("data analysis/future climate/")

load("data analysis/Domain.ena.rda")
ENA.name=c("ME","NH","VT","NY","MA","RI","CT","NJ","PA","DE","MD","MI","OH","IN","IL","WI","WV","VA","NC","TN","KY","SC","GA","AL","MS","FL")
ENA <- read_sf('data analysis/cb_2018_us_state_20m/cb_2018_us_state_20m.shp')%>%filter(STUSPS%in%ENA.name)%>%as(., "Spatial")
# Use Variable inflation factor (VIF) to select the best variables for your model.
library(usdm)	
for (i in clim.list){	
	rst.10m <- brick(paste0("data analysis/future climate/",i))
	ENA.rst=ENA%>%raster::mask(rst.10m,.)%>%raster::crop(extent(ENA))#%>%aggregate(fact=1.5)
	Domain.ena=raster::as.data.frame(ENA.rst,xy=TRUE)#
	#remove ocean cells
	Domain.ena=na.omit(Domain.ena)	
	#select VIF<10 envars
	v=vifstep(Domain.ena[,-c(1:2)],th=10)
	var.list=attributes(v)$results[,1]
	ENA.rst2=ENA.rst[[var.list]]
	Domain.ena=Domain.ena%>%left_join(domain.ena[,c("x","y","domain.ID")],by=c("x","y"))
	print (c(i,dim(Domain.ena)))
		
	outfile=gsub("wc2.1_10m_bioc_","",i);outfile=gsub(".tif","",outfile)
	re=list(ENA.rst2,Domain.ena)
	save(re,file=paste0("data analysis/future climate/",outfile,".rda"))
	save()
}

# load sp dis
library(dplyr)
toload=list.files("data analysis/future climate/")%>%.[grep(".rda",.)]
#rangesize=rangeData%>%group_by(Tip_Label)%>%summarize(rangesize=length(domain.ID))
#splist=rangesize%>%filter(rangesize==2)%>%as.data.frame()%>%.[,1]
#for (i in 1:3) 	p_df=get.maxent(splist[i],re[[1]],rangeData)	
# conduct Maxent	
get.maxent=function(p,ENA.rst2){
	# ENA.rst2=re[[1]]
	# j=splist[1]
	# Set up java memory 
	# p=rangeData[[1]]
	#options(java.parameters = "- Xmx16g") # increase memory that can be used
	require(dplyr);require(ENMeval)
	#p=rangeData%>%filter(Tip_Label%in%j)%>%dplyr::select(x,y)
	if (nrow(p)<3) p_df=data.frame(p[,2:3],predic=1,Species=p[,1],type="real") else {
	eval1 <- try({ENMeval::ENMevaluate(occ = p[,2:3], 
								  env = ENA.rst2,
								  tune.args = list(fc = c("L","Q"), rm = 1:raster::nlayers(ENA.rst2)), # test the feature classes L = linear and Q = quadratic
								  partitions = "block",#
								  parallel = FALSE,
								  algorithm = 'maxent.jar')}, silent = TRUE)
	if (class(eval1)== "try-error"){
		p_df=data.frame(p[,2:3],predic=1,Species=p[,1],type="real")
	}else{		
	# evaldis <- dismo::maxent(x = ENA.rst2, p = p,args = c("-X", 20,"jackknife"))#20% reserved for testing
                                  # "randomseed", "randomtestpoints=25"), # randomly set aside 25% of the sample records for testing.
                         # removeDuplicates = TRUE)
						 
	### Inspect the results
	### Identify the best model
	#### selecting models with the lowest average test omission rate and 
	#### the highest average validation AUC and lowest AICc
	results <- eval.results(eval1)
	opt.seq <- results %>% 
            dplyr::filter(or.10p.avg == min(or.10p.avg,na.rm=T)) %>% 
            dplyr::filter(auc.val.avg == max(auc.val.avg,na.rm=T))
	if (nrow(opt.seq)>1) opt.seq = opt.seq	%>% dplyr::filter(delta.AICc == min(delta.AICc,na.rm=T))
	#plot(eval1@predictions[[which (eval1@results$delta.AICc == 0) ]])
	### Subset model	
	mod.seq <- eval.models(eval1)[[opt.seq$tune.args[1]]]
	#predic <- maxent::predict(mod.seq, ENA.rst2) 
	#p_df <-  as.data.frame(predic, xy = TRUE)%>%na.omit()
	env=raster::as.data.frame(ENA.rst2,xy=TRUE)%>%na.omit()
	predic <- try(enmSdmX::predictMaxEnt(mod.seq, env[,-c(1:2)]))
	if (class(predic)== "try-error"){
		p_df=data.frame(p[,2:3],predic=1,Species=p[,1],type="real")
	}else{	
		p_df <- cbind(env,predic)%>%dplyr::select(x,y,predic)%>%filter(predic>0.1)
		if (nrow(p_df)==0) p_df=data.frame(p[,2:3],predic=1,Species=p[,1],type="real") else p_df$Species=unique(p[,1]);p_df$type="sdm"
	}
	}
	}
	return(p_df)#write.csv(p_df,paste("sdm.map/p_df",j,reg,"csv",sep="."))	
}		

# for (j in splist){
# p_df=get.maxent(j,ENA.rst2,rangeData)
# }
library(parallel)
no_cores <-40#detectCores() - 1
for (i in toload){
	load(paste0("data analysis/future climate/",i))
	rangeData=get(load("data analysis/NativeRangeFragments.plants.ENA.clean.RData"))%>%left_join(re[[2]],by="domain.ID")%>%na.omit()%>%
		dplyr::select(Tip_Label,x,y)%>%split(.,.[,1])
	#rangeData=	rangeData[8041:10616]
	print(i)
	nstart=1	
	#if (i==toload[5]) nstart=8041 else nstart=1	
	seql=c(seq(nstart,length(rangeData),no_cores*67),length(rangeData))#30
	seql2=lapply(2:length(seql),function(i,seql){if (i<length(seql)){c(seql[i-1],(seql[i]-1))
		} else {
		c(seql[i-1],seql[i])
		}},seql)
	
	for (j in 1:length(seql2)){
		mycl <- makePSOCKcluster(no_cores);		
		p_df=parLapply(cl=mycl,rangeData[seql2[[j]][1]:seql2[[j]][2]],get.maxent,re[[1]])
		#p_df=parLapply(cl=mycl,rangeData,get.maxent,re[[1]])
		stopCluster(mycl)
		save(p_df,file=paste0("data analysis/SDM/sdm_",seql2[[j]][2],"_",i));
		#save(p_df,file=paste0("data analysis/SDM/sdm_10616_",i));
		rm(p_df);gc()
		print(j);				
	}	
	
}


### Visualize
toload.sdm=list.files("data analysis/SDM/")
load(paste0("data analysis/SDM/",toload.sdm[1]))
library(ggplot2);library(dplyr)
p=rangeData[[j]]
prd=re[[2]]%>%left_join(p_df,by=c("x","y"))#%>%filter(predic>=0.9)
### Plot
ggplot() + geom_tile(data = prd, aes(x = x, y = y, fill = predic))+
	geom_point(data= p,mapping = aes(x = x, y = y),col='red', cex=1,alpha=0.8) +
	coord_quickmap() +  theme_bw() + 
	scale_fill_gradientn(colours = viridis::viridis(99),na.value = "black")	
	
##combine all the predicted distribution within each clim scenerio
load("data analysis/Domain.ena.rda")
require(dplyr)
domain.ena=domain.ena%>%dplyr::select(domain.ID,x,y)
# clim.model=c("ACCESS-CM2","EC-Earth3-Veg","MPI-ESM1-2-HR","UKESM1-0-LL")
# period=c("2021-2040","2041-2060","2061-2080")
# ssp=c("ssp126","ssp245","ssp370","ssp585")
file.list=list.files("data analysis/future climate/")%>%.[grep(".rda",.)]%>%gsub(".rda","",.)
#cal.comm=function(name,domain.ena){
	for (name in file.list){
	#name=file.list[1]
	require(dplyr)
	get.bidis=function(p_df,thes,domain.ena){
		p_df2=p_df[p_df$predic>=thes,]
		p_df2=p_df2%>%left_join(domain.ena,by=c("x","y"))%>%na.omit()%>%select(Species,domain.ID,predic)%>%distinct()
		colnames(p_df2)[1]="Tip_Label"
		return(p_df2)
	}
	spdis.list=list.files("data analysis/SDM/")%>%.[grep(name,.)]		
	spdis=c()
	for (i in spdis.list){
		spdis0=get(load(paste0("data analysis/SDM/",i)))
		spdis.t=do.call(rbind,lapply(spdis0,get.bidis,thes=0.1,domain.ena))
		spdis=rbind(spdis,spdis.t)
	}
	save(spdis,file=paste0("data analysis/SDM.spdis/sdm.spdis_",name,".rda"))
	rm(spdis,spdis0,spdis.t);gc()
	print(name)
	}
#}
# library(parallel)
# no_cores <-detectCores() - 1
# mycl <- makePSOCKcluster(no_cores);	
# parSapply(cl=mycl,file.list,cal.comm,domain.ena)
# stopCluster(mycl)

## caculate future CommSimilarity
#species similarity
library(dplyr)
file.list=list.files("data analysis/SDM.spdis/")%>%gsub(".rda","",.)%>%gsub("sdm.spdis_","",.)
get.sim.sp=function(name){
	for(name in file.list){
	require(dplyr)
	nCells=8490	
	rangeData=data.frame(domain.ID=1:nCells)%>%left_join(get(load(paste0("data analysis/SDM.spdis/sdm.spdis_",name,".rda"))),by="domain.ID")
	speciesByCell<-split(rangeData$Tip_Label,rangeData$domain.ID)
	predictByCell<-split(rangeData$predic,rangeData$domain.ID)
	rm(rangeData);gc()
	CommSimilarity<-matrix(ncol=nCells,nrow=nCells)
	for(i in 1:nCells){ # loop over each cell
	 cellIspecies<-na.omit(speciesByCell[[i]]) # species in cell i 			  
	  if(length(cellIspecies)>0){
	    focalN<-sum(predictByCell[[i]])
		for(j in 1:nCells){
		  if (i==j) CommSimilarity[i,j]<-NA else{		
		  cellJspecies<-na.omit(speciesByCell[[j]])
		  if(length(cellJspecies)>0){		  
			share=which(cellIspecies%in%intersect(cellIspecies,cellJspecies))
			CommSimilarity[i,j]<-sum(predictByCell[[i]][share])/focalN			
		  }else{
			CommSimilarity[i,j]<-0 # if the richness of cell j is zero, set simmilarity to zero
		  }
		 }
		}
	  }else{
		CommSimilarity[i,]<-0 # if the richness of cell i is zero, set simmilarity to zero
	  }
	  print(i)
	}
	rownames(CommSimilarity)=colnames(CommSimilarity)=1:8490
	save(CommSimilarity,file=paste0("data analysis/SDM.sp.simm/CommSimilarity_",name,".rda"))
	rm(speciesByCell,predictByCell);gc()
	}	
}
no_cores <-24#detectCores() - 1
mycl <- makePSOCKcluster(no_cores);	
parSapply(cl=mycl,file.list,get.sim.sp)
stopCluster(mycl)

#climate similarity
load("data analysis/Domain.ena.rda")
cal.climdis=function(name,domain.ena){
	require(dplyr)
	load(paste0("data analysis/future climate/",name,".rda"))
	bioclim=domain.ena%>%dplyr::select(-x,-y,-starts_with("bio"))%>%left_join(re[[2]],by="domain.ID")#%>%dplyr::select(domain.ID,ele,starts_with("bio"))
	miss=bioclim%>%filter(is.na(bio01))
	imp=function(i,bioclim){
		cells.row=(miss[i,]$row-5):(miss[i,]$row+5)
		cells.col=(miss[i,]$col-5):(miss[i,]$col+5)
		impute=bioclim%>%filter(row%in%cells.row&col%in%cells.col)%>%dplyr::select(starts_with("bio"))%>%colMeans(na.rm =TRUE)
		return(impute)
	}
	miss.imp=do.call(rbind,lapply(1:nrow(miss),imp,bioclim))
	miss.imp=cbind(miss%>%dplyr::select(-starts_with("bio")),miss.imp)
	bioclim2=rbind(bioclim%>%filter(!is.na(bio01)),miss.imp)
	bioclim2=bioclim2[match(bioclim$domain.ID,bioclim2$domain.ID),]%>%dplyr::select(domain.ID,ele,starts_with("bio"))
	pca=ade4::dudi.pca(bioclim2[,-c(1:2)],center = TRUE,scale = T,scannf = FALSE,nf = 5)
	EnvironmentData=cbind(scale(bioclim2[,-1]),pca$li)  
	rownames(EnvironmentData)=bioclim2$domain.ID
	clim=dist(EnvironmentData[,c("Axis1","Axis2","Axis3","Axis4","Axis5","ele")])%>%as.matrix()%>%scales::rescale();diag(clim)=NA
	save(clim, file=paste0("data analysis/SDM.clim.simm/CommSimilarity.clim.",name,".rda"))
	}
no_cores <-detectCores() - 1
mycl <- makePSOCKcluster(no_cores);	
parSapply(cl=mycl,file.list,cal.climdis,domain.ena)
stopCluster(mycl)	

########################################################################################################
### calculate expected alien range size for different ER thresholds under future climate and obtain the best ERQ
########################################################################################################
library(dplyr);
load("data analysis/connectedGraph.ena.rda")
#inputER<-c(0.1,0.2,0.3,seq(0.4,1,0.01))
inputER<-c(0.1,0.3,0.4,0.5,seq(0.55,0.98,0.01))
#RangeFragments<-get(load("data analysis/AlienRangeFragments.plants.ENA.RData"))
load("data analysis/AlienRangeFragments.alltaxa.ENA.RData")
splist=read.csv("data analysis/US-RIIS.csv")%>%select(kingdom,phylum,class,scientificName)%>%filter(scientificName%in%names(RangeFragments.ena.inv))
taxa=c("Plantae","Chordata","Arthropoda","Other animal")

toload<-list.files("data analysis/SDM.sp.simm/")	
library(parallel)
no_cores <- detectCores() - 1
for (mat in 1:9){#length(toload)){
	for (taxa.i in taxa){
	print(taxa.i)
	mycl <- makePSOCKcluster(no_cores);
	mat.name=gsub('CommSimilarity_', '', toload[mat])%>%gsub('.rda', '', .)
	simmat=get(load(paste0("data analysis/SDM.sp.simm/",toload[mat])))
	print(mat.name)
	parSapply(cl=mycl,X=inputER,get.simrange,taxa.i,splist,RangeFragments.ena.inv,simmat,g1cc,mat.name)
	stopCluster(mycl)
	rm(simmat);gc()
	}
}

mycl <- makePSOCKcluster(no_cores);
mat.list<-list.files("C:/Users/yunpeng.liu/OneDrive - University of Florida/UF postdoc/disease macroecology/Lovell2021/data analysis/SDM.sp.simm/")%>%gsub('.rda', '', .)%>%gsub('CommSimilarity_', '', .)
for (taxa.i in taxa) parSapply(cl=mycl,mat.list[1:9],get.ERQ,taxa.i,inputER)
stopCluster(mycl)

library(ggpubr);library(dplyr)
rmse=c();
for (taxa.i in taxa){
for (mat in mat.list[1:9]){
	load(paste0("data analysis/simulation.future.predict/",taxa.i,"/BestERvalues_",mat,".rda"))
	rmse=rbind(rmse,data.frame(taxa=taxa.i,model=mat,input[,c("inputER","rmseQ","rmseQspecies")]))
}
}
stat.rmse=rmse%>%group_by(taxa,model)%>%summarize(bestERQ=inputER[which(rmseQ==min(rmseQ))],bestERQspecies=inputER[which(rmseQspecies==min(rmseQspecies))])%>%
	group_by(taxa)%>%summarize(bestERQ=round(mean(bestERQ),2),bestERQspecies=round(mean(bestERQspecies),2))

 # taxa         bestERQ bestERQspecies
# 1 Arthropoda      0.87           0.84
# 2 Chordata        0.68           0.78
# 3 Other animal    0.77           0.91
# 4 Plantae         0.75           0.87

###################################################################################
### calculate expected alien range size under future climate changes using the best ERQ
###################################################################################
library(dplyr);
load("data analysis/connectedGraph.ena.rda")
#toload=c(list.files("data analysis/SDM.sp.simm/"),list.files("data analysis/SDM.clim.simm/"))%>%.[grep("CommSimilarity.",.)]%>%gsub(".rda","",.)
toload=list.files("C:/Users/yunpeng.liu/OneDrive - University of Florida/UF postdoc/disease macroecology/Lovell2021/data analysis/SDM.sp.simm/")%>%gsub(".rda","",.)%>%.[grep("ACCESS-CM2",.)]
source("data analysis/Functions.R")
load("data analysis/AlienRangeFragments.alltaxa.ENA.RData")
splist=read.csv("data analysis/US-RIIS.csv")%>%select(kingdom,phylum,class,scientificName)%>%filter(scientificName%in%names(RangeFragments.ena.inv))
taxa=c("Chordata","Arthropoda","Other animal","Plantae")
nRandomStarts<-100

get.sdm.predict=function(mat,toload,AlienfragmentTable,RangeFragments,g1cc,taxa.i){
#for (mat in 1:length(toload)){
	#source("Functions.R")
	require(dplyr)
	source("data analysis/Functions.R")
	mat.name=gsub('CommSimilarity.', '', toload[mat])%>%gsub('.rda', '', .)%>%gsub('CommSimilarity.', '', .)
	load(paste0("C:/Users/yunpeng.liu/OneDrive - University of Florida/UF postdoc/disease macroecology/Lovell2021/data analysis/SDM.sp.simm/CommSimilarity_",mat.name,".rda"))
	inputER=ifelse(taxa.i=="Arthropoda",0.84,ifelse(taxa.i=="Chordata",0.78,ifelse(taxa.i=="Other animal",0.91,0.87)))
	print(mat.name)
	output<-runConnectedComponentsModel(threshold=inputER,nRandomStarts=100,
										  fragmentTable=AlienfragmentTable,
										  RangeFragments=RangeFragments,
										  CommSimilarity=CommSimilarity,
										  g1cc=g1cc)
	#save the simulated ranges
	fragmentTable<-output[[1]]
	simRanges<-output[[2]]	
	save(simRanges,fragmentTable,file=paste("data analysis/SDM.bestER/PredictedER",taxa.i,mat.name,"rda",sep=".")) 
	rm(output,fragmentTable,simRanges)
	gc()
#}
}

library(parallel)
no_cores <- detectCores() - 1
for (taxa.i in taxa){
	if(taxa.i=="Plantae"){
		RangeFragments=subset(RangeFragments.ena.inv,names(RangeFragments.ena.inv)%in%splist[splist$kingdom=="Plantae","scientificName"])
	}else{
		if (taxa.i%in%c("Chordata","Arthropoda")) RangeFragments=subset(RangeFragments.ena.inv,names(RangeFragments.ena.inv)%in%splist[splist$phylum==taxa.i,"scientificName"])
		if (taxa.i%in%"Other animal") RangeFragments=subset(RangeFragments.ena.inv,names(RangeFragments.ena.inv)%in%splist[splist$kingdom=="Animalia"&(!splist$phylum%in%c("Chordata","Arthropoda")),"scientificName"])
	}

	fragmentTable1=data.frame(fragBinomial=unique(names(RangeFragments)),fragSpecies=1:length(unique(names(RangeFragments))))
	fragmentTable2=data.frame(fragNum=1:length(RangeFragments),fragSize=unlist(lapply(RangeFragments,length)),fragBinomial=names(RangeFragments))
	AlienfragmentTable=fragmentTable2%>%left_join(fragmentTable1,by="fragBinomial")
	print(taxa.i)

	mycl <- makePSOCKcluster(no_cores);	
	parSapply(cl=mycl,X=1:length(toload),get.sdm.predict,toload,AlienfragmentTable,RangeFragments,g1cc,taxa.i)  
	stopCluster(mycl)
	gc()
}

########################################################################################  #####
# generate matrix of invasion probabilities for each species in each cell using best ER model of each species
########################################################################################  
library(dplyr)
load("data analysis/Domain.ena.rda")
load("data analysis/AlienRangeFragments.alltaxa.ENA.RData")
splist=read.csv("data analysis/US-RIIS.csv")%>%select(kingdom,phylum,class,scientificName)%>%filter(scientificName%in%names(RangeFragments.ena.inv))
taxa=c("Chordata","Arthropoda","Other animal","Plantae")
inputER=ifelse(taxa.i=="Arthropoda",0.84,ifelse(taxa.i=="Chordata",0.78,ifelse(taxa.i=="Other animal",0.91,0.87)))
toload=list.files("data analysis/SDM.bestER/")%>%.[grep(taxa.i,.)]%>%gsub(".rda","",.)%>%gsub(paste("PredictedER",taxa.i,"ACCESS-CM2_",sep="."),"",.)
# lab=c("ssp126_2040","ssp126_2060","ssp126_2080",
	# "ssp245_2040","ssp245_2060","ssp245_2080",
	# "ssp370_2040","ssp370_2060","ssp370_2080",
	# "ssp585_2040","ssp585_2060","ssp585_2080")
get.invasion.prob=function(i,domain.ena,AlienfragmentTable,fragmentTable,simRanges,nRandomStarts=100){
 focRowAlien<-which(AlienfragmentTable$fragSpecies==i)
 nFrags<-length(focRowAlien)
 InvasionProbs<-matrix(0,ncol=nrow(domain.ena),nrow=nFrags)
 for(x in 1:nFrags){
    if(fragmentTable$RandomStart[focRowAlien[x]]==1){
      fragSize<-AlienfragmentTable$fragSize[focRowAlien[x]]
      InvasionProbs[x,as.numeric(names(simRanges[[focRowAlien[x]]]))]<-as.numeric(simRanges[[focRowAlien[x]]])/min(nRandomStarts,fragSize)
    }else{
      InvasionProbs[x,simRanges[[focRowAlien[x]]]]<-1
    }
  }
  InvasionProbVec<-1-matrixStats::colProds(1-InvasionProbs)
  InvasionProbability<-InvasionProbVec   
  return(InvasionProbability)
}
		
if(taxa.i=="Plantae"){
		RangeFragments=subset(RangeFragments.ena.inv,names(RangeFragments.ena.inv)%in%splist[splist$kingdom=="Plantae","scientificName"])
	}else{
		if (taxa.i%in%c("Chordata","Arthropoda")) RangeFragments=subset(RangeFragments.ena.inv,names(RangeFragments.ena.inv)%in%splist[splist$phylum==taxa.i,"scientificName"])
		if (taxa.i%in%"Other animal") RangeFragments=subset(RangeFragments.ena.inv,names(RangeFragments.ena.inv)%in%splist[splist$kingdom=="Animalia"&(!splist$phylum%in%c("Chordata","Arthropoda")),"scientificName"])
	}
fragmentTable1=data.frame(fragBinomial=unique(names(RangeFragments)),fragSpecies=1:length(unique(names(RangeFragments))))
fragmentTable2=data.frame(fragNum=1:length(RangeFragments),fragSize=unlist(lapply(RangeFragments,length)),fragBinomial=names(RangeFragments))
AlienfragmentTable=fragmentTable2%>%left_join(fragmentTable1,by="fragBinomial")
nSpecies<-length(unique(AlienfragmentTable$fragSpecies))

# Current alien richness and invasion probs
load(paste0("data analysis/simulation.ena.predict/",taxa.i,"/PredictedRangeSize_",inputER,".phylo2.rda"))   #simRanges; fragmentTable
InvasionProbabilityMatrix<-do.call(rbind,lapply(1:nSpecies,get.invasion.prob,domain.ena,AlienfragmentTable,fragmentTable,simRanges))
InvasionProbs=colSums(InvasionProbabilityMatrix);names(InvasionProbs)=1:nrow(domain.ena)
InvasionRichness<-as.numeric(table(unlist(RangeFragments)));names(InvasionRichness)=1:nrow(domain.ena)

mf=function(i,RangeFragments) {return(data.frame(species=names(RangeFragments)[i],domain.ID=RangeFragments[[i]]))}
rangeData=do.call(rbind,lapply(1:length(RangeFragments),mf,RangeFragments))
rangeData$val=1;
CurrentRichnessMatrix=rangeData%>% tidyr::pivot_wider(names_from = species, values_from = val,values_fill = 0) %>% left_join(data.frame(domain.ID=domain.ena$domain.ID),.,by="domain.ID")%>% 
	tibble::column_to_rownames("domain.ID")%>% t() %>%.[unique(AlienfragmentTable$fragBinomial),]
CurrentRichnessMatrix[is.na(CurrentRichnessMatrix)]=0

#Future invasion prob
invprob=c();range.stat.all=c()
for (name in 1:length(toload)){
	load(paste0("data analysis/SDM.bestER/",paste("PredictedER",taxa.i,"ACCESS-CM2_",sep="."),toload[name],".rda"))
	InvasionProbabilityMatrix.future<-matrix(0,ncol=nrow(domain.ena),nrow=nSpecies)
	range.stat=c()
	for(i in 1:nSpecies){ # loop over each species  
	  focRowAlien<-which(AlienfragmentTable$fragSpecies==i)
	  nFrags<-length(focRowAlien)
	  # for each alien range fragment calculate the probability of invading each cell
	  InvasionProbs<-matrix(0,ncol=nrow(domain.ena),nrow=nFrags)
	  for(x in 1:nFrags){
		if(fragmentTable$RandomStart[focRowAlien[x]]==1){
		  fragSize<-AlienfragmentTable$fragSize[focRowAlien[x]]
		  InvasionProbs[x,as.numeric(names(simRanges[[focRowAlien[x]]]))]<-as.numeric(simRanges[[focRowAlien[x]]])/min(nRandomStarts,fragSize)
		}else{
		  InvasionProbs[x,simRanges[[focRowAlien[x]]]]<-1
		}
	  }
	  
	  # the overall probability of an alien species invading a cell (a cell can be invaded from multiple fragments)
	  InvasionProbVec<-1-matrixStats::colProds(1-InvasionProbs)
	  
	  # # set invasion probabilities of cells that are in the species' native range to zero
	  # focRowNative<-which(NativefragmentTable$fragSpecies==i)
	  # focNativeRange<-unique(unlist(NativeRangeFragments[focRowNative]))
	  # InvasionProbVec[focNativeRange]<-0
	  InvasionProbabilityMatrix.future[i,]<-InvasionProbVec

	  range.diff=InvasionProbabilityMatrix.future[i,]-CurrentRichnessMatrix[i,]
	  range.diff.prob=InvasionProbabilityMatrix.future[i,]-InvasionProbabilityMatrix[i,]
	  tmp=rbind(data.frame(type="range.shift",value=(sum(InvasionProbabilityMatrix.future[i,])-sum(CurrentRichnessMatrix[i,]))/sum(CurrentRichnessMatrix[i,])*100),	  	  
		data.frame(type="range.gain",value=sum(range.diff[range.diff>0])/sum(CurrentRichnessMatrix[i,])*100),
	    data.frame(type="range.loss",value=-sum(range.diff[range.diff<0])/sum(CurrentRichnessMatrix[i,])*100),
		data.frame(type="range.shift.prob",value=(sum(InvasionProbabilityMatrix.future[i,])-sum(InvasionProbabilityMatrix[i,]))/sum(InvasionProbabilityMatrix[i,])*100),
		data.frame(type="range.gain.prob",value=sum(range.diff.prob[range.diff.prob>0])/sum(InvasionProbabilityMatrix[i,])*100),
	    data.frame(type="range.loss.prob",value=-sum(range.diff.prob[range.diff.prob<0])/sum(InvasionProbabilityMatrix[i,])*100))
	  tmp$species=unique(AlienfragmentTable$fragBinomial)[i]
	  range.stat=rbind(range.stat,tmp)	  
	}
	rownames(InvasionProbabilityMatrix.future)=unique(AlienfragmentTable$fragBinomial)
	colnames(InvasionProbabilityMatrix.future)=domain.ena$domain.ID
	range.stat$scenerio=lab[name]
	invprob=cbind(invprob,colSums(InvasionProbabilityMatrix.future))
	range.stat.all=rbind(range.stat.all,range.stat)
	#print(name)
}
colnames(invprob)=lab
invprob2=cbind(invprob,current=colSums(InvasionProbabilityMatrix))

AlienRichness<-table(unlist(RangeFragments))
domain.ena$CurrentRichness<-0
domain.ena$CurrentRichness[match(as.numeric(names(AlienRichness)),domain.ena$domain.ID)]<-as.numeric(AlienRichness)
domain=cbind(domain.ena,invprob2)
domain[is.na(domain)]=0
domain.fl=domain%>%filter(x>=-88&x<=-79.5&y>=24&y<=31)

library(ggpubr)
#boxplot1: range expansion for each species compared to current
ggplot(range.stat.all%>%filter(grepl("shift.prob",type)), aes(x=scenerio, y=value)) + 
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
	theme_bw()+ coord_flip()

range.stat2=range.stat.all%>%filter(grepl("shift.prob",type))%>%group_by(type,scenerio)%>%summarize(mean=mean(value),se=plotrix::std.error(value),median=median(value))
ggplot(range.stat2, aes(x=scenerio, y=mean)) + 	
	geom_errorbar(aes(ymin=mean-se, ymax=mean+se),size=1,width=0.2,show.legend=FALSE,alpha=0.6)+
	geom_point(size=4,shape=21,color="black",fill="#619CFF") +coord_flip() +
	geom_vline(xintercept=3.5,col="#F8766D",size=1,linetype="longdash")+
	geom_vline(xintercept=6.5,col="#F8766D",size=1,linetype="longdash")+
	geom_vline(xintercept=9.5,col="#F8766D",size=1,linetype="longdash")+
	theme_bw()+
	labs(y="Predicted range shift \n(Future - present, %)",x="Climate change scenerios")#+ facet_wrap(~type, nrow = 1,scales="free_x")
	
#boxplot2
p=c()
for (i in lab){
tmp=data.frame(scenerio=i,diff=domain.fl[,i]-domain.fl[,"CurrentRichness"],diff.predic=domain.fl[,i]-domain.fl[,"current"])
p=rbind(p,tmp)
}
ggplot(p, aes(x=scenerio, y=diff)) + 
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))

#map
data.p=domain.ena[,c("domain.ID","x","y")];rownames(data.p)=data.p$domain.ID
data.p=cbind(data.p,InvasionProbs,InvasionRichness)
ma=ceiling(max(apply(data.p[,-c(1:3)],2,max,na.rm=T)))

the=theme_bw()+theme(legend.background=element_rect(fill='transparent'),
		legend.text=element_text(size=10),
		legend.title=element_text(face="bold",size=10))

p=list();ii=0
for (i in toload[1:12]){
	ii=ii+1
	pdata=domain[,c("x","y",i,"current","CurrentRichness")]
	colnames(pdata)[3]="value"
	pdata=pdata%>%filter(value>1&current>1)
	pdata$value=(pdata$value-pdata$current)#/pdata$current
	p[[ii]]=ggplot() + geom_tile(data = pdata, aes(x = x, y = y, fill = value)) +
	  coord_quickmap() +  labs(x="",y="",fill=lab[ii])+
	  scale_fill_gradientn(colours = viridis::viridis(99),na.value = "lightgray")+the#+theme(legend.position=c(0.8,0.25))	
}  
ggarrange(p[[1]],p[[2]],p[[3]],p[[4]],
	p[[5]],p[[6]],p[[7]],p[[8]],
	p[[9]],p[[10]],p[[11]],p[[12]],
	nrow=4,ncol=3,labels="auto",font.label = list(size = 20))	

#plot richness
p1=ggplot() +  geom_tile(data = data.p, aes(x = x, y = y, fill = InvasionRichness)) +
  coord_quickmap() +  labs(x="",y="",fill="InvasionRichness")+
  scale_fill_gradientn(limits=c(0,ma),colours = viridis::viridis(99),na.value = "lightgray")+the
p2=ggplot() +  geom_tile(data = data.p, aes(x = x, y = y, fill = InvasionProbs)) +
  coord_quickmap() +  labs(x="",y="",fill="InvasionProbs")+
  scale_fill_gradientn(limits=c(0,ma),colours = viridis::viridis(99),na.value = "lightgray")+the
ggarrange(p1,p2,nrow=1,ncol=2,labels="auto",font.label = list(size = 20))

########################################################################################  #####
# shifts in climate vs. shifts in invasion riskcvs. shifts in native species richness ###############################
########################################################################################  
	#plot patterns of similarity
	
	
###code graveyard ----	

##phylo and trait similarity
source("data analysis/read_Phylogenetic_Tree.R")
library(ape);library(dplyr)	
tree0=get(load("data analysis/Phylo.us.Rdata"))
rangeData=spdis%>%filter(Tip_Label%in%tree0$tip)
tree=keep.tip(tree0,tip=unique(rangeData$Tip_Label))

d=rangeData[,1:2];d$val=1;rownames(d)=NULL
d=d %>% tidyr::pivot_wider(names_from = Tip_Label, values_from = val,values_fill = 0) %>% tibble::column_to_rownames("domain.ID")
d <- as.matrix(d)
for (i in 1:ncol(d)) d[,i] <- as.logical(d[,i])
ancs <- ancestor(tree=tree, tip=tree$tip.label, include.tip = TRUE)

ii <- 1;
for (i in 1:(nrow(d)-1)) {
	for (j in (i+1):nrow(d)) {
		comm1 <- colnames(d)[which(as.logical(d[i,]))]
		comm2 <- colnames(d)[which(as.logical(d[j,]))]
		comm <- list(comm1, comm2);names(comm)=paste(i,j,sep="_")
		all.p <- phyBeta(comm,tree,ancs,method="sim2")# modified simpson like taxonmic beta
		pBeta.mat[i,j]=all.p[[1]]
		pBeta.mat[j,i]=all.p[[2]]
		ii <- ii + 1
		}
	print(i)
	}
save(pBeta.mat, file="data analysis/CommSimilarity.ENA.phylo2.rda")

## functional similarity
## Measure functional beta
library(dplyr)	
load("data analysis/traits.ENA.RData")
source("data analysis/MVNH_functions.R") 
rangeData=spdis%>%filter(Tip_Label%in%rownames(trait_mat3))
d=rangeData[,1:2];d$val=1;rownames(d)=NULL
d=d %>% tidyr::pivot_wider(names_from = Tip_Label, values_from = val,values_fill = 0) %>% tibble::column_to_rownames("domain.ID")
for (i in 1:ncol(d)) d[,i] <- as.logical(d[,i])
d <- as.matrix(d)
Beta.turnover <- matrix(NA, nrow=nrow(d), ncol=nrow(d));rownames(Beta.turnover) <- colnames(Beta.turnover) <- rownames(d)

ii <- 1;
for (i in 1:(nrow(d)-1)) {
	for (j in (i+1):nrow(d)) {
		comm1 <- colnames(d)[which(as.logical(d[i,]))]
		comm2 <- colnames(d)[which(as.logical(d[j,]))]
		db1=subset(trait_mat3,rownames(trait_mat3)%in%comm1)
		db2=subset(trait_mat3,rownames(trait_mat3)%in%comm2)
		re=MVNH_dissimilarity(db1,db2,var.names = colnames(trait_mat3))
		Beta.turnover[i,j] <- re$M[1]#Mahalanobis distance between the two niche centroids
		ii <- ii + 1
		}
	print(i)
	}

mf=function(pBeta.mat){
	pBeta.mat[is.na(pBeta.mat)]=0
	psim.t=pBeta.mat+t(pBeta.mat)
	diag(psim.t)=NA
	return(psim.t)
}
Beta.turnover=mf(Beta.turnover);Beta.nest=mf(Beta.nest)
save(Beta.turnover, file="data analysis/CommSimilarity.ENA.trait.turnover.rda")

## simulation
library(data.table);library(dplyr)
load("data analysis/Domain.ena.rda")
EnvironmentData=domain.ena%>%select(domain.ID)
toload=c(list.files("data analysis/SDM.sp.simm/"),list.files("data analysis/SDM.clim.simm/"))%>%.[grep("CommSimilarity.",.)]%>%gsub(".rda","",.)
RangeFragments<-get(load("data analysis/AlienRangeFragments.plants.ENA.RData"))
mf.clim=function(i,clim) i[i%in%clim]
fragmentTable1=data.frame(fragBinomial=unique(names(RangeFragments)),fragSpecies=1:length(unique(names(RangeFragments))))
fragmentTable2=data.frame(fragNum=1:length(RangeFragments),fragSize=unlist(lapply(RangeFragments,length)),fragBinomial=names(RangeFragments))
fragmentTable=fragmentTable2%>%left_join(fragmentTable1,by="fragBinomial")
mf=function(i,RangeFragments) {return(data.frame(Tip_Label=names(RangeFragments)[i],domain.ID=RangeFragments[[i]],fragNum=i))}
IntroductionsFiltered=do.call(rbind,lapply(1:length(RangeFragments),mf,RangeFragments))

get.sim=function(REP,fragmentTable,RangeFragments,EnvironmentData,IntroductionsFiltered,adjCells.ena,CommSimilarity,mu,cene,mat.name,inputProbSurface=TRUE){
	source("Functions.R")
	#source("Code/2.Functions/Functions.R")
	simResults<-SimulateSpread(Range="Alien",fragmentTable=fragmentTable,RangeFragments=RangeFragments,
									 EnvironmentData=EnvironmentData,IntroductionsFiltered=IntroductionsFiltered,adjCells=adjCells.ena,
									 SpeciesOptimum=NULL,CommSimilarity=CommSimilarity,plotIt=FALSE,envVar=c(),mu=mu,
									 seed="Random",niche="ClimateMatching",inputProbSurface=TRUE)
	names(simResults)=names(RangeFragments)		
    outfile<-paste0("data analysis/SDM.",cene,"/ER_sdm_",mat.name,"_",mu,"_",REP,".rda")
	save(simResults,file=outfile)
}

library(parallel)
no_cores <- 50#detectCores() - 1
### simulate environmental resistance model with different mu values
for (mat in 1:length(toload)){
	mycl <- makePSOCKcluster(no_cores);
	mat.name=toload[mat]
	if (grepl("clim.",mat.name)==0) {
		simmat=get(load(paste0("data analysis/SDM.sp.simm/",mat.name,".rda")))	
		simmat=1-simmat;
		mu=3	
		cene=strsplit(mat.name,"_")[[1]][2]
			
	}else{
		simmat=get(load(paste0("data analysis/SDM.clim.simm/",mat.name,".rda")))
		mu=10	
		cene=strsplit(mat.name,"_")[[1]][1]%>%gsub("CommSimilarity.clim.","",.)	
		cene=paste0("clim.",cene)		
	}
	print(mat.name)
	parSapply(cl=mycl,X=1:200,get.sim,fragmentTable,RangeFragments,EnvironmentData=EnvironmentData,IntroductionsFiltered,adjCells.ena,
		CommSimilarity=simmat,mu=mu,cene=cene,mat.name)  
	stopCluster(mycl)
	rm(simmat);gc()
  }  
		
# plot observed range
	windows()
	library(scales)
	tmp=get(load("data analysis/simulation/Alien_Random_ClimateResistance_ele_30_17.rda"))
	currspeciesname=unique(names(RangeFragments))[45]
	plot(EnvironmentData$x,EnvironmentData$y,cex=0.5)
	focalDom<-EnvironmentData%>%filter(domain.ID%in%unlist(subset(RangeFragments,names(RangeFragments) %in% currspeciesname)))
    points(focalDom$x,focalDom$y,col="green",pch=16,cex=0.5)  
	focalDom2<-EnvironmentData%>%filter(domain.ID%in%unlist(subset(tmp,names(tmp) %in% currspeciesname)))
    points(focalDom2$x,focalDom2$y,col=alpha("red", 0.5),pch=16,cex=0.5)
