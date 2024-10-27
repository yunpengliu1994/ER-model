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
####get distributions ################################
##############################################################
## step 1.2 --- data downloading of invasive species
library(dplyr)
#splist=read.csv("data analysis/horizon_scan_100_plant_list.csv")%>%filter(!Dis%in%"US native")
#Global Register of Introduced and Invasive Species - United States (Contiguous) (ver.2.0, 2022)
splist=read.csv("data analysis/US-RIIS.csv")%>%.[grep("invasive",.[,"degreeOfEstablishment"]),]%>%filter(kingdom%in%c("Animalia","Plantae"))#invasive species lists obtained from ISSG
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
ele.us <- raster("data analysis/wc2.1_10m_bio/wc2.1_10m_elev.tif")#raster map downloading from worldclim.
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
	if (us) dat.all=dat.all%>%filter(longitude>=-124.7258&(longitude<=-66.94989&latitude>=24.49813 & latitude<=49.38436))%>%na.omit()
	
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
#########get distribution data from BONAP ####################
##############################################################
#load map and intersect with grid cells
library(raster);library(dplyr);library(sf);
load("data analysis/Domain.ena.rda")
ele.us <- raster("data analysis/wc2.1_10m_bio/wc2.1_10m_elev.tif")#%>%aggregate(fact=3)
ENA.name=c("ME","NH","VT","NY","MA","RI","CT","NJ","PA","DE","MD","MI","OH","IN","IL","WI","WV","VA","NC","TN","KY","SC","GA","AL","MS","FL")
ENA.state=read_sf('data analysis/cb_2014_us_county_20m/cb_2014_us_state_20m.shp')%>%filter(STUSPS%in%ENA.name)%>%as.data.frame%>%.[,"STATEFP"]
ENA <- read_sf('data analysis/cb_2014_us_county_20m/cb_2014_us_county_20m.shp')%>%filter(!STATEFP%in%ENA.state)%>%
	st_make_valid()%>%  mutate(area.county = st_area(.) %>% as.numeric())%>% select("GEOID","NAME","area.county")
grids=rasterFromXYZ(domain.ena[,c("x","y","domain.ID")],res=res(ele.us),crs=st_crs(ele.us))%>%stars::st_as_stars()%>%st_as_sf()%>% 
	st_set_crs(st_crs(ENA))%>%mutate(area.grids = st_area(.) %>% as.numeric())
int <- st_intersection(grids, ENA)%>%mutate(area.int = st_area(.) %>% as.numeric(),area.por = area.int/area.grids*100,area.por2= area.int/area.county*100) 
#save(int,file="data analysis/cb_2014_us_county_20m/intersec_with_10min_grids.rda")

load("data analysis/cb_2014_us_county_20m/intersec_with_10min_grids.rda")
int2=int%>%filter(area.por>=25)
length(unique(int2$domain.ID))
miss.grid=domain.ena%>%filter(!domain.ID%in%unique(int2$domain.ID))
int3=int%>%filter(domain.ID%in%miss.grid$domain.ID)#these two county are within a grid
int.fin=rbind(int2,int3)%>%as.data.frame%>%dplyr::select(domain.ID,GEOID)%>%distinct
int.fin$GEOID=paste0("X",int.fin$GEOID)
spdis=data.table::fread("data analysis/BONAP data/county_data.csv")%>%filter(V1%in%int.fin$GEOID)%>%as.data.frame
nodis=colSums(spdis[,-1])
spdis=spdis[,c("V1",names(nodis)[nodis>0])]
rownames(spdis)=spdis[,1];spdis=spdis[,-1]
dis.bonap=tidyr::pivot_longer(spdis,everything(),names_to = "Species", values_to = "dis")%>%arrange(Species) 
dis.bonap$GEOID=rep(rownames(spdis),ncol(spdis))
dis.bonap=dis.bonap%>%filter(dis>0)%>%left_join(int.fin,by="GEOID",relationship =  "many-to-many")%>%dplyr::select(Species,domain.ID)%>%distinct
dis.bonap$Species=gsub("_", " ",dis.bonap$Species)
load("data analysis/introducedData.ENA.Rdata")
RangeFragments<-get(load("data analysis/NativeRangeFragments.plants.ENA.RData"))
mf=function(i,RangeFragments) {return(data.frame(Species=names(RangeFragments)[i],domain.ID=RangeFragments[[i]]))}
rangeData=do.call(rbind,lapply(1:length(RangeFragments),mf,RangeFragments))%>%rbind(dis.bonap)%>%distinct%>%
	left_join(introduce.us,by=c("Species"="taxon_name","domain.ID"))

#correct species names
trenames=rangeData%>%dplyr::select(Species)%>%distinct()
library(TNRS)	
res=TNRS(taxonomic_names = trenames$Species,sources = "wfo")%>%filter((!Taxonomic_status%in%"No opinion")&(!Name_matched_rank%in%"genus")&Genus_score==1&Overall_score>=0.9)%>%
	dplyr::select(Name_submitted,Taxonomic_status,Accepted_name_id,Accepted_name,Source)
rangeData=rangeData%>%left_join(res,by=c("Species"="Name_submitted"))
rangeData[is.na(rangeData$Accepted_name),"Accepted_name"]=rangeData[is.na(rangeData$Accepted_name),"Species"]

## intergrate invasive species records from pointed data and from new gridded data based on BONAP,POWO, etc.
## then merge fragsments which distance lower than 18 cells
RangeFragments<-get(load("data analysis/AlienRangeFragments.alltaxa.ENA.RData"))
sp=unique(names(RangeFragments))
res=TNRS(taxonomic_names = sp,sources = "wfo")%>%filter((!Taxonomic_status%in%"No opinion")&(!Name_matched_rank%in%"genus")&Genus_score==1&Overall_score>=0.9)%>%
	dplyr::select(Name_submitted,Taxonomic_status,Accepted_name_id,Accepted_name)
namecheck=data.frame(name=names(RangeFragments),Acc=res[match(names(RangeFragments),res$Name_submitted),"Accepted_name"])
namecheck[is.na(namecheck$Acc),"Acc"]=namecheck[is.na(namecheck$Acc),"name"]
names(RangeFragments)=namecheck$Acc

mf=function(i,RangeFragments) {return(data.frame(Species=names(RangeFragments)[i],domain.ID=RangeFragments[[i]],fragNum=i))}
rangeData.inv=do.call(rbind,lapply(1:length(RangeFragments),mf,RangeFragments))%>%.[!duplicated(.[,c('Species','domain.ID')]),]
newdata=rangeData%>%filter(Accepted_name%in%names(RangeFragments))%>%left_join(rangeData.inv,by=c("Accepted_name"="Species","domain.ID"))%>%
	filter(is.na(fragNum))%>%select(Accepted_name,domain.ID)%>%distinct()
# newfrags=setNames(split(newdata$domain.ID, seq(nrow(newdata))),newdata$Species)
# frags=c(RangeFragments,newfrags)
sp.new=unique(newdata$Accepted_name)
domain=domain.ena%>%select(x,y,domain.ID,row,col)
arrange.frag=function(i,RangeFragments,newdata,domain){
	require(dplyr)
#for (i in sp.new){
	#i=sp.new[1]
	cell.i=newdata%>%filter(Accepted_name==i)%>%.[,"domain.ID"]
	range.i=subset(RangeFragments,names(RangeFragments)==i)
	# add newdata to the exiting fragments
	get.range.new=function(range.i.j,domain,cell.i){
		cell.limit=domain%>%filter(domain.ID%in%range.i.j)%>%summarize(row.max=max(row)+9,row.min=min(row)-9,col.max=max(col)+9,col.min=min(col)-9)
		cell.range=domain%>%filter(row>=cell.limit$row.min&row<=cell.limit$row.max&col>=cell.limit$col.min&col<=cell.limit$col.max)
		range.i.j.new=c(range.i.j,cell.i[cell.i%in%cell.range$domain.ID])
		return(range.i.j.new)
	}
	range.i.new=lapply(range.i,get.range.new,domain,cell.i)

	# adding new fragments
	cell.other=cell.i[!cell.i%in%unlist(range.i.new)]
	if (length(cell.other)==0)  input = range.i.new else{
		if(length(cell.other)==1) input = c(list(cell.other),range.i.new) else{
			mf=function(n,domain){
				cell.limit=domain%>%filter(domain.ID%in%n)%>%summarize(row.max=max(row)+9,row.min=min(row)-9,col.max=max(col)+9,col.min=min(col)-9)
				cell.range=domain%>%filter(row>=cell.limit$row.min&row<=cell.limit$row.max&col>=cell.limit$col.min&col<=cell.limit$col.max)
				return(cell.range$domain.ID)
			}
			cell.other.range=lapply(cell.other,mf,domain)
			names(cell.other.range)=cell.other
			repeat {
				 tbl <- table(unlist(cell.other.range))
				 if (length(tbl[tbl > 1]) == 0) { break }
				 idx <- which(sapply(seq_len(length(cell.other.range)), function(i) {
					any(cell.other.range[[i]] == names(tbl[tbl > 1])[1])
				  }))
				  newvec <- list(sort(unique(unlist(cell.other.range[idx]))))
				  names(newvec)=paste(names(cell.other.range)[idx],collapse="_")
				  cell.other.range <- c(cell.other.range, newvec)[-idx]
				}
			new.frags=lapply(names(cell.other.range),function(x) stringr::str_split(x,"_")%>%unlist);
			input=c(new.frags,range.i.new)
		}
	}

	# Merging Listed Vectors that share Elementsand then make identify connected sets of cells (i.e., fragments)
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
	names(input)=rep(i,length(input))
	return(input)
}

library(parallel)
no_cores <- detectCores()-1
mycl <- makePSOCKcluster(no_cores);
RangeFragments.new=do.call(c,parLapply(cl=mycl,X=sp.new,arrange.frag,RangeFragments,newdata,domain))

other.taxa = subset(RangeFragments,!names(RangeFragments)%in%sp.new)
sp.other=unique(names(other.taxa))
reduce.frag=function(i,other.taxa,domain){
	require(dplyr)
	#i=sp.other[1]
	range.i=subset(other.taxa,names(other.taxa)==i)	
	#range.i[[3]]=c(1,2,3,range.i[[1]][1:3])
	mf=function(n,domain){
		cell.limit=domain%>%filter(domain.ID%in%n)%>%summarize(row.max=max(row)+9,row.min=min(row)-9,col.max=max(col)+9,col.min=min(col)-9)
		cell.range=domain%>%filter(row>=cell.limit$row.min&row<=cell.limit$row.max&col>=cell.limit$col.min&col<=cell.limit$col.max)
		return(cell.range$domain.ID)
	}
	range.i.extent=lapply(range.i,mf,domain)
	names(range.i.extent)=1:length(range.i.extent)
	repeat {
		tbl <- table(unlist(range.i.extent))
		if (length(tbl[tbl > 1]) == 0) { break }
		idx <- which(sapply(seq_len(length(range.i.extent)), function(i) {
			any(range.i.extent[[i]] == names(tbl[tbl > 1])[1])
		}))
		newvec <- list(sort(unique(unlist(range.i[idx]))))
		#names(newvec)=paste(names(range.i.extent)[idx],collapse="_")		
		range.i <- c(range.i, newvec)[-idx]
		range.i.extent <- range.i.extent[-idx]
	}
	names(range.i)=rep(i,length(range.i))
	return(range.i)
}
other.taxa.new=do.call(c,parLapply(cl=mycl,X=sp.other,reduce.frag,other.taxa,domain))

RangeFragments.v2=c(RangeFragments.new,other.taxa.new)
stopCluster(mycl)
save(RangeFragments.v2,file="data analysis/AlienRangeFragments.alltaxa.ENA.new.RData")

## combine introduced information in BONAP and removed non-native records
introduced.sp=read.csv("data analysis/BONAP data/North American Species List 2014 Nativity.csv")
res=TNRS(taxonomic_names = introduced.sp$AcceptedName,sources = "wfo")%>%filter((!Taxonomic_status%in%"No opinion")&(!Name_matched_rank%in%"genus")&Genus_score==1&Overall_score>=0.9)%>%
	dplyr::select(Name_submitted,Taxonomic_status,Accepted_name_id,Accepted_name)
introduced.sp2=	introduced.sp%>%left_join(res,by=c("AcceptedName"="Name_submitted"))
introduced.sp2[is.na(introduced.sp2$Accepted_name),"Accepted_name"]=introduced.sp2[is.na(introduced.sp2$Accepted_name),"AcceptedName"]
introduced.sp2=introduced.sp2%>%dplyr::select(Nativity.status,Accepted_name)%>%distinct
inv.list=unique(c(names(RangeFragments.v2),sp))
rangeData.native=rangeData%>%filter((!Accepted_name%in%inv.list)&(!is.na(Taxonomic_status)))%>%left_join(introduced.sp2,by="Accepted_name",relationship = "many-to-many")%>%
	filter(is.na(introduced)&(is.na(Nativity.status)|Nativity.status!="I"))%>%dplyr::select(domain.ID,Accepted_name,Accepted_name_id,Source)%>%distinct
save(rangeData.native,file="data analysis/NativeRangeFragments.plants.ENA.clean.new.RData")#8524 spp.

##############################################################
#########overlay distribution data with BONAP ################
##############################################################
library(dplyr);library(sf);
load("data analysis/Domain.ena.rda")
load("data analysis/cb_2014_us_county_20m/intersec_with_10min_grids.rda")
load("data analysis/AlienRangeFragments.alltaxa.ENA.new.RData")
load("data analysis/NativeRangeFragments.plants.ENA.clean.new.RData")
int2=int%>%filter(area.por>=25)
length(unique(int2$domain.ID))
miss.grid=domain.ena%>%filter(!domain.ID%in%unique(int2$domain.ID))
int3=int%>%filter(domain.ID%in%miss.grid$domain.ID)#these two county are within a grid
int.fin=rbind(int2,int3)%>%as.data.frame%>%dplyr::select(domain.ID,GEOID)%>%distinct
int.fin$GEOID=paste0("X",int.fin$GEOID)

bonap=read.csv("data analysis/0BONAP county for 100 species 4-15-24 for Yunpeng Liu.csv")
bonap$FIPS=as.character(bonap$FIPS)
bonap=bonap%>%left_join(int.fin,by=c("FIPS"="GEOID"),relationship =  "many-to-many")%>%na.omit

sp.list=unique(bonap$Scientific.Name)#71 spp. occured in ena
rangeData.native=rangeData.native%>%filter(Accepted_name%in%sp.list)
inv.dis=subset(RangeFragments.v2,names(RangeFragments.v2)%in%sp.list)
mf=function(i,dis)data.frame(domain.ID=as.numeric(dis[[i]]),Accepted_name=names(dis)[i])
rangeData.inv=do.call(rbind,lapply(1:length(inv.dis),mf,inv.dis))
rangeData=rbind(rangeData.native,rangeData.inv)%>%left_join(domain.ena[,c("x","y","domain.ID")],by="domain.ID")
sp.list=unique(rangeData$Accepted_name)
bonap=bonap%>%filter(Scientific.Name%in%sp.list)%>%left_join(domain.ena,by="domain.ID")%>%select(domain.ID,Scientific.Name,x,y)%>%distinct
colnames(bonap)[2]="Accepted_name"
bonap$Source="BONAP"
rangeData$Source="Mydata"
dat.all=rbind(bonap,rangeData)
#Plot
library(ggplot2)
theme=theme_bw()+
		theme(axis.text = element_text(size=8,color='black'),
			axis.title = element_blank(),
			strip.text=element_text(size=12,color='black',face="italic"),
		strip.background = element_rect(fill = 'white', colour = 'black', linewidth = rel(2), linetype = 2),
		panel.background = element_rect(fill = '#619CFF', colour = 'black'))	
for(n in 1:6){
nstart=1+12*(n-1)
nend=ifelse(n==6,70,12+12*(n-1))
print(c(nstart,nend))
windows(width=12, height=9)
dat=dat.all%>%filter(Accepted_name%in%sp.list[nstart:nend])
ggplot(dat) + 
	geom_tile(aes(x = x, y = y),data = domain.ena,fill = "white")+
	geom_tile(aes(x = x, y = y),data = dat%>%filter(Source=="BONAP"),fill = "orange",alpha=0.5) +
	geom_tile(aes(x = x, y = y),data = dat%>%filter(Source=="Mydata"),fill = "blue",alpha=0.5) +
	facet_wrap(~Accepted_name)+	theme
	
}

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
# occur_points=get(load("data analysis/occ_data.plants/occ.us.plusBIEN.RData"))%>%filter(longitude>=-92.83333&longitude<=-66.83333&latitude>=24.16667 & latitude<=48.16667&(!Sci_name%in%""))
# splist1=sort(unique(occur_points$Sci_name))
# splist2=read.csv("data analysis/horizon_scan_100_plant_list.csv")%>%filter(!Dis%in%"US native")
# splist=unique(c(splist2$acc_name,splist1))
load("data analysis/NativeRangeFragments.plants.ENA.clean.new.RData")
splist=unique(rangeData.native$Accepted_name)
sp.in.ena=trp.sp.powo%>%filter(taxon_name%in%splist)
# seql=c(seq(1,nrow(sp.in.ena),9990),nrow(sp.in.ena))
# seql2=lapply(2:length(seql),function(i,seql){if (i<length(seql)){c(seql[i-1],(seql[i]-1))
	# } else {
	# c(seql[i-1],seql[i])
	# }},seql)
# for(i in 1:length(seql2)){
	# spname=paste(paste(sp.in.ena$AccSpeciesID[seql2[[1]][1]:seql2[[1]][2]],collapse=","),",",sep="")
	# write.table(spname,paste("data analysis/try.query.name",i,"txt",sep="."),quote=FALSE,row.names = FALSE,col.names = FALSE)
# }

#ref:Kattge, J, Bönisch, G, Díaz, S, et al. TRY plant trait database – enhanced coverage and open access. Glob Change Biol. 2020; 26: 119– 188. https://doi.org/10.1111/gcb.14904
require(data.table);library(dplyr)
TRYdata <- rbind(fread("data analysis/TRY data/30740.txt", header = T, sep = "\t", dec = ".", quote = "", data.table = T),
	fread("data analysis/TRY data/30421.txt", header = T, sep = "\t", dec = ".", quote = "", data.table = T),
	fread("data analysis/TRY data/32167.txt", header = T, sep = "\t", dec = ".", quote = "", data.table = T))%>%
	distinct	
TRYdata2=TRYdata%>%filter(!is.na(TraitID))%>%
	select(AccSpeciesID,TraitID,TraitName,DataName,OriglName,OrigValueStr,OrigUnitStr,StdValue,UnitName,Reference)%>%
	distinct()%>%left_join(sp.in.ena[,c("AccSpeciesID","taxon_name")],by="AccSpeciesID")%>%
	filter(!is.na(taxon_name))%>%distinct()
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
sp=data.frame(Spname=unique(alldata2$Spname))
spnames <- rWCVPdata::wcvp_names
splist.mat = sp %>% left_join(spnames,by=c("Spname"="taxon_name"),relationship ="many-to-many") %>% filter(taxon_status%in%c("Accepted","Synonym")) %>% 
	dplyr::select(Spname,plant_name_id,accepted_plant_name_id,taxon_rank,taxon_status)
sp.unmat=sp %>% filter(!Spname%in%unique(splist.mat$Spname)) 
write.csv(sp.unmat$Spname,"sp.unmat.csv")

res=read.csv("tnrs_result.csv")%>%filter((!Taxonomic_status%in%"No opinion")&(!Name_matched_rank%in%"genus")&Genus_score==1&Overall_score>=0.9)%>%
	dplyr::select(Name_submitted,Taxonomic_status,Accepted_name_id,Accepted_name)
sp.mat=splist.mat[order(splist.mat$taxon_status),]%>% .[!duplicated(.[,'Spname']),] %>% dplyr::select('Spname','taxon_status','accepted_plant_name_id')%>%
	left_join(spnames[,c('plant_name_id','taxon_name')],by=c('accepted_plant_name_id'='plant_name_id'))
colnames(sp.mat)=colnames(res)
spname.cor=rbind(sp.mat,res)%>%select('Name_submitted','Accepted_name')%>%distinct

alldata2=alldata2%>%left_join(spname.cor,by=c("Spname"='Name_submitted'))
colnames(alldata2)[4]="taxon_name"
traitdata=rbind(alldata2[,c("taxon_name","StdValue","trait")],alldata)%>%distinct()%>%na.omit()
traitdata$StdValue=as.numeric(traitdata$StdValue)
traitdata=traitdata%>%na.omit()%>%group_by(taxon_name,trait)%>%summarize(StdValue=mean(StdValue))
save(traitdata,file="data analysis/traits.all.RData")

#imputation missing values and filter ENA species
trait_mat0=traitdata %>% tidyr::pivot_wider(names_from = taxon_name, values_from = StdValue) %>% tibble::column_to_rownames("trait") %>%t()
load("data analysis/NativeRangeFragments.plants.ENA.clean.new.RData")
trait_mat0=subset(trait_mat0,rownames(trait_mat0)%in%unique(rangeData.native$Accepted_name))#7723 out of 9380 spp.
#imputation missing values
library(Rphylopars);library(ape);
tree0=get(load("data analysis/Phylo.us.imputate.Rdata"))
tree0$tip.label=gsub("_"," ",tree0$tip.label)
trait_mat=subset(trait_mat0,rownames(trait_mat0)%in%tree0$tip)[,-5]#for root dep. too many NAs (6172 out of 7723 spp.)
tree=keep.tip(tree0,tip=rownames(trait_mat))

trait_data=as.data.frame(trait_mat)
trait_data$species=rownames(trait_mat)
trait_data=trait_data[,c("species",colnames(trait_mat))]
trait_mat1=phylopars(trait_data ,tree,
    pheno_error = TRUE,phylo_correlated = TRUE,pheno_correlated = TRUE)
trait_mat2=	trait_mat1$anc_recon[rownames(trait_mat),]
trait_mat.nophy=subset(trait_mat0[,-5],!rownames(trait_mat0)%in%rownames(trait_mat2))
trait_mat3=scale(rbind(trait_mat2,trait_mat.nophy))
colnames(trait_mat3)=c("SLA","LN","WD","Hmax")
save(trait_mat3,file="data analysis/traits.ENA.RData")

################################################################################################################
# use native ranges to generate a community similarity matrix for calculating environmental resistance
################################################################################################################
## ER1: species assemblages ----
# generate a list where each element is a vector of the species in each grid cell
get.sim.sp=function(rangeData){
	speciesByCell<-list()
	speciesByCellsubset<-split(rangeData$Accepted_name,rangeData$domain.ID)
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
load("data analysis/NativeRangeFragments.plants.ENA.clean.new.RData")
CommSimilarity.sp=get.sim.sp(rangeData.native)
save(CommSimilarity.sp,file="/blue/matthewthomas1/yunpeng.liu/data analysis/CommSimilarity/ENA.sp.wfo.rda")

# #species with phylo infromation
# library(ape);library(dplyr)	
# #load("data analysis/Phylo.tipnames.POWO.Rdata")
# tree0=get(load("data analysis/Phylo.us.Rdata"))
# rangeData.phy=rangeData%>%filter(Tip_Label%in%tree0$tip)
# CommSimilarity.sp2=get.sim.sp(rangeData.phy)
# save(CommSimilarity.sp2,file=paste0("/blue/matthewthomas1/yunpeng.liu/data analysis/CommSimilarity/ENA.sp.rda"))

################################################################################################################
# use native ranges to generate a phylogenetic similarity matrix for calculating environmental resistance
################################################################################################################
## ER2: phylogenetic similarity ----
##using POWO to correct names
library(ape);library(data.table);library(dplyr)
tre=read.tree("data analysis/FromAo_Smith_100treesV2/ALLMB.tre")		#obtained:Nov 7,2023
trenames=data.frame(tip=sort(tre$tip))
trenames$Taxon_name=gsub("_", " ",trenames$tip)
spnames <- rWCVPdata::wcvp_names
splist.mat = trenames %>% left_join(spnames,by=c("Taxon_name"="taxon_name"),relationship = "many-to-many") %>% filter(taxon_status%in%c("Accepted","Synonym")) %>% 
	select(Taxon_name,accepted_plant_name_id,taxon_rank,taxon_status)

sp.unmat=trenames %>% filter(!Taxon_name%in%splist.mat$Taxon_name) 
library(TNRS)	
res=TNRS(taxonomic_names = sp.unmat$Taxon_name,sources = "wcvp")%>%filter((!Taxonomic_status%in%"No opinion")&(!Name_matched_rank%in%"genus")&Genus_score==1&Overall_score>=0.9)%>%
	dplyr::select(Name_submitted,Taxonomic_status,Accepted_name_id,Accepted_name)
sp.unmat2=sp.unmat%>% filter(!Taxon_name%in%res$Name_submitted) 
res2=TNRS(taxonomic_names = sp.unmat$Taxon_name,sources = "wfo")%>%filter((!Taxonomic_status%in%"No opinion")&(!Name_matched_rank%in%"genus")&Genus_score==1&Overall_score>=0.9)%>%
	dplyr::select(Name_submitted,Taxonomic_status,Accepted_name_id,Accepted_name)
	
sp.mat=splist.mat[order(splist.mat$taxon_status),] %>%dplyr::select('Taxon_name','taxon_status','accepted_plant_name_id')%>%
	left_join(spnames[,c('plant_name_id','taxon_name')],by=c('accepted_plant_name_id'='plant_name_id'))
colnames(sp.mat)=colnames(res)
spname.cor=rbind(sp.mat,res,res2)%>%select(-Accepted_name_id)%>%.[order(.[,"Taxonomic_status"]),]%>%.[!duplicated(.[,'Accepted_name']),]
trenames.powo=trenames%>%left_join(spname.cor,by=c("Taxon_name"="Name_submitted"))
trenames.powo[is.na(trenames.powo$Accepted_name),"Accepted_name"]=trenames.powo[is.na(trenames.powo$Accepted_name),"Taxon_name"]
tre$tip.label=trenames.powo[match(tre$tip.label,trenames.powo$tip),"Accepted_name"]

## imputate missing species
library(phytools);library(dplyr)	
#load("data analysis/Phylo.us.Rdata")#6894 matched
#load("data analysis/Phylo.us.park.Rdata")#3714 matched
load("data analysis/NativeRangeFragments.plants.ENA.clean.new.ori.RData")#correct names using wcvp
tre$tip.label=gsub(" ","_",tre$tip.label)
rangeData.native$Accepted_name=gsub(" ","_",rangeData.native$Accepted_name)

mf=function(i){re=ifelse (i[1]=="×",i[2],i[1]); return(re)}
tre.genus=data.frame(tip=tre$tip.label,genus=do.call(c,lapply(strsplit(tre$tip.label,split="_"),mf)))
sp.mat=rangeData.native%>%filter(Accepted_name%in%tre$tip)%>%select(Accepted_name)%>%distinct
sp.mat$genus=do.call(rbind,lapply(strsplit(sp.mat$Accepted_name,split="_"),mf))

sp.miss=rangeData.native%>%filter(!Accepted_name%in%tre$tip)%>%select(Accepted_name)%>%distinct
sp.miss$genus=do.call(rbind,lapply(strsplit(sp.miss$Accepted_name,split="_"),mf))
sp.add=sp.miss%>%filter(genus%in%tre.genus$genus)#1903 out of 2361 speices could be adding to the exiting genus
genus.add=sp.add$genus[!sp.add$genus%in%sp.mat$genus]%>%unique
keep.sp=tre.genus%>%filter(genus%in%genus.add|tip%in%sp.mat$Accepted_name)
tree=keep.tip(tre,tip=keep.sp$tip)%>%force.ultrametric()
for (i in 1:nrow(sp.add)) {tree=add.species.to.genus(tree,sp.add$Accepted_name[i], where="root");print(i)}
tree2=keep.tip(tree,tip=c(sp.mat$Accepted_name,sp.add$Accepted_name))
save(tree2,file="data analysis/Phylo.us.imputate.Rdata")
	
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


## caculate phylo simm
load("data analysis/Phylo.us.imputate.Rdata")
load("data analysis/NativeRangeFragments.plants.ENA.clean.new.RData")
tree2$tip.label=gsub("_"," ",tree2$tip.label)
rangeData=rangeData.native%>%filter(Accepted_name%in%tree2$tip.label)%>%distinct

d=rangeData[,c("Accepted_name","domain.ID")];d$val=1;rownames(d)=NULL
d=d %>% distinct %>% tidyr::pivot_wider(names_from = Accepted_name, values_from = val,values_fill = 0) %>% tibble::column_to_rownames("domain.ID")
#phylobeta=betapart::phylo.beta.pair(d, tree, index.family="sorensen")
## distribution data
for (i in 1:ncol(d)) d[,i] <- as.logical(d[,i])
d <- as.matrix(d)
pBeta.mat <- matrix(NA, nrow=nrow(d), ncol=nrow(d))
rownames(pBeta.mat) <- colnames(pBeta.mat) <- rownames(d)
ancs <- ancestor(tree=tree2, tip=tree2$tip.label, include.tip = TRUE)

for (i in 1:(nrow(d)-1)) {
	for (j in (i+1):nrow(d)) {
		comm1 <- colnames(d)[which(as.logical(d[i,]))]
		comm2 <- colnames(d)[which(as.logical(d[j,]))]
		comm <- list(comm1, comm2);names(comm)=paste(i,j,sep="_")
		all.p <- phyBeta(comm,tree2,ancs,method="sim2")# modified simpson like taxonmic beta
		pBeta.mat[i,j]=all.p[[1]]
		pBeta.mat[j,i]=all.p[[2]]
		}
	print(i)
	}
save(pBeta.mat, file="/blue/matthewthomas1/yunpeng.liu/data analysis/CommSimilarity/ENA.phylo.rda")

## functional similarity
## Measure functional beta
library(dplyr)	
load("data analysis/traits.ENA.RData")
source("data analysis/MVNH_functions.R") ## #this function comes from Lu et al. (2021) Methods in Ecology and Evolution, 12(10), 1953-1968. https://github.com/lvmuyang/MVNH
load("data analysis/NativeRangeFragments.plants.ENA.clean.new.RData")
rangeData=rangeData.native%>%filter(Accepted_name%in%rownames(trait_mat3))%>%distinct

d=rangeData[,1:2];d$val=1;rownames(d)=NULL
d=d %>% tidyr::pivot_wider(names_from = Accepted_name, values_from = val,values_fill = 0) %>% tibble::column_to_rownames("domain.ID")
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
		db1=subset(trait_mat3,rownames(trait_mat3)%in%comm1)%>%na.omit
		db2=subset(trait_mat3,rownames(trait_mat3)%in%comm2)%>%na.omit
		if (nrow(db1)>0&nrow(db2)>0){
		re=MVNH_dissimilarity(db1,db2,var.names = colnames(trait_mat3))#this function comes from Lu et al. (2021) Methods in Ecology and Evolution, 12(10), 1953-1968.
		Beta.turnover[i,j] <- re$M[1]#Mahalanobis distance between the two niche centroids
		Beta.nest[i,j] <- re$D[1]#determinant ratio component measuring the difference of the niche volumes		
		}else{
		Beta.turnover[i,j] <- NA
		Beta.nest[i,j] <- NA
		}
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
save(Beta.turnover, file="/blue/matthewthomas1/yunpeng.liu/data analysis/CommSimilarity/ENA.trait.turnover.rda")
save(Beta.traits, file="/blue/matthewthomas1/yunpeng.liu/data analysis/CommSimilarity/ENA.trait.rda")

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
human=human%>%ade4::dudi.pca(center = TRUE,scale = T,scannf = FALSE,nf = 4)
soil=envir.ena%>%dplyr::select(contains("_"))%>%ade4::dudi.pca(center = TRUE,scale = T,scannf = FALSE,nf = 5)
clim.soil=envir.ena%>%dplyr::select(ele,starts_with("bio"),contains("_"))%>%ade4::dudi.pca(center = TRUE,scale = T,scannf = FALSE,nf = 5)
	
clim.d=dist(clim$li)%>%as.matrix()%>%scales::rescale();diag(clim.d)=NA
soil.d=dist(soil$li)%>%as.matrix()%>%scales::rescale();diag(soil.d)=NA
human.d=dist(human$li)%>%as.matrix()%>%scales::rescale();diag(human.d)=NA
clim.soil.d=dist(clim.soil$li)%>%as.matrix()%>%scales::rescale();diag(clim.soil.d)=NA
clim.human.d=dist(cbind(clim$li,human$li))%>%as.matrix()%>%scales::rescale();diag(clim.human.d)=NA
clim.soil.human.d=dist(cbind(clim.soil$li,human$li))%>%as.matrix()%>%scales::rescale();diag(clim.soil.human.d)=NA
#p=rasterFromXYZ(cbind(envir.ena[,c("x","y")],clim.soil.d[,1]),res=res(ele.us),crs=crs(ele.us));plot(p)
envir.d=list(clim.d,soil.d,human.d,clim.human.d,clim.soil.d,clim.soil.human.d);names(envir.d)=c("clim","soil","human","clim.human","clim.soil","clim.soil.human")
#envir.d=list(clim.soil.pop.d,human.d,clim.human.d);names(envir.d)=c("clim.soil.pop","human","clim.human")
rm(envir.ena,clim,human,soil,clim.soil,clim.d,human.d,soil.d,clim.soil.d,clim.human.d,clim.soil.human.d)
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

RangeFragments=get(load("data analysis/AlienRangeFragments.alltaxa.ENA.new.RData"))
fragmentTable1=data.frame(fragBinomial=unique(names(RangeFragments)),fragSpecies=1:length(unique(names(RangeFragments))))
fragmentTable2=data.frame(fragNum=1:length(RangeFragments),fragSize=unlist(lapply(RangeFragments,length)),fragBinomial=names(RangeFragments))
fragmentTable=fragmentTable2%>%left_join(fragmentTable1,by="fragBinomial")
rm(fragmentTable1,fragmentTable2,RangeFragments.v2)

mf=function(i,RangeFragments) {return(data.frame(Tip_Label=names(RangeFragments)[i],domain.ID=RangeFragments[[i]],fragNum=i))}
IntroductionsFiltered=do.call(rbind,lapply(1:length(RangeFragments),mf,RangeFragments))

get.sim=function(REP,fragmentTable,RangeFragments,EnvironmentData,IntroductionsFiltered,adjCells.ena,CommSimilarity,mu,mat.name,inputProbSurface=TRUE){
	source("data analysis/Functions.R")
	#source("Code/2.Functions/Functions.R")
	#this function comes from Lovell et al. (2021).,Nature Ecology & Evolution, 5(3), 322-329 https://doi.org/10.6084/m9.figshare.13270406
	simResults<-SimulateSpread(Range="Alien",fragmentTable=fragmentTable,RangeFragments=RangeFragments,
									 EnvironmentData=EnvironmentData,IntroductionsFiltered=IntroductionsFiltered,adjCells=adjCells.ena,
									 SpeciesOptimum=NULL,CommSimilarity=CommSimilarity,plotIt=FALSE,envVar=c(),mu=mu,
									 seed="Random",niche="ClimateMatching",inputProbSurface=inputProbSurface)
	names(simResults)=names(RangeFragments)		
    outfile<-paste0("/blue/matthewthomas1/yunpeng.liu/data analysis/simulation.ena2/",mat.name,"/ER_",mat.name,"_",mu,"_",REP,".rda")
	save(simResults,file=outfile)
}

library(parallel)
no_cores <- 50#detectCores() - 1
### simulate environmental resistance model with different mu values
muval<-sort(c(0.5,1.5,3,5,10,30,100),decreasing =TRUE)
#muval=c(10,5,0.5,30)
mat.list=c(names(envir.d),"sp","trait","phylo")	
#mat.list=c("clim.human","clim.soil.pop","null")

gc() 
for (mat in 1:length(mat.list)){
	mat.name=mat.list[mat]
	print(mat.name)
	if (mat.name=="null"){
		mycl <- makePSOCKcluster(no_cores);
		parSapply(cl=mycl,X=1:100,get.sim,fragmentTable,RangeFragments,EnvironmentData=domain.ena,IntroductionsFiltered,adjCells.ena,
			CommSimilarity=NULL,mu=NA,mat.name="null",inputProbSurface=FALSE)
		stopCluster(mycl)	
	}else{
		if (mat.name%in%names(envir.d)){		
			simmat=envir.d[[mat.name]]
		}else{
			simmat=get(load(paste0("/blue/matthewthomas1/yunpeng.liu/data analysis/CommSimilarity/ENA.",mat.name,".rda")))
		}
		if (grepl("sp",mat.name)+grepl("sp2",mat.name)==0) simmat=1-simmat
		if (max(simmat,na.rm=T)>1) simmat=scales::rescale(simmat)
		#nstarts=1;toload=1:100
		nstarts=ifelse(mat.name=="trait",3,1)
	for(i in nstarts:length(muval)){
		if (nstarts>1){
			toload=list.files(paste0("/blue/matthewthomas1/yunpeng.liu/data analysis/simulation.ena2/",mat.name))%>%
					.[grep(paste0("ER_",mat.name,"_",muval[i],"_"),.)]%>%gsub(paste0("ER_",mat.name,"_",muval[i],"_"),"",.)%>%gsub(".rda","",.)%>%as.numeric();
			toload=c(1:100)[-toload]
		}else{toload=c(1:100)}
		
		mycl <- makePSOCKcluster(no_cores);		
		parSapply(cl=mycl,X=toload,get.sim,fragmentTable,RangeFragments,EnvironmentData=domain.ena,IntroductionsFiltered,adjCells.ena,
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
mat.list=c("null","clim","soil","human","clim.human","clim.soil","clim.soil.human","sp","phylo","trait")	
RangeFragments=get(load("data analysis/AlienRangeFragments.alltaxa.ENA.new.RData"))
fragmentTableA=data.frame(fragBinomial=unique(names(RangeFragments)),fragSpecies=1:length(unique(names(RangeFragments))))
fragmentTableB=data.frame(fragNum=1:length(RangeFragments),fragSize=unlist(lapply(RangeFragments,length)),fragBinomial=names(RangeFragments))
AlienfragmentTable=fragmentTableB%>%left_join(fragmentTableA,by="fragBinomial")
get.prob=function(i,toload,mat){
		 get.simdis=function(j,i){
			re=get(load(paste0("/blue/matthewthomas1/yunpeng.liu/data analysis/simulation.ena2/",mat,"/",j)))[[i]]  
			#re=get(load(paste0("data analysis/simulation.ena2/",mat,"/",j)))[[i]]  
		   return(re)
		  }
		  simOutput.t<-do.call(c,lapply(toload,get.simdis,i))
		  simOutput.freq=table(simOutput.t)/length(toload)
		  return(simOutput.freq)
		  rm(re,simOutput.t);gc()
	  }
intersectRanges<-function(x,y){return(length(intersect(x,y)))}	  
library(parallel)
no_cores <- 50#detectCores() - 1
muval<-sort(c(0.5,1.5,3,5,10,30,100),decreasing =TRUE)
nFrag=nrow(AlienfragmentTable);nReps=100;
overlapArrayMuVals=list();
for (mat in mat.list){
	if (mat=="null"){
		overlapArrayMuVals[[mat]]=array(dim=c(nFrag,1,nReps))
		toload<-list.files(paste0("/blue/matthewthomas1/yunpeng.liu/data analysis/simulation.ena2/",mat,"/"))
		toload=gtools::mixedsort(sort(toload))
		#tmp=get(load(paste0("data analysis/simulation.ena/",mat,"/",toload[1])))
		#AlienfragmentTable=AlienfragmentTable.all%>%filter(fragBinomial%in%names(tmp))
		print(c(mat,length(toload)))
		
		mycl <- makePSOCKcluster(no_cores);
		InvasionProb<-parLapply(cl=mycl,X=1:nFrag,get.prob,toload,mat)  
		names(InvasionProb)<-AlienfragmentTable$fragNum	
		outfile<-paste("/blue/matthewthomas1/yunpeng.liu/data analysis/InvasionProbabilites.ena2/Alien_Random",mat,"rda",sep=".")
		save(InvasionProb,file=outfile)	 
		stopCluster(mycl);
		rm(InvasionProb);gc()
		
		#Fragments=subset(RangeFragments,names(RangeFragments)%in%names(tmp))	
		for(i in 1:length(toload)){
			load(paste0("/blue/matthewthomas1/yunpeng.liu/data analysis/simulation.ena2/",mat,"/",toload[i]))
			overlap<-mapply(intersectRanges,RangeFragments,simResults)/AlienfragmentTable$fragSize
			overlapArrayMuVals[[mat]][,1,i]<-overlap
		}				
	}else{
		overlapArrayMuVals[[mat]]=array(dim=c(nFrag,length(muval),nReps))		
		for(x in 1:length(muval)){
			SimFiles<-list.files(paste0("/blue/matthewthomas1/yunpeng.liu/data analysis/simulation.ena2/",mat,"/"))	
			toload<-SimFiles[grep(paste0("ER_",mat,"_",muval[x],"_"),SimFiles)]
			toload=gtools::mixedsort(sort(toload))
			print(c(mat,muval[x],length(toload)))
			
			mycl <- makePSOCKcluster(no_cores);		
			InvasionProb<-parLapply(cl=mycl,X=1:nFrag,get.prob,toload,mat)  
			names(InvasionProb)<-AlienfragmentTable$fragNum	
			outfile<-paste("/blue/matthewthomas1/yunpeng.liu/data analysis/InvasionProbabilites.ena2/Alien_Random",mat,muval[x],"rda",sep=".")
			#outfile<-paste("data analysis/InvasionProbabilites.ena/Alien_Random",mat,muval[x],"rda",sep=".")
			save(InvasionProb,file=outfile)	  	
			stopCluster(mycl);
			rm(InvasionProb);gc()
			
			for(i in 1:length(toload)){
			  load(paste0("/blue/matthewthomas1/yunpeng.liu/data analysis/simulation.ena2/",mat,"/",toload[i]))
			  #load(paste0("data analysis/simulation.ena2/",mat,"/",toload[i]))
			  overlap<-mapply(intersectRanges,RangeFragments,simResults)/AlienfragmentTable$fragSize
			  overlapArrayMuVals[[mat]][,x,i]<-overlap		 
			}	
		  }
		}
	  save(overlapArrayMuVals,file=paste0("/blue/matthewthomas1/yunpeng.liu/data analysis/InvasionProbabilites.ena2/ModelOverlapScores.rda"))
}
#save(overlapArrayMuVals,file=paste0("data analysis/InvasionProbabilites.ena/ModelOverlapScores.human.rda"))	  
########################################################################################
# calculate species level accuracy (i.e. overlap) scores
########################################################################################
library(dplyr);library(matrixStats)
load("/blue/matthewthomas1/yunpeng.liu/data analysis/InvasionProbabilites.ena2/ModelOverlapScores.rda")
RangeFragments.ena.inv=get(load("data analysis/AlienRangeFragments.alltaxa.ENA.new.RData"))

# calculate species accuracy scores for each model
source("data analysis/Functions.R")
muval<-sort(c(0.5,1.5,3,5,10,30,100),decreasing =TRUE)
mat.list=names(overlapArrayMuVals)	
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
		scoreMatrix<-calculateSpeciesOverlapScores(FragTable=AlienfragmentTable,overlapArray=overlapArrayMuVals[[mat]])##this function comes from Lovell et al. (2021).,Nature Ecology & Evolution, 5(3), 322-329
		OScoresMuvalsAlienMean[[mat]][,1]<-rowMeans(scoreMatrix)
		OScoresMuvalsAlienMedian[[mat]][,1]<-rowQuantiles(scoreMatrix,probs=c(0.5))
		OScoresMuvalsAlienSE[[mat]][,1]<-apply(scoreMatrix,1,se)	  
	}else {		
		for(x in 1:length(muval)){
			AlienfragmentTable$FocalModel<-x
			scoreMatrix<-calculateSpeciesOverlapScores(FragTable=AlienfragmentTable,overlapArray=overlapArrayMuVals[[mat]]) #this function comes from Lovell et al. (2021)., https://doi.org/10.6084/m9.figshare.13270406
			OScoresMuvalsAlienMean[[mat]][,x]<-rowMeans(scoreMatrix)
			OScoresMuvalsAlienMedian[[mat]][,x]<-rowQuantiles(scoreMatrix,probs=c(0.5))
			OScoresMuvalsAlienSE[[mat]][,x]<-apply(scoreMatrix,1,se)
			  #print(x)
		}
	}	
}
save(OScoresMuvalsAlienMean,OScoresMuvalsAlienMedian,OScoresMuvalsAlienSE,
	file=paste0("/blue/matthewthomas1/yunpeng.liu/data analysis/InvasionProbabilites.ena2/SpeciesModelOverlapScores.rda"))

########################################################################################
## model comparison
########################################################################################
library(ggpubr);library(dplyr);library(tidyr)
load("data analysis/InvasionProbabilites.ena2/SpeciesModelOverlapScores.rda")
load("data analysis/AlienRangeFragments.alltaxa.ENA.new.RData")
#mat.list=c("null","clim","clim.pca","phylo","phylo2","sp","sp2","trait","trait.turnover")	
mat.list=c("null","clim","soil","human","clim.human","clim.soil","clim.soil.human","sp","phylo","trait")	
muval<-sort(c(0.5,1.5,3,5,10,30,100),decreasing =TRUE)
splist=read.csv("data analysis/US-RIIS.csv")%>%select(kingdom,phylum,class,scientificName,habitat)%>%filter(scientificName%in%names(RangeFragments.v2))%>%
	filter(!phylum%in%c("Rhodophyta","Chlorophyta","Charophyta"))
	#filter(!class%in%c("Teleostei","Ascidiacea","Malacostraca","Maxillopoda","Ostracoda","Branchiopoda"))
# aquatic=rbind(splist[grep("Aquatic",splist$habitat),],splist[grep("Marine",splist$habitat),])
# splist=splist%>%filter(!scientificName%in%aquatic$scientificName)
taxa=c("Arthropoda","Chordata","Other animal","Tracheophyta")#other animal incl. Mollusca,Annelida,Nematoda, etc.
splist%>%group_by(kingdom,phylum)%>%tally()

##step 1: statistic
re=c()
for (taxa.i in taxa){
	spt=splist[splist$phylum==taxa.i,"scientificName"]
	if (taxa.i%in%"Other animal") spt=splist[splist$kingdom=="Animalia"&(!splist$phylum%in%c("Chordata","Arthropoda")),"scientificName"]
	
	null.mean=OScoresMuvalsAlienMean[["null"]]%>%subset(rownames(.)%in%spt)
	null.mean=null.mean*100
	null.se=OScoresMuvalsAlienSE[["null"]]%>%subset(rownames(.)%in%spt)
	for (mat in mat.list[-1]){
		for (mu in 1:length(muval)){
		er.mean=as.matrix(OScoresMuvalsAlienMean[[mat]][,mu])%>%subset(rownames(.)%in%spt);colnames(er.mean)=muval[mu]
		er.mean=er.mean*100
		er.se=as.matrix(OScoresMuvalsAlienSE[[mat]][,mu])%>%subset(rownames(.)%in%spt);colnames(er.se)=muval[mu]
		se=function(x) sd(x, na.rm = TRUE)/length(na.omit(x))
		
		#how many species does ER outperform random dispersal?
		diffScore<-er.mean-null.mean
		a=length(which(diffScore>0))/length(diffScore)*100

		# how many species can we reject the null model?
		diffScore<-(er.mean-2*er.se)-(null.mean+2*null.se)
		b=length(diffScore[diffScore>0])/length(diffScore)*100
		
		#paired t test		
		sta=t.test(er.mean, null.mean,paired=T)
		tmp=data.frame(taxa=taxa.i,model=mat,mu=muval[mu],a,b,mean.s=mean(er.mean),se=se(er.mean),t=sta$statistic,df=sta$parameter,p=sta$p.value,null.mean=mean(null.mean),null.se=se(null.mean))
		re=rbind(re,tmp)
		}
	}
}
re%>%select(taxa,null.mean,null.se)%>%distinct

re2=re%>%group_by(taxa,model)%>%summarize(a=max(a),b=max(b),se.s=se[which(mean.s==max(mean.s))],mean.s=max(mean.s),
	t=t[which(p==min(p))],df=df[which(p==min(p))],p=min(p))%>%as.data.frame()%>%distinct
Px=ifelse(re2$p<0.001,"***",ifelse(re2$p<0.01,"**",ifelse(re2$p<0.05,"*","ns")))
re2=cbind(re2,Px)
#write.csv(re2,"data analysis/model.perform2.csv")
re.mu=re%>%group_by(model,taxa)%>%reframe(a=mu[which(a==max(a))],b=mu[which(b==max(b))],mean.s=mu[which(mean.s==max(mean.s))],diff.with.null=mu[which(diff.with.null==max(diff.with.null))])%>%as.data.frame()
#write.csv(re.mu,"data analysis/re.mu.csv")

## step2: compare different ER models
range.size=data.frame(fragSize=unlist(lapply(RangeFragments.v2,length)),Species=names(RangeFragments.v2))%>%group_by(Species)%>%summarize(rangesize=sum(fragSize))
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

score2=score
score2$score.mean=round(score2$score.mean,2)
score2=score2%>%group_by(Species)%>%reframe(bestmodel=model[which(score.mean==max(score.mean))],bestmu=mu[which(score.mean==max(score.mean))])
save(score2,file="data analysis/score2.rda")#the best model for each species

load("data analysis/score2.rda")
dat=c()
for (taxa.i in taxa){
	sp=splist[splist$phylum==taxa.i,"scientificName"]
	if (taxa.i%in%"Other animal") sp=splist[splist$kingdom=="Animalia"&(!splist$phylum%in%c("Chordata","Arthropoda")),"scientificName"]
	score3=score2%>%filter(Species%in%sp)
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

re2$model=factor(re2$model,levels=mat.list[-1])	
ggplot(re2, aes(x=model, y=a*100)) + 
		geom_bar(stat = "identity")+
		#annotate("text",x=7,y=0.95,size=8,label=taxa.i,color="red")+
		labs(x="Type of models",y="proportion of species outperform null model (%)")+   
		facet_wrap(~ taxa,scales="free_y")+theme
	
#use the same mu for each model	
# dat=c()
# for (taxa.i in taxa){
	# sp=splist[splist$phylum==taxa.i,"scientificName"]
	# if (taxa.i%in%"Other animal") sp=splist[splist$kingdom=="Animalia"&(!splist$phylum%in%c("Chordata","Arthropoda")),"scientificName"]
	# score2=c()
	# for (mat in mat.list){
		# if (mat=="null") tmp=score%>%filter(Species%in%sp&is.na(mu)) else{
			# bestmu=unique(re.mu[re.mu$model==mat&re.mu$taxa==taxa.i,"a"])
			# tmp=score%>%filter(Species%in%sp&model==mat&mu==bestmu)}
		# score2=rbind(score2,tmp)
	# }
	# score3=score2
	# score3$score.mean=round(score3$score.mean,2)
	# score3=score3%>%group_by(Species)%>%reframe(bestmodel=model[which(score.mean==max(score.mean))])
	# nModelsTied<-score3%>%group_by(Species)%>%tally()
	# Weights<-1/nModelsTied$n
	# WeightsVec<-rep(Weights,nModelsTied$n)
	# BestModelsTable<-cbind(score3,WeightsVec)%>%group_by(bestmodel)%>%summarize(spnum=sum(WeightsVec))
	# BestModelsTable$bestmodel=factor(BestModelsTable$bestmodel,levels=mat.list)
	# dat=rbind(dat,data.frame(taxa=taxa.i,BestModelsTable))
# }

	
#compare different models
# score2=score2%>%left_join(splist,by=c("Species"="scientificName"))
# ggplot(score2,aes(x = rangesize,y=score.mean,fill=model))+
		# geom_ribbon(aes(x = rangesize,ymin=score.mean-score.se, ymax=score.mean+score.se,fill=model),alpha=0.3,show.legend=FALSE)+
		# geom_line(aes(x=rangesize, y=score.mean,colour=model),size=0.8,linetype=1,alpha=0.8)+
		# # scale_color_discrete(labels = c("Random dispersal", "Species similarity","Phylogenetic similarity","Trait distance","Trait turnover",
			# # "Climatic similarity","Soil similarity","Human activity","Cliamte and soil", "Climate, soil and human activity"))+
		# # geom_errorbar(aes(ymin=score.mean-score.se, ymax=score.mean+score.se,color=model),linewidth=1,width=150,show.legend=FALSE,alpha=0.6)+
		# # geom_point(size=3,shape=21,color="black",alpha=0.6) + 
		# #geom_smooth(aes(x = rangesize,y=score.mean,color=model),data=score2,method = "loess",size=1.5,show.legend=FALSE,se =TRUE,span=0.8,linetype=1,alpha=0.5)+
		# annotate("text",x=2000,y=0.95,size=8,label=taxa.i,color="red")+
		# facet_wrap(~ phylum,scales="free")+
		# labs(x="Range sizes of species",y="Model accuracy scores",colour="Model") +theme+theme(legend.position =c(0.4,0.8))		

############################################################
######### map predicted ranges for different ER model #########
############################################################
library(dplyr);
load("data analysis/domain.ena.rda")
mat.list=c("null","clim","soil","human","clim.human","clim.soil","clim.soil.human","sp","phylo","trait")	
re.mu=read.csv("data analysis/re.mu.csv")[,-1]
load("data analysis/AlienRangeFragments.alltaxa.ENA.new.RData")
AlienfragmentTable=data.frame(fragNum=1:length(RangeFragments.v2),fragSize=unlist(lapply(RangeFragments.v2,length)),fragBinomial=names(RangeFragments.v2))
splist=read.csv("data analysis/US-RIIS.csv")%>%select(kingdom,phylum,class,scientificName,habitat)%>%filter(scientificName%in%names(RangeFragments.v2))%>%
	filter(!phylum%in%c("Rhodophyta","Chlorophyta","Charophyta"))
taxa=c("Arthropoda","Chordata","Other animal","Tracheophyta")#other animal incl. Mollusca,Annelida,Nematoda, etc.
#richness of invasive species based on model simulation
get.rich=function(taxa.i,splist,AlienfragmentTable,re.mu,mat=NULL){
	require(matrixStats);require(dplyr)
	sp=splist[splist$phylum==taxa.i,"scientificName"]
	if (taxa.i%in%"Other animal") sp=splist[splist$kingdom=="Animalia"&(!splist$phylum%in%c("Chordata","Arthropoda")),"scientificName"]	
	
	if (is.null(mat)){
		#using the best model of each species
		score=get(load("data analysis/score2.rda"))%>%filter(Species%in%sp&(!bestmodel%in%"null"))%>%
			 .[order(.$bestmu,decreasing=T),]%>%select(-bestmu)%>%#.[order(match(.$bestmodel, mat.list)),]%>%
			  .[!duplicated(.[,'Species']),]		
		models=unique(score$bestmodel)
		ER.model=c()
		for (mat in models){
			sp2=score%>%filter(bestmodel==mat)%>%as.data.frame()%>%.[,"Species"]
			bestmu=re.mu[re.mu$model==mat&re.mu$taxa==taxa.i,"a"][1]
			ER.model.ori=get(load(paste("data analysis/InvasionProbabilites.ena2/Alien_Random",mat,bestmu,"rda",sep=".")))
			names(ER.model.ori)=AlienfragmentTable$fragBinomial
			ER.model.t=subset(ER.model.ori,names(ER.model.ori)%in%sp2)
			ER.model=c(ER.model,ER.model.t)
		}
	}else{	
		#using the same ER model for all species
		bestmu=unique(re.mu[re.mu$model==mat&re.mu$taxa==taxa.i,"a"])[1]
		if (mat%in%"null") ER.model=get(load(paste("data analysis/InvasionProbabilites.ena2/Alien_Random",mat,"rda",sep="."))) else {
				ER.model=get(load(paste("data analysis/InvasionProbabilites.ena2/Alien_Random",mat,bestmu,"rda",sep=".")))}
		names(ER.model)=AlienfragmentTable$fragBinomial
		ER.model=subset(ER.model,names(ER.model)%in%sp)
	}
	mf=function(x)data.frame(domain.ID=names(x),Freq=as.numeric(x))	
	rich.ER=do.call(rbind,lapply(ER.model,mf))
	rich.ER[,1]=as.numeric(rich.ER[,1])	
	rich.ER=rich.ER%>%group_by(domain.ID)%>%summarize(ER=sum(Freq))
	colnames(rich.ER)[2]=taxa.i	
	return(rich.ER)
}

#using sp-based ER model
rich.ER=list()
for (taxa.i in taxa) rich.ER[[taxa.i]]=get.rich(taxa.i,splist,AlienfragmentTable,re.mu,mat="soil")

#using the best model of each species
rich.ER=list()
for (taxa.i in taxa) {rich.ER[[taxa.i]]=get.rich(taxa.i,splist,AlienfragmentTable,re.mu);print(taxa.i)}

#actural richness of invasive species
rich.current=list()
for (taxa.i in taxa) {
sp=splist[splist$phylum==taxa.i,"scientificName"]
if (taxa.i%in%"Other animal") sp=splist[splist$kingdom=="Animalia"&(!splist$phylum%in%c("Chordata","Arthropoda")),"scientificName"]	
RangeFragments=subset(RangeFragments.v2,names(RangeFragments.v2)%in%sp)
AlienRichness<-table(unlist(RangeFragments))
rich.current[[taxa.i]]=data.frame(domain.ID=domain.ena$domain.ID,richness=as.numeric(AlienRichness)[match(domain.ena$domain.ID,as.numeric(names(AlienRichness)))])
colnames(rich.current[[taxa.i]])[2]=paste0("current_",taxa.i)
}

domain=domain.ena[,c("domain.ID","x","y")]
domain=c(list(domain),rich.ER,rich.current) %>% purrr::reduce(left_join,by="domain.ID")
#ma=ceiling(max(apply(domain[,-c(1:3)],2,max,na.rm=T)))

#plot richness
library(ggpubr)
get.p=function(taxa.i,domain,type,cor.only=FALSE){
the=theme_bw()+theme(legend.background=element_rect(fill='transparent'),
		legend.text=element_text(size=12),
		legend.title=element_text(face="bold",size=10))
		#legend.position=c(0.75,0.2))
datp=domain%>%select(x,y,contains(taxa.i))
colnames(datp)[c(3:4)]=c("pred","currt")
if (cor.only){
re=data.frame(taxa=taxa.i,r=round(cor(na.omit(datp[,c(3:4)])),2)[1,2])
return(re)
}else{
ma=ceiling(max(apply(datp[,c(3:4)],2,max,na.rm=T)))

if (type=="observed"){	
p=ggplot() +  geom_tile(data = datp, aes(x = x, y = y, fill = currt)) +
  coord_quickmap() +  labs(x="",y="",fill="Richness")+
  scale_fill_gradientn(colours = viridis::viridis(99),na.value = "lightgray",limits=c(0,ma))+the+
  annotate("text",x=-73,y=30,size=3.5,label=taxa.i,color="black") 
}
if (type=="predict"){
p=ggplot() +  geom_tile(data = datp, aes(x = x, y = y, fill = pred)) +
  coord_quickmap() +  labs(x="",y="",fill="Richness")+
  scale_fill_gradientn(colours = viridis::viridis(99),na.value = "lightgray",limits=c(0,ma))+the+
  annotate("text",x=-73,y=30,size=3.5,label=taxa.i,color="black") 
p=p+annotate("text",x=-72,y=25,size=5,label=paste0("r = ",round(cor(na.omit(datp[,c(3:4)])),2)[1,2])) 
  
}  
return(p)
}
}
#plt=lapply(taxa,get.p,domain,"observed")
plt=lapply(taxa,get.p,domain,"predict")
ggarrange(plt[[1]],plt[[2]],plt[[3]],plt[[4]],labels="auto",font.label = list(size = 20))

#obtain cor for each model
re=c()
for (mol in mat.list[-c(1,2,3,8)]){
rich.ER=list()
for (taxa.i in taxa) rich.ER[[taxa.i]]=get.rich(taxa.i,splist,AlienfragmentTable,re.mu,mat=mol)
domain=domain.ena[,c("domain.ID","x","y")]
domain=c(list(domain),rich.ER,rich.current) %>% purrr::reduce(left_join,by="domain.ID")
tmp=do.call(rbind,lapply(taxa,get.p,domain,"predict",cor.only=TRUE))
re=rbind(re,data.frame(mat=mol,tmp))
}
r.all=read.csv("data analysis/r.of.invasion.richness.csv")
r.all$Model=factor(r.all$Model,levels=mat.list)
ggplot(r.all,aes(x=Model, y=Pearson.r)) +
    geom_bar(position="dodge", stat="identity",show.legend=F,fill="lightgray",color="black")+
	#coord_cartesian(ylim=c(0.35,0.9))+
   # geom_errorbar(aes(x=type1, ymin=value*100-se*196, ymax=value*100+se*196,colour=type2), width=0.4,  alpha=0.9, size=1.3,position=position_dodge(.9),show.legend=F)+
	labs(x="Type of models",y="Correlation betweeen predicted and observed invasion richness")+
	#scale_x_discrete(labels=c("<1%" = "<1%","<5%"="1-5%","<10%"="5-10%",">=10%"=">=10%"))+	
	facet_wrap(~Taxon, scales="free_y",ncol=1)+
	theme_classic()+theme(axis.text = element_text(size=12,color='black'),
		#axis.text.x =element_text(angle=15),
		axis.title = element_text(size=12,color='black'),
		strip.text = element_text(colour = 'black', face="italic",size = rel(1.2)), 
		strip.background = element_rect(fill = 'white', colour = 'darkgray', size = rel(2), linetype = 2)
		)
#####################################################################
### calculate expected alien range size for different ER thresholds
#####################################################################
library(dplyr);
inputER<-c(0.1,0.3,0.4,seq(0.45,0.8,0.01),0.85,0.9,0.98)
load("data analysis/AlienRangeFragments.alltaxa.ENA.new.RData")
load("data analysis/connectedGraph.ena.rda")

get.simrange=function(x,RangeFragments,simmat,g1cc,mat.name){
	require(dplyr)
	fragmentTable1=data.frame(fragBinomial=unique(names(RangeFragments)),fragSpecies=1:length(unique(names(RangeFragments))))
	fragmentTable2=data.frame(fragNum=1:length(RangeFragments),fragSize=unlist(lapply(RangeFragments,length)),fragBinomial=names(RangeFragments))
	FragmentTable=fragmentTable2%>%left_join(fragmentTable1,by="fragBinomial")
	#for(x in 1:inputER){  
	source("data analysis/Functions.R")
	output<-runConnectedComponentsModel(threshold=x,nRandomStarts=500,fragmentTable=FragmentTable,RangeFragments=RangeFragments,CommSimilarity=simmat,g1cc=g1cc) #this function comes from Lovell et al. (2021).,Nature Ecology & Evolution, 5(3), 322-329
	#save the simulated ranges
	fragmentTable<-output[[1]]
	simRanges<-output[[2]]
	save(simRanges,fragmentTable,file=paste0("/blue/matthewthomas1/yunpeng.liu/data analysis/simulation.ena.predict2/PredictedRangeSize_",x,".",mat.name,".rda")) 
	rm(fragmentTable1,fragmentTable2,output,RangeFragments,fragmentTable,simRanges);gc()
	}
	
library(parallel)
no_cores <- 45#detectCores() - 1	
#mat.list=c("sp","clim","soil","phylo2","human","trait","clim.soil","clim.soil.pop","trait.turnover")	
simmat=get(load(paste0("/blue/matthewthomas1/yunpeng.liu/data analysis/CommSimilarity/ENA.sp.rda")))
mycl <- makePSOCKcluster(no_cores);
parSapply(cl=mycl,X=inputER,get.simrange,RangeFragments.v2,simmat,g1cc,mat.name="sp")
stopCluster(mycl)
gc()	
#}

# fit quantile regressions between observed and predicted range sizes for different ER thresholds and, for quantile egressions, tau values 
splist=read.csv("data analysis/US-RIIS.csv")%>%select(kingdom,phylum,class,scientificName)%>%filter(scientificName%in%names(RangeFragments.v2))
taxa=c("Arthropoda","Chordata","Other animal","Tracheophyta")	
# different ER thresholds and tau values to explore
get.BestERvalues=function(x,mat,taxa,splist){
	require(Metrics);require(quantreg);require(dplyr)
	load(paste0("/blue/matthewthomas1/yunpeng.liu/data analysis/simulation.ena.predict2/PredictedRangeSize_",x,".",mat,".rda"))	  
	fragmentTable$fragSizeLOG<-log10(fragmentTable$fragSize)
	fragmentTable$predfragSizeLOG<-log10(fragmentTable$PredictedFragSize)
	input=c()
	for (taxa.i in taxa){
		sp=splist[splist$phylum==taxa.i,"scientificName"]
		if (taxa.i%in%"Other animal") sp=splist[splist$kingdom=="Animalia"&(!splist$phylum%in%c("Chordata","Arthropoda")),"scientificName"]	
		fragmentTable.t=fragmentTable%>%filter(fragBinomial%in%sp)
		underpredicted=length(which(fragmentTable.t$fragSize>fragmentTable.t$PredictedFragSize))/nrow(fragmentTable.t)*100
		
		fragmentTable.t.sp=fragmentTable.t%>%group_by(fragBinomial)%>%summarize(fragSizeLOG=log10(sum(fragSize)),predfragSizeLOG=log10(sum(PredictedFragSize)))%>%as.data.frame
		underpredicted.sp=length(which(fragmentTable.t.sp$fragSizeLOG>fragmentTable.t.sp$predfragSizeLOG))/nrow(fragmentTable.t.sp)*100		  
		mf=function(tau.t,fragtab){
			if (tau.t==1) rqfit<-lm(fragSizeLOG ~ predfragSizeLOG, data = fragmentTable) else {
			rqfit <- rq(fragSizeLOG ~ predfragSizeLOG, data = fragtab,tau=tau.t,method="fn")}
			predicted<-predict(rqfit,newdata =data.frame(predfragSizeLOG=fragtab$predfragSizeLOG))
			return(rmse(fragtab$predfragSizeLOG, predicted)) 
		 }		  
		 tmp=data.frame(taxa=taxa.i,model=mat,inputER=x,tau=c(0.5,0.9,0.95,0.99,1),
			rmseQ=do.call(rbind,lapply(c(0.5,0.9,0.95,0.99,1),mf,fragmentTable.t)),
			underpredicted=underpredicted,
			rmseQspecies=do.call(rbind,lapply(c(0.5,0.9,0.95,0.99,1),mf,fragmentTable.t.sp)),
			underpredicted.sp=underpredicted.sp)
		input=rbind(input,tmp)
	 }	    
	return(input)	
}

mycl <- makePSOCKcluster(no_cores);
input=do.call(rbind,parLapply(cl=mycl,X=inputER,get.BestERvalues,"sp",taxa,splist))
stopCluster(mycl)
save(input,file="/blue/matthewthomas1/yunpeng.liu/data analysis/simulation.ena.predict2/BestERvalues.sp.rda")

load("data analysis/simulation.ena.predict2/BestERvalues.sp.rda")
#get optimal tau and ER thres values
# #based on fragment
# stat.ERQ=input%>%group_by(taxa,tau)%>%reframe(inputER=inputER[which(rmseQ==min(rmseQ))],underpred=underpredicted[which(rmseQ==min(rmseQ))])%>%
	# group_by(taxa)%>%reframe(tau.opt=tau[which(underpred==min(underpred))],bestERQ=inputER[which(underpred==min(underpred))],underpred.por=min(underpred))%>%as.data.frame
 
#based on species
stat.ERQ=input%>%group_by(taxa,tau)%>%reframe(inputER=inputER[which(rmseQspecies==min(rmseQspecies))],underpred=underpredicted[which(rmseQspecies==min(rmseQspecies))])%>%
	group_by(taxa)%>%reframe(tau.opt=tau[which(underpred==min(underpred))],bestERQ=inputER[which(underpred==min(underpred))],underpred.por=min(underpred))%>%as.data.frame

## plot root mean square error of different ERs
library(ggpubr);library(dplyr)
rmse=c();
for (taxa.i in taxa) rmse=rbind(rmse,input%>%filter(taxa==taxa.i,tau==stat.ERQ[stat.ERQ$taxa==taxa.i,"tau.opt"]))
theme=theme_bw()+theme(axis.text = element_text(size=15,color='black'),
			axis.title = element_text(size=18,color='black'),
			legend.background=element_rect(fill='transparent'),
			legend.text=element_text(face="bold",size=15),
			legend.title=element_text(face="bold",size=18))
rmse$taxa=factor(rmse$taxa,levels=taxa)
p=ggplot(rmse,aes(x = inputER,y=rmseQspecies,colour=taxa))+
	geom_line(linewidth=1.2,linetype=1,alpha=0.8)+
	#scale_color_discrete(labels = c("Null model", "Climatic distance","Species composition","Phylogenetic composition","Trait distance","Trait turnover"))+
	labs(x="Threshold of Environmental resistance",y="Root mean square error",colour="Taxa") +theme+theme(legend.position =c(0.2,0.75))
for (i in 1:length(taxa)) {
	bestERQ=stat.ERQ[stat.ERQ$taxa==taxa[i],"bestERQ"]
	p=p+geom_vline(xintercept=bestERQ,col=scales::hue_pal()(4)[i],linewidth=1.2,linetype="longdash",alpha=0.7)+
		annotate("text", x=bestERQ-0.03 , y=2.5+i*0.1 ,size=5,label=bestERQ,color=scales::hue_pal()(4)[i])
}

## generate matrix of invasion probabilities using the best ER model #### 
load("data analysis/Domain.ena.rda")
get.InvasionProb=function(mat.name,taxa,splist,stat.ERQ,domain.ena,RangeFragments.ena.inv,current=TRUE){
	nRandomStarts<-500
	InvasionProbabilityMatrix=list()
	for (taxa.i in taxa){
		sp=splist[splist$phylum==taxa.i,"scientificName"]
		if (taxa.i%in%"Other animal") sp=splist[splist$kingdom=="Animalia"&(!splist$phylum%in%c("Chordata","Arthropoda")),"scientificName"]	
		nSpecies<-length(sp)
		bestERQ=stat.ERQ[stat.ERQ$taxa==taxa.i,"bestERQ"]
		if(current) load(paste0("data analysis/simulation.ena.predict2/PredictedRangeSize_",bestERQ,".",mat.name,".rda"))	
		if(!current) load(paste0("data analysis/SDM.simulation.ena.predict/PredictedRangeSize_",bestERQ,"_",mat.name,".rda"))
		InvasionProbabilityMatrix[[taxa.i]]<-matrix(0,ncol=nrow(domain.ena),nrow=nSpecies);
		colnames(InvasionProbabilityMatrix[[taxa.i]])=domain.ena$domain.ID;
		rownames(InvasionProbabilityMatrix[[taxa.i]])=sp	
		for(i in 1:nSpecies){ # loop over each species  
		  focRowAlien<-which(fragmentTable$fragBinomial==sp[i])
		  nFrags<-length(focRowAlien)
		   # for each alien range fragment calculate the probability of invading each cell
		  InvasionProbs<-matrix(0,ncol=nrow(domain.ena),nrow=nFrags)
		  for(x in 1:nFrags){
			fragSize<-fragmentTable$fragSize[focRowAlien[x]]
			InvasionProbs[x,as.numeric(names(simRanges[[focRowAlien[x]]]))]<-as.numeric(simRanges[[focRowAlien[x]]])/min(nRandomStarts,fragSize)    
		  }	  
		  # the overall probability of an alien species invading a cell (a cell can be invaded from multiple fragments)
		  InvasionProbVec<-1-matrixStats::colProds(1-InvasionProbs)
		  
		  # set invasion probabilities of cells that are already invaded to 1
		  cells=subset(RangeFragments.ena.inv,names(RangeFragments.ena.inv)==sp[i])%>%unlist%>%unique%>%as.numeric%>%sort
		  InvasionProbVec[cells]<-1	  
		  InvasionProbabilityMatrix[[taxa.i]][i,]<-InvasionProbVec	  
		}	
	}
	return(InvasionProbabilityMatrix)
}		

InvasionProbabilityMatrix<-get.InvasionProb("sp",taxa,splist,stat.ERQ,domain.ena,RangeFragments.ena.inv=RangeFragments.v2)
InvasionProb<-do.call(cbind,lapply(InvasionProbabilityMatrix,colSums))
	
#devided by range size
#rangesize=AlienfragmentTable%>%group_by(fragBinomial)%>%summarize(ranges=sum(fragSize))%>%as.data.frame()

#plot richness
library(ggpubr)
the=theme_bw()+theme(legend.background=element_rect(fill='transparent'),
		legend.text=element_text(size=10),
		legend.title=element_text(face="bold",size=10))		
plt=list()
for (taxa.i in taxa){
domain=cbind(domain.ena[,c("x","y")],InvasionProb[,taxa.i])
colnames(domain)[3]="InvasionProb"
plt[[taxa.i]]=ggplot() +  geom_tile(data = domain, aes(x = x, y = y, fill = InvasionProb)) +
  coord_quickmap() +  labs(x="",y="",fill="Potential \ninvasion\nrichness")+
  annotate("text",x=-72.5,y=30,size=4.3,label=taxa.i,color="black")+
  scale_fill_gradientn(colours = viridis::viridis(99),na.value = "lightgray")+the#+theme(legend.position=c(0.8,0.25))
}
ggarrange(plt[[1]],plt[[2]],plt[[3]],plt[[4]],labels="auto",font.label = list(size = 20))	
############################################################################################################
####get future climate data from worldclim and project future distributions and conduct SDM ################################
############################################################################################################
#impute missing climate value
library(dplyr);
load("data analysis/Domain.ena.rda")
toload=list.files("data analysis/future climate/")	
cal.climdis=function(name,domain.ena,RangeFragments.ena.inv,glcc,no_cores){
	require(dplyr)
	#caculate future climate similarity
	load(paste0("data analysis/future climate/",name))
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
	save()
	}
	
library(parallel)
no_cores <- detectCores() - 1	
for (name in toload) {cal.climdis(name,domain.ena,RangeFragments.ena.inv,glcc,no_cores);print(name)}

library(sf);library(dplyr);library(raster)
clim.list=list.files("data analysis/future climate/")
load("data analysis/Domain.ena.rda")
ENA.name=c("ME","NH","VT","NY","MA","RI","CT","NJ","PA","DE","MD","MI","OH","IN","IL","WI","WV","VA","NC","TN","KY","SC","GA","AL","MS","FL")
ENA <- read_sf('data analysis/cb_2018_us_state_20m/cb_2018_us_state_20m.shp')%>%filter(STUSPS%in%ENA.name)%>%as(., "Spatial")
# Use Variable inflation factor (VIF) to select the best variables for your model.
library(usdm)	
for (i in c(clim.list,"current")){	
	if (i=="current") {
		v=vifstep(domain.ena%>%dplyr::select(starts_with("bio")),th=10)
		var.list=attributes(v)$results[,1]
		ENA.rst2=rasterFromXYZ(domain.ena[,c("x","y",var.list)],res=res(ENA.rst),crs=crs(ENA.rst))	
		outfile="current"
		bioclim=domain.ena%>%dplyr::select(c("x","y",starts_with("bio")))
	} else{
		rst.10m <- brick(paste0("data analysis/future climate/",i))
		ENA.rst=ENA%>%raster::mask(rst.10m,.)%>%raster::crop(extent(ENA))#%>%aggregate(fact=1.5)
		Domain.ena=raster::as.data.frame(ENA.rst,xy=TRUE)#
		Domain.ena=na.omit(Domain.ena)	#remove ocean cells	
		bioclim=domain.ena%>%dplyr::select(-starts_with("bio"))%>%left_join(Domain.ena,by=c("x","y"))#%>%dplyr::select(domain.ID,ele,starts_with("bio"))
		miss=bioclim%>%filter(is.na(bio01))
		imp=function(i,bioclim){
			cells.row=(miss[i,]$row-5):(miss[i,]$row+5)
			cells.col=(miss[i,]$col-5):(miss[i,]$col+5)
			impute=bioclim%>%filter(row%in%cells.row&col%in%cells.col)%>%dplyr::select(starts_with("bio"))%>%colMeans(na.rm =TRUE)
			return(impute)
		}
		miss.imp=do.call(rbind,lapply(1:nrow(miss),imp,bioclim))
		miss.imp=cbind(miss%>%dplyr::select(-starts_with("bio")),miss.imp)
		bioclim=rbind(bioclim%>%filter(!is.na(bio01)),miss.imp)
		
		#select VIF<10 envars
		v=vifstep(bioclim[,c(12:30)],th=10)
		var.list=attributes(v)$results[,1]
		ENA.rst2=rasterFromXYZ(bioclim[,c("x","y",var.list)],res=res(ENA.rst),crs=crs(ENA.rst))	
		outfile=gsub("wc2.1_10m_bioc_","",i);outfile=gsub(".tif","",outfile)	
	}		
	re=list(ENA.rst2,bioclim)
	save(re,file=paste0("data analysis/future climate/",outfile,".rda"))
	print(i)
}

# load sp dis
rm(list=ls());gc()
library(dplyr)
load("data analysis/NativeRangeFragments.plants.ENA.clean.new.RData")	
load("data analysis/Domain.ena.rda")
rangeData=rangeData.native%>%left_join(domain.ena[,c("x","y","domain.ID")],by="domain.ID")%>%dplyr::select(Accepted_name,x,y)%>%split(.,.[,1])		
toload=list.files("/blue/matthewthomas1/yunpeng.liu/data analysis/future climate/")%>%.[grep(".rda",.)]
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
	load(paste0("/blue/matthewthomas1/yunpeng.liu/data analysis/future climate/",i))
	#rangeData=	rangeData[8041:10616]
	print(i)
	seql=c(seq(1,length(rangeData),no_cores*67),length(rangeData))#30
	seql2=lapply(2:length(seql),function(i,seql){if (i<length(seql)){c(seql[i-1],(seql[i]-1))
		} else {
		c(seql[i-1],seql[i])
		}},seql)
	if (i==toload[1]) nstart=3 else nstart=1	
	if (i==toload[7]) nend=2 else nend=length(seql2)
	#nstart=1	
	for (j in nstart:nend){
		mycl <- makePSOCKcluster(no_cores);		
		p_df=parLapply(cl=mycl,rangeData[seql2[[j]][1]:seql2[[j]][2]],get.maxent,re[[1]])
		#p_df=parLapply(cl=mycl,rangeData,get.maxent,re[[1]])
		stopCluster(mycl)
		save(p_df,file=paste0("/blue/matthewthomas1/yunpeng.liu/data analysis/SDM/sdm_",seql2[[j]][2],"_",i));
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
	
###  combine all the predicted distribution within each clim scenerio ------
# calculate threshold and converted the probability of occurrence of each species × grid-cell combination into binary maps. 
# We adopted the threshold that maximized the specificity and the sensitivity, using the rocr package (Sing et al., 2005; v.1.0-11). 
# If the probability of a species occurrence in a given grid cell was higher than this threshold, then the species was considered to be present in that grid cell. 
# ref:shijia new phy, https://github.com/Shijia818/Phenology-informed-SDM/blob/main/RCode.R
library(ROCR) 
current=list.files("/blue/matthewthomas1/yunpeng.liu/data analysis/SDM/")%>%.[grep("current",.)]
p_current=c()
for (i in current){
load(paste0("/blue/matthewthomas1/yunpeng.liu/data analysis/SDM/",i))
p_current=c(p_current,p_df)
}
#load("data analysis/sdm_2680_ACCESS-CM2_ssp245_2061-2080.rda")
spnames <- names(p_current)
thres=c()
for (i in spnames){	
	rangeData[[i]]$dis=1
	dat=p_current[[i]]%>%left_join(rangeData[[i]],by=c("x","y"))
	dat[is.na(dat)]=0
	if (sum(dat$predic-dat$dis)==0){
		thres=rbind(thres,data.frame(Species=i,thres=0))
	}else{
		pred <- ROCR::prediction(dat$predic,dat$dis)
		ss <- ROCR::performance(pred, "sens", "spec")
		threshold <- ss@alpha.values[[1]][is.finite(ss@alpha.values[[1]])][which.max(ss@x.values[[1]]+ss@y.values[[1]])]
		thres=rbind(thres,data.frame(Species=i,thres=threshold)) 
	}
	#print(which(spnames==i))
}
save(thres,file="data analysis/SDM/thres.rda")

# combine all the predicted distribution within each clim scenerio
load("data analysis/Domain.ena.rda")
load("/blue/matthewthomas1/yunpeng.liu/data analysis/SDM/thres.rda")
require(dplyr)
# clim.model=c("ACCESS-CM2","EC-Earth3-Veg","MPI-ESM1-2-HR","UKESM1-0-LL")
# period=c("2021-2040","2041-2060","2061-2080")
# ssp=c("ssp126","ssp245","ssp370","ssp585")
file.list=list.files("/blue/matthewthomas1/yunpeng.liu/data analysis/future climate/")%>%.[grep(".rda",.)]%>%gsub(".rda","",.)
#cal.comm=function(name,domain.ena){
for (name in file.list){
	#name=file.list[1]
	require(dplyr)
	get.bidis=function(p_df,thres,domain.ena){
		thes=thres[thres$Species%in%unique(p_df$Species),"thres"]
		p_df2=p_df[p_df$predic>=thes,]
		if (nrow(p_df2)==0) {p_df2=p_df;thes=0}
		p_df2=p_df2%>%left_join(domain.ena,by=c("x","y"))%>%na.omit()%>%select(Species,domain.ID,,x,y,predic)%>%distinct()		
		if (thes>0) p_df2$predic=1		
		colnames(p_df2)[1]="Tip_Label"		
		return(p_df2)
	}
	spdis.list=list.files("/blue/matthewthomas1/yunpeng.liu/data analysis/SDM/")%>%.[grep(name,.)]		
	spdis=c()
	for (i in spdis.list){
		spdis0=get(load(paste0("/blue/matthewthomas1/yunpeng.liu/data analysis/SDM/",i)))
		spdis.t=do.call(rbind,lapply(spdis0,get.bidis,thres,domain.ena))
		spdis=rbind(spdis,spdis.t)
	}
	save(spdis,file=paste0("/blue/matthewthomas1/yunpeng.liu/data analysis/SDM/sdm.spdis_",name,".rda"))
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
get.sim.sp=function(name){
	#for(name in file.list){
	require(dplyr)
	nCells=8490	
	rangeData=data.frame(domain.ID=1:nCells)%>%left_join(get(load(paste0("/blue/matthewthomas1/yunpeng.liu/data analysis/SDM/sdm.spdis_",name,".rda"))),by="domain.ID")
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
	save(CommSimilarity,file=paste0("/blue/matthewthomas1/yunpeng.liu/data analysis/SDM/CommSimilarity_",name,".rda"))
	rm(speciesByCell,predictByCell);gc()
	#}	
}
no_cores <-45#detectCores() - 1
mycl <- makePSOCKcluster(no_cores);	
parSapply(cl=mycl,file.list,get.sim.sp)
stopCluster(mycl)

########################################################################################################
### calculate expected alien range size for different ER thresholds under future climate and obtain the best ERQ
########################################################################################################
library(dplyr);
load("data analysis/connectedGraph.ena.rda")
#inputER<-c(0.1,0.2,0.3,seq(0.4,1,0.01))
#inputER<-c(0.3,seq(0.35,0.65,0.01),0.7)
inputER<-c(0.4,0.5,0.52,0.57)
#RangeFragments<-get(load("data analysis/AlienRangeFragments.plants.ENA.RData"))
load("data analysis/AlienRangeFragments.alltaxa.ENA.new.RData")
splist=read.csv("data analysis/US-RIIS.csv")%>%select(kingdom,phylum,class,scientificName)%>%filter(scientificName%in%names(RangeFragments.v2))
taxa=c("Chordata","Arthropoda","Other animal","Tracheophyta")#other animal incl. Mollusca,Annelida,Nematoda, etc.

get.simrange=function(x,RangeFragments,simmat,g1cc,mat.name){
	require(dplyr)
	fragmentTable1=data.frame(fragBinomial=unique(names(RangeFragments)),fragSpecies=1:length(unique(names(RangeFragments))))
	fragmentTable2=data.frame(fragNum=1:length(RangeFragments),fragSize=unlist(lapply(RangeFragments,length)),fragBinomial=names(RangeFragments))
	FragmentTable=fragmentTable2%>%left_join(fragmentTable1,by="fragBinomial")
	#for(x in 1:inputER){  
	source("data analysis/Functions.R")
	output<-runConnectedComponentsModel(threshold=x,nRandomStarts=500,fragmentTable=FragmentTable,RangeFragments=RangeFragments,CommSimilarity=simmat,g1cc=g1cc)#this function comes from Lovell et al. (2021).,Nature Ecology & Evolution, 5(3), 322-329
	#save the simulated ranges
	fragmentTable<-output[[1]]
	simRanges<-output[[2]]
	save(simRanges,fragmentTable,file=paste0("/blue/matthewthomas1/yunpeng.liu/data analysis/SDM.simulation.ena.predict/PredictedRangeSize_",x,"_",mat.name,".rda")) 
	rm(fragmentTable1,fragmentTable2,output,RangeFragments,fragmentTable,simRanges);gc()
	}
	
toload<-list.files("/blue/matthewthomas1/yunpeng.liu/data analysis/SDM/")%>%.[grep("CommSimilarity_",.)]
library(parallel)
no_cores <- 45#detectCores() - 1
for (mat in toload[-c(1:13)]){#length(toload)){
	mycl <- makePSOCKcluster(no_cores);
	simmat=get(load(paste0("/blue/matthewthomas1/yunpeng.liu/data analysis/SDM/",mat)))
	print(mat)
	mat.name=mat%>%gsub("CommSimilarity_","",.)%>%gsub(".rda","",.)
	parSapply(cl=mycl,X=inputER,get.simrange,RangeFragments.v2,simmat,g1cc,mat.name)	
	stopCluster(mycl)
	rm(simmat);gc()	
}
########################################################################################################
## generate matrix of invasion probabilities using the best ER model #### 
########################################################################################################
library(dplyr);
load("data analysis/Domain.ena.rda")
load("data analysis/AlienRangeFragments.alltaxa.ENA.new.RData")
splist=read.csv("data analysis/US-RIIS.csv")%>%select(kingdom,phylum,class,scientificName)%>%filter(scientificName%in%names(RangeFragments.v2))
rangesize0=do.call(c,lapply(RangeFragments.v2,length))
rangesize0=data.frame(sp=names(rangesize0),range=rangesize0)%>%group_by(sp)%>%summarize(rangesize=sum(range))
rangesize=rangesize0$rangesize;names(rangesize)=rangesize0$sp
taxa=c("Chordata","Arthropoda","Other animal","Tracheophyta")	
# load("data analysis/simulation.ena.predict2/BestERvalues.sp.rda")
# stat.ERQ=input%>%group_by(taxa,tau)%>%reframe(inputER=inputER[which(rmseQspecies==min(rmseQspecies))],underpred=underpredicted[which(rmseQspecies==min(rmseQspecies))])%>%
	# group_by(taxa)%>%reframe(tau.opt=tau[which(underpred==min(underpred))],bestERQ=inputER[which(underpred==min(underpred))],underpred.por=min(underpred))%>%as.data.frame
stat.ERQ=data.frame(taxa=c("Chordata","Arthropoda","Other animal","Tracheophyta"),bestERQ=c(0.4,0.57,0.52,0.5))	
get.InvasionProb=function(mat.name,taxa,splist,stat.ERQ,domain.ena,RangeFragments.ena.inv,current=TRUE){
	nRandomStarts<-500
	InvasionProbabilityMatrix=list()
	for (taxa.i in taxa){
		sp=splist[splist$phylum==taxa.i,"scientificName"]
		if (taxa.i%in%"Other animal") sp=splist[splist$kingdom=="Animalia"&(!splist$phylum%in%c("Chordata","Arthropoda")),"scientificName"]	
		nSpecies<-length(sp)
		bestERQ=stat.ERQ[stat.ERQ$taxa==taxa.i,"bestERQ"]
		if (current)load(paste0("/blue/matthewthomas1/yunpeng.liu/data analysis/simulation.ena.predict2/PredictedRangeSize_",bestERQ,".",mat.name,".rda"))
		if (!current)load(paste0("/blue/matthewthomas1/yunpeng.liu/data analysis/SDM.simulation.ena.predict/PredictedRangeSize_",bestERQ,"_",mat.name,".rda"))	
		InvasionProbabilityMatrix[[taxa.i]]<-matrix(0,ncol=nrow(domain.ena),nrow=nSpecies);
		colnames(InvasionProbabilityMatrix[[taxa.i]])=domain.ena$domain.ID;
		rownames(InvasionProbabilityMatrix[[taxa.i]])=sp	
		for(i in 1:nSpecies){ # loop over each species  
		  focRowAlien<-which(fragmentTable$fragBinomial==sp[i])
		  nFrags<-length(focRowAlien)
		   # for each alien range fragment calculate the probability of invading each cell
		  InvasionProbs<-matrix(0,ncol=nrow(domain.ena),nrow=nFrags)
		  for(x in 1:nFrags){
			fragSize<-fragmentTable$fragSize[focRowAlien[x]]
			InvasionProbs[x,as.numeric(names(simRanges[[focRowAlien[x]]]))]<-as.numeric(simRanges[[focRowAlien[x]]])/min(nRandomStarts,fragSize)    
		  }	  
		  # the overall probability of an alien species invading a cell (a cell can be invaded from multiple fragments)
		  InvasionProbVec<-1-matrixStats::colProds(1-InvasionProbs)
		  
		  # set invasion probabilities of cells that are already invaded to 1
		  cells=subset(RangeFragments.ena.inv,names(RangeFragments.ena.inv)==sp[i])%>%unlist%>%unique%>%as.numeric%>%sort
		  InvasionProbVec[cells]<-1	  
		  InvasionProbabilityMatrix[[taxa.i]][i,]<-InvasionProbVec	  
		}	
	}
	return(InvasionProbabilityMatrix)
}		
## caculate number of range expand in each climate scenerio ---
InvasionProbabilityMatrix.current.sdm<-get.InvasionProb("current",taxa,splist,stat.ERQ,domain.ena,RangeFragments.v2,current=FALSE)
extent.current.sdm=lapply(InvasionProbabilityMatrix.current.sdm,rowSums)
#toload<-list.files("/blue/matthewthomas1/yunpeng.liu/data analysis/SDM/")%>%.[grep("CommSimilarity_",.)]

lab=c("ssp126_2021-2040","ssp245_2021-2040","ssp370_2021-2040","ssp585_2021-2040",
	"ssp126_2041-2060","ssp245_2041-2060","ssp370_2041-2060","ssp585_2041-2060",
		"ssp126_2061-2080","ssp245_2061-2080","ssp370_2061-2080","ssp585_2061-2080")
select.clim=lab[-c(1,5,9)]#c("ssp126_2021-2040","ssp370_2021-2040","ssp585_2021-2040")
model=c("EC-Earth3-Veg","MPI-ESM1-2-HR")
model.list=c(paste(rep(model[1],length(lab)),select.clim,sep="_"),paste(rep(model[2],length(lab)),select.clim,sep="_"))

range.stat.all=c()
for (mat.name in model.list){#length(toload)){
	#mat.name=mat%>%gsub(".rda","",.)%>%gsub("CommSimilarity_","",.)	
	InvasionProbabilityMatrix.future<-get.InvasionProb(mat.name,taxa,splist,stat.ERQ,domain.ena,RangeFragments.ena.inv=RangeFragments.v2,current=FALSE)
	extent.future.sdm=lapply(InvasionProbabilityMatrix.future,rowSums)
	for (i in 1:4) {
		range.diff=extent.future.sdm[[i]]-extent.current.sdm[[i]]
		range.sign=ifelse(range.diff>0,"expand",ifelse(range.diff<0,"shrink",NA))
		range.compare=data.frame(ranges=rangesize[names(range.diff)],range.diff,range.sign)%>%na.omit
		range.compare$range.sign=factor(range.compare$range.sign,levels=c("expand","shrink"))
		comp=t.test(ranges~range.sign,data=range.compare)
		range.stat=range.compare%>%group_by(range.sign)%>%summarize(mean=mean(ranges),se=sd(ranges)/sqrt(length(ranges)),Richness=length(ranges))%>%
			cbind(data.frame(scenerio=mat.name,taxa=names(extent.future.sdm)[i],t=comp$statistic,df=comp$parameter,p=comp$p.value),.)
		range.stat.all=rbind(range.stat.all,range.stat)			
	}	
}
save(range.stat.all,file="/blue/matthewthomas1/yunpeng.liu/data analysis/range.stat2.rda")
#Plot
library(ggpubr);library(dplyr)
load("data analysis/range.stat2.rda")
head(range.stat.all)
model.list2=paste(rep(model[1],length(lab)),select.clim,sep="_")
ranges=range.stat.all%>%filter(scenerio%in%model.list2)%>% select(scenerio,taxa,range.sign,Richness)#%>% tidyr::pivot_longer(cols=expand:shrink,names_to = "Type", values_to = "Richness")
ranges$scenerio=ranges$scenerio%>%gsub(paste0(model[1],"_"),"",.)
range.stat=	range.stat.all %>% filter(scenerio%in%model.list2) %>% select(scenerio,taxa,range.sign,mean,se)
range.stat$scenerio=range.stat$scenerio%>%gsub(paste0(model[1],"_"),"",.)

p1=ggplot(ranges,aes(y=scenerio, x=Richness,fill=range.sign)) +
    geom_bar(position="dodge", stat="identity")+
	labs(y="Future climate scenerio",x="Number of species",fill="Range size change")+
	facet_wrap(~taxa, scales="free_x",nrow=1)+
	theme_classic()+theme(axis.text = element_text(size=12,color='black'),
		#axis.text.x =element_text(angle=60),
		axis.title = element_text(size=12,color='black'),
		#strip.text = element_blank(),
		strip.text = element_text(colour = 'black', face="italic",size = rel(1.2)), 
		strip.background = element_rect(fill = 'white', colour = 'darkgray', size = rel(2), linetype = 2),
		legend.text=element_text(size=12),
		legend.title=element_text(face="bold",size=12))
p2=ggplot(range.stat,aes(y=scenerio, x=mean,fill=range.sign)) +
    geom_point(shape=21,color="black",size=3,show.legend=F) +
	geom_errorbar(aes(y=scenerio, xmin=mean-se, xmax=mean+se,color=range.sign), width=0.4, alpha=0.8, size=1,show.legend=F)+
	labs(y="Future climate scenerio",x="Realized range sizes",fill="Range size change")+
	facet_wrap(~taxa, scales="free_x",nrow=1)+
	theme_classic()+theme(axis.text = element_text(size=12,color='black'),
		#axis.text.x =element_text(angle=60),
		axis.title = element_text(size=12,color='black'),
		#strip.text = element_blank(),
		strip.text = element_text(colour = 'black', face="italic",size = rel(1.2)), 
		strip.background = element_rect(fill = 'white', colour = 'darkgray', size = rel(2), linetype = 2),
		legend.text=element_text(size=12),
		legend.title=element_text(face="bold",size=12))
ggarrange(p1,p2,ncol=1,align="v",labels="auto",hjust=0,vjust=0.8)
	
## compare range extent of each species ---
#current invasion extent
InvasionProbabilityMatrix.current<-get.InvasionProb("sp",taxa,splist,stat.ERQ,domain.ena,RangeFragments.v2)
extent.current=lapply(InvasionProbabilityMatrix.current,rowSums)
#current invasion extent based on SDM
InvasionProbabilityMatrix.current.sdm<-get.InvasionProb("current",taxa,splist,stat.ERQ,domain.ena,RangeFragments.v2,current=FALSE)
extent.current.sdm=lapply(InvasionProbabilityMatrix.current,rowSums)

#toload<-list.files("/blue/matthewthomas1/yunpeng.liu/data analysis/SDM/")%>%.[grep("CommSimilarity_",.)]
for (mat in toload){#length(toload)){
	mat.name=mat%>%gsub(".rda","",.)%>%gsub("CommSimilarity_","",.)	
	InvasionProbabilityMatrix.future<-get.InvasionProb(mat.name,taxa,splist,stat.ERQ,domain.ena,RangeFragments.ena.inv=RangeFragments.v2,current=FALSE)
	InvasionProb<-do.call(cbind,lapply(InvasionProbabilityMatrix.future,colSums))

	#compare range extent of each species
	range.diff=list()
	for (taxa.i in taxa){
		range.current=InvasionProbabilityMatrix.current[[taxa.i]]
		range.future=InvasionProbabilityMatrix.future[[taxa.i]]
		range.diff[[taxa.i]]=(range.future-range.current)	
	}
	Invasion.shift.grid<-do.call(cbind,lapply(range.diff,colSums))
	
	#sp level change
	get.shift=function(i,range.diff,extent.current) {
		shifts=rowSums(range.diff[[i]])/extent.current[[i]]*100	
		se=function(x) sd(x)/sqrt(length(x))
		return(data.frame(shifts.mean=mean(shifts),shifts.se=se(shifts)))
	}
	InvasionProb.shift.sp=do.call(rbind,lapply(taxa,get.shift,range.diff,extent.current))
	InvasionProb.shift.sp$taxa=taxa;
	InvasionProb.shift.sp$clim.scenerio=mat.name%>%gsub("ACCESS-CM2_","",.)
	Invasion.shift.sp=InvasionProb.shift.sp
	
	#current based on SDM
	range.diff.sdm=list()
	for (taxa.i in taxa){
		range.current=InvasionProbabilityMatrix.current.sdm[[taxa.i]]
		range.future=InvasionProbabilityMatrix.future[[taxa.i]]
		range.diff.sdm[[taxa.i]]=(range.future-range.current)	
	}
	Invasion.shift.grid.sdm<-do.call(cbind,lapply(range.diff.sdm,colSums))
	
	InvasionProb.shift.sp.sdm=do.call(rbind,lapply(taxa,get.shift,range.diff.sdm,extent.current.sdm))
	InvasionProb.shift.sp.sdm$taxa=taxa;
	InvasionProb.shift.sp.sdm$clim.scenerio=mat.name%>%gsub("ACCESS-CM2_","",.)
	Invasion.shift.sp.sdm=InvasionProb.shift.sp.sdm
	
	
	output=list(Invasion.shift.grid,Invasion.shift.sp,Invasion.shift.grid.sdm,Invasion.shift.sp.sdm,InvasionProb)
	save(output,file=paste0("data analysis/Invasion.shift/Invasion.shift.",mat.name,".rda"))
	print(mat)
	# #plot future invasion richness
	# library(ggpubr)
	# the=theme_bw()+theme(legend.background=element_rect(fill='transparent'),
			# legend.text=element_text(size=10),
			# legend.title=element_text(face="bold",size=10))	
	# InvasionProb<-do.call(cbind,lapply(InvasionProbabilityMatrix,colSums))
	# plt=list()
	# for (taxa.i in taxa){
	# domain=cbind(domain.ena[,c("x","y")],InvasionProb[,taxa.i])
	# colnames(domain)[3]="InvasionProb"
	# plt[[taxa.i]]=ggplot() +  geom_tile(data = domain, aes(x = x, y = y, fill = InvasionProb)) +
	  # coord_quickmap() +  labs(x="",y="",fill="Potential \ninvasion\nrichness")+
	  # annotate("text",x=-72,y=30,size=5,label=taxa.i,color="black")+
	  # scale_fill_gradientn(colours = viridis::viridis(99),na.value = "lightgray")+the#+theme(legend.position=c(0.8,0.25))
	# }
	# ggarrange(plt[[1]],plt[[2]],plt[[3]],plt[[4]],labels="auto",font.label = list(size = 20))	
}
		
library(dplyr)	
toload<-list.files("data analysis/Invasion.shift/")
range.stat=c();range.stat.sdm=c()
for (mat in toload){
	load(paste0("data analysis/Invasion.shift/",mat))
	model=strsplit(mat,"_")[[1]][1] %>% gsub("Invasion.shift.","",.)
	range.stat=rbind(range.stat,data.frame(model=model,output[[2]]))
	range.stat.sdm=rbind(range.stat.sdm,data.frame(model=model,output[[4]]))
}
save(range.stat,file="data analysis/Invasion.shift/range.stat.rda")
save(range.stat.sdm,file="data analysis/Invasion.shift/range.stat.sdm.rda")

## plot species level changes
library(ggpubr);library(dplyr)
#load("data analysis/Invasion.shift/range.stat.rda")
load("data analysis/Invasion.shift/range.stat.sdm.rda")
model.list=unique(range.stat.sdm$model)
lab=c("ssp126_2021-2040","ssp245_2021-2040","ssp370_2021-2040","ssp585_2021-2040",
	"ssp126_2041-2060","ssp245_2041-2060","ssp370_2041-2060","ssp585_2041-2060",
		"ssp126_2061-2080","ssp245_2061-2080","ssp370_2061-2080","ssp585_2061-2080")
select.clim=lab[-c(1,5,9)]#c("ssp126_2021-2040","ssp370_2021-2040","ssp585_2021-2040")
p=list()
for (i in 1:length(model.list)){
dat=range.stat.sdm%>%filter(model==model.list[i])
dat$clim.scenerio=gsub(paste0(model.list[i],"_"),"",dat$clim.scenerio)
dat$clim.scenerio=factor(dat$clim.scenerio,levels=lab)
dat=dat%>%filter(clim.scenerio%in%select.clim)

p[[i]]=ggplot(dat, aes(x=clim.scenerio, y=shifts.mean)) + 	
	geom_errorbar(aes(ymin=shifts.mean-shifts.se, ymax=shifts.mean+shifts.se),linewidth=1,width=0.2,show.legend=FALSE,alpha=0.6)+
	# geom_vline(xintercept=4.5,col="lightgray",linewidth=1,linetype="longdash")+
	# geom_vline(xintercept=8.5,col="lightgray",linewidth=1,linetype="longdash")+
	geom_hline(yintercept=0,col="#F8766D",linewidth=1,linetype="longdash")+
	geom_point(size=3,shape=21,color="black",fill="#619CFF") +coord_flip() +	
	theme_bw()+
	labs(y="Shift in predicted extent\n(Future - present, %)",x="Climate change scenerios")+ 
	facet_wrap(~taxa,scales="fixed",nrow=1)+
	theme(axis.text = element_text(size=10,color='black',angle=0,hjust=1),
		axis.title = element_text(size=15,color='black',angle=0),		
		strip.text=element_text(size=12,color='black',face="italic"),
		strip.background = element_rect(fill = 'white', colour = 'black', size = rel(2), linetype = 2))
}
p[[2]]
p[[3]]
#plot range shift under a certain clim scenerio
the=theme_bw()+theme(legend.background=element_rect(fill='transparent'),
			legend.text=element_text(size=10),
			legend.title=element_text(face="bold",size=10))		
taxa=c("Chordata","Arthropoda","Other animal","Tracheophyta")
load("data analysis/Domain.ena.rda")
load(paste0("data analysis/Invasion.shift/",toload[16]))

plt=c()
for (taxa.i in taxa){
	domain=cbind(domain.ena[,c("x","y")],output[[3]][,taxa.i])
	colnames(domain)[3]="InvasionProb"
	plt[[taxa.i]]=ggplot() +  geom_tile(data = domain, aes(x = x, y = y, fill = InvasionProb)) +
	  coord_quickmap() +  labs(x="",y="",fill="Changes in\npotential\ninvasion\nrichness")+
	   annotate("text",x=-72,y=30,size=5,label=taxa.i,color="black")+
	  scale_fill_gradient2(low = "red",mid = "white",high = "blue",midpoint = 0, na.value = "darkgray")+the+
	theme(panel.background = element_rect(fill = '#619CFF', colour = 'black'))
}

#plot invasion extent under a certain clim scenerio
plt=list()
for (taxa.i in taxa){
domain=cbind(domain.ena[,c("x","y")],output[[5]][,taxa.i])
colnames(domain)[3]="InvasionProb"
plt[[taxa.i]]=ggplot() +  geom_tile(data = domain, aes(x = x, y = y, fill = InvasionProb)) +
  coord_quickmap() +  labs(x="",y="",fill="Potential \ninvasion\nrichness")+
  annotate("text",x=-72,y=30,size=5,label=taxa.i,color="black")+
  scale_fill_gradientn(colours = viridis::viridis(99),na.value = "lightgray")+the#+theme(legend.position=c(0.8,0.25))
}
ggarrange(plt[[1]],plt[[2]],plt[[3]],plt[[4]],labels="auto",font.label = list(size = 20))	
		
#plot richness
library(ggpubr)
load(paste0("data analysis/Invasion.shift/",toload[16]))

the=theme_bw()+theme(legend.background=element_rect(fill='transparent'),
		legend.text=element_text(size=10),
		legend.title=element_text(face="bold",size=10))		
plt=list()
for (taxa.i in taxa){
domain=cbind(domain.ena[,c("x","y")],output[[5]][,taxa.i])
colnames(domain)[3]="InvasionProb"
plt[[taxa.i]]=ggplot() +  geom_tile(data = domain, aes(x = x, y = y, fill = InvasionProb)) +
  coord_quickmap() +  labs(x="",y="",fill="Future\ninvasion\nrichness")+
 annotate("text",x=-73,y=30,size=4,label=taxa.i,color="black")+
  scale_fill_gradientn(colours = viridis::viridis(99),na.value = "lightgray")+the#theme(legend.position=c(0.8,0.25))
  
}
p2=ggarrange(plt[[1]],plt[[2]],plt[[3]],plt[[4]],labels=c("b","c","d","e"),font.label = list(size = 20))	
ggarrange(p[[2]],p2,labels=c("a",""),font.label = list(size = 20))	

#	plot richness v2
EC=toload[16:24]
MP=toload[28:36]
plt=list()
for (taxa.i in taxa){
for (clim in EC){
	name=paste0(taxa.i,clim)
	load(paste0("data analysis/Invasion.shift/",clim))
	domain=cbind(domain.ena[,c("x","y")],output[[5]][,taxa.i])
	colnames(domain)[3]="InvasionProb"
	plt[[name]]=ggplot() +  geom_tile(data = domain, aes(x = x, y = y, fill = InvasionProb)) +
	  coord_quickmap() +  labs(x="",y="",fill="Future invasion richness")+
	 #annotate("text",x=-73,y=30,size=4,label=taxa.i,color="black")+
	  scale_fill_gradientn(colours = viridis::viridis(99),na.value = "lightgray")+the#theme(legend.position=c(0.8,0.25))
} 
}
ggarrange(plt[[1]],plt[[2]],plt[[3]],plt[[4]],plt[[5]],plt[[6]],plt[[7]],plt[[8]],plt[[9]],labels="auto",font.label = list(size = 20),common.legend=T)
ggarrange(plt[[10]],plt[[11]],plt[[12]],plt[[13]],plt[[14]],plt[[15]],plt[[16]],plt[[17]],plt[[18]],labels="auto",font.label = list(size = 20),common.legend=T)
ggarrange(plt[[19]],plt[[20]],plt[[21]],plt[[22]],plt[[23]],plt[[24]],plt[[25]],plt[[26]],plt[[27]],labels="auto",font.label = list(size = 20),common.legend=T)
ggarrange(plt[[28]],plt[[29]],plt[[30]],plt[[31]],plt[[32]],plt[[33]],plt[[34]],plt[[35]],plt[[36]],labels="auto",font.label = list(size = 20),common.legend=T)