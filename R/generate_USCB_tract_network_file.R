###disable scientific notation###
options(scipen = 999)

###load packages###
library(data.table)

###for reading in dbf###
library(foreign)
library(sf)
library(censusapi)
library(igraph)
#library(censusxy)
library(jsonlite)
data.table::setDTthreads(1)


generate_USCB_tract_network_file <- function(FIPS_dt, USCB_TIGER.path, omit.park_openspace=TRUE, omit.unpopulated=TRUE, use.bridges=TRUE, geo.year="2010", ADDR_dt=NULL, output.path=NULL, output.file_name=NULL){

	if(as.character(geo.year)=="2010") {
		f_year <- "2019"
	} else {
		f_year <- "2022"
		geo.year <- '2020'
	}
	
	FIPS.dt <- copy(as.data.table(FIPS_dt))

	###pad state and county codes with leading zeros###
	FIPS.dt[,state := sprintf("%02d", as.numeric(state))]
	FIPS.dt[,county := sprintf("%03d", as.numeric(county))]
	
	FIPS.dt <- unique(FIPS.dt[,c("state","county"),with=FALSE])
	
	old.wd <- getwd()

	setwd(USCB_TIGER.path)
	
	
	#######################################################################
	###consolidate county-level edge shapefiles into a single data.table###
	#######################################################################
	
	edge.files <- paste0("tl_",f_year,"_",FIPS.dt$state,FIPS.dt$county,"_edges")
	
	###check if files are where they should be###
	ff <- list.files(file.path(USCB_TIGER.path,"EDGES"))
	
	if(length(edge.files[!edge.files %in% ff]) > 0){
		stop("\nThe following files are missing:\n",paste(edge.files[!edge.files %in% ff],collapse="\n"))
	}
	
	edges.dt <- rbindlist(lapply(edge.files,function(j){
	
		###read in shapefile###
		temp.sf <- sf::st_read(file.path("EDGES",j), j, stringsAsFactors = F, quiet=T)
		temp.sf$edge_len <- st_length(temp.sf)
		
		###convert to data.table###
		return(as.data.table(st_drop_geometry(temp.sf)))

	}), use.names=TRUE, fill=TRUE)
	
	#######################################################################
	###consolidate county-level face shapefiles into a single data.table###
	#######################################################################
	
	face.files <- paste0("tl_",f_year,"_",FIPS.dt$state,FIPS.dt$county,"_faces")
	
	###check if files are where they should be###
	ff <- list.files(file.path(USCB_TIGER.path,"FACES"))
	
	if(length(face.files[!face.files %in% ff]) > 0){
		stop("\nThe following files are missing:\n",paste(face.files[!face.files %in% ff],collapse="\n"))
	}
	
	faces.dt <- rbindlist(lapply(face.files,function(j){

		###read in shapefile###
		temp.sf <- sf::st_read(file.path("FACES",j), j, stringsAsFactors = F, quiet=T)
		temp.sf$poly_area <- as.numeric(st_area(temp.sf))
		
		###convert to data.table###
		return(as.data.table(st_drop_geometry(temp.sf)))

	}), use.names=TRUE, fill=TRUE)
	
	
	if(as.character(geo.year)=="2010") {
		faces.dt[,USCB_tract_this := paste0(STATEFP10,COUNTYFP10,TRACTCE10)]
		
		faces.dt[,USCB_block_this := paste0(STATEFP10,COUNTYFP10,TRACTCE10,BLOCKCE10)]
	} else {
		faces.dt[,USCB_tract_this := paste0(STATEFP20,COUNTYFP20,TRACTCE20)]
		
		faces.dt[,USCB_block_this := paste0(STATEFP20,COUNTYFP20,TRACTCE20,BLOCKCE20)]
	}
	
	######################################
	###pull block-level population data###
	######################################
	
	mycensuskey <-"2ca0b2830ae4835905efab6c35f8cd2b3f570a8a"
	
	if(as.character(geo.year)=="2010") {
		my.survey <- "dec/sf1"
		my.vars <- c("P001001")
		my.vintage <- 2010
	} else {
		my.survey <- "dec/pl"
		my.vars <- c("P1_001N")
		my.vintage <- 2020
	}

	api.data_cb <- rbindlist(lapply(1:nrow(FIPS.dt), function(x) as.data.table(getCensus(name = my.survey,
		vintage = my.vintage,
		key = mycensuskey,
		vars = my.vars,
		region = "block:*",
		regionin = paste0("state:",FIPS.dt[x]$state,"+county:",FIPS.dt[x]$county)))))

	api.data_cb[,USCB_block_this := paste0(state,county,tract,block)]
	api.data_cb[,USCB_tract_this := paste0(state,county,tract)]
	setnames(api.data_cb,my.vars,c("USCB_pop"))
	
	
	##################################################
	###load Topological Faces / Area Landmark files###
	##################################################
	
	facesal.files <- paste0("tl_",f_year,"_",unique(FIPS.dt$state),"_facesal")
	
	###check if files are where they should be###
	ff <- list.files(file.path(USCB_TIGER.path,"FACESAL"))
	
	if(length(facesal.files[!facesal.files %in% ff]) > 0){
		stop("\nThe following files are missing:\n",paste(facesal.files[!facesal.files %in% ff],collapse="\n"))
	}
	
	
	facesal.dt <- rbindlist(lapply(facesal.files,function(j){

		###read in dbf###
		return(as.data.table(read.dbf(file.path(file.path("FACESAL",j),paste0(j,'.dbf')))))

	}), use.names=TRUE, fill=TRUE)
	
	
	##############################
	###load Area Landmark files###
	##############################
	
	arealm.files <- paste0("tl_",f_year,"_",unique(FIPS.dt$state),"_arealm")
	
	###check if files are where they should be###
	ff <- list.files(file.path(USCB_TIGER.path,"AREALM"))
	
	if(length(arealm.files[!arealm.files %in% ff]) > 0){
		stop("\nThe following files are missing:\n",paste(arealm.files[!arealm.files %in% ff],collapse="\n"))
	}
	
	arealm.dt <- rbindlist(lapply(arealm.files,function(j){

		###read in shapefile###
		temp.sf <- sf::st_read(file.path("AREALM",j), j, stringsAsFactors = F, quiet=T)
		#temp.sf$poly_area <- as.numeric(st_area(temp.sf))
		
		###convert to data.table###
		return(as.data.table(st_drop_geometry(temp.sf)))

	}), use.names=TRUE, fill=TRUE)
	
	
	###restrict landmarks list to parks, zoos, cemeteries and other unpopulated open spaces###
	p.os <- c("K2180", "K2181", "K2183", "K2184", "K2185", "K2186", "K2187", "K2188", "K2189", "K2190", "K2582", "K2586", "K2564", "K2451", "K2456", "K2457", "K2561", "K2564")
	
	arealm.dt <- arealm.dt[(MTFCC %in% p.os)]
	
	###calculate proportion of census block that is park or open space###
	arealm.dt <- merge(arealm.dt, facesal.dt, by="AREAID")
	
	arealm.dt <- merge(arealm.dt, faces.dt[LWFLAG !="P", c("TFID","USCB_tract_this","USCB_block_this","poly_area"),with=FALSE], by="TFID")
	
	arealm.dt <- merge(arealm.dt[,.(poly_area=sum(poly_area)),by=list(USCB_block_this)], faces.dt[LWFLAG != "P",.(poly_area=sum(poly_area)), by=USCB_block_this], by="USCB_block_this")
	
	arealm.dt[,prop.area := round(poly_area.x/poly_area.y,2)]

	arealm.dt <- merge(arealm.dt, api.data_cb[,.(USCB_pop=sum(USCB_pop)),by=list(USCB_block_this)], by="USCB_block_this",all.x=TRUE)
	
	#############################################
	###generate table of census blocks to omit###
	#############################################
	
	omit.dt <- merge(faces.dt[LWFLAG !="P", .(poly_area=sum(poly_area)),by=list(USCB_block_this,USCB_tract_this)], api.data_cb[,c("USCB_block_this","USCB_pop"), with=FALSE], by="USCB_block_this", all.x=TRUE)
	
	omit.dt[,USCB_pop := ifelse(is.na(USCB_pop),0,USCB_pop)]

	###deal with errors where unpopulated parks/open spaces have USCB_pop > 0 (e.g., Central Park, Bronx Zoo)###
	omit.dt[,USCB_pop := ifelse(USCB_block_this %in% unique(arealm.dt[USCB_pop > 0 & prop.area==1]$USCB_block_this), 0, USCB_pop)]
	
	omit.dt <- merge(omit.dt[USCB_pop==0,.(poly_area=sum(poly_area)),by=list(USCB_tract_this)],omit.dt[,.(poly_area=sum(poly_area)),by=list(USCB_tract_this)], by="USCB_tract_this", all.y=TRUE)

	omit.dt[,prop.area := round(poly_area.x/poly_area.y,2)]
	
	omit.dt <- omit.dt[prop.area==1]
	
	omit.dt[,type := ifelse(USCB_tract_this %in% unique(substr(arealm.dt[USCB_pop > 0 & prop.area==1]$USCB_block_this,1,11)),"park_openspace","unpopulated")]
	
	
	###################################################
	###generate a table of neighboring census blocks###
	###################################################
	
	b.dt <- merge(merge(unique(edges.dt[,c("TLID","TFIDL","TFIDR","edge_len"),with=FALSE]),unique(faces.dt[,c("TFID","USCB_block_this","LWFLAG"),with=FALSE]),by.x="TFIDL",by.y="TFID"), unique(faces.dt[,c("TFID","USCB_block_this","LWFLAG"),with=FALSE]),by.x="TFIDR",by.y="TFID")

	b.dt <- b.dt[USCB_block_this.x != USCB_block_this.y]

	b.dt[,USCB_block_this.1 := ifelse(as.numeric(USCB_block_this.x) > as.numeric(USCB_block_this.y), USCB_block_this.y, USCB_block_this.x)]
	b.dt[,LWFLAG.1 := ifelse(as.numeric(USCB_block_this.x) > as.numeric(USCB_block_this.y), LWFLAG.y, LWFLAG.x)]

	b.dt[,USCB_block_this.2 := ifelse(as.numeric(USCB_block_this.x) > as.numeric(USCB_block_this.y), USCB_block_this.x, USCB_block_this.y)]
	b.dt[,LWFLAG.2 := ifelse(as.numeric(USCB_block_this.x) > as.numeric(USCB_block_this.y), LWFLAG.x, LWFLAG.y)]

	b.dt2 <- b.dt[,.(edge_len=sum(edge_len)),by=list(USCB_block_this.1,LWFLAG.1,USCB_block_this.2,LWFLAG.2)]

	b.dt2 <- b.dt2[LWFLAG.1 != "P" & LWFLAG.2 != "P"]
	
	#############################################################
	###determine number of total number of neighbors per block###
	#############################################################

	###all relationships###
	zz1 <- b.dt2[,.(tot=.N),by=USCB_block_this.1]
	setnames(zz1,c("USCB_block_this.1"),c("USCB_block_this"))

	zz2 <- b.dt2[,.(tot=.N),by=USCB_block_this.2]
	setnames(zz2,c("USCB_block_this.2"),c("USCB_block_this"))

	block.dt1 <- rbindlist(list(zz1,zz2),use.names=TRUE,fill=TRUE)

	###inter county relationships###
	zz1 <- b.dt2[substr(USCB_block_this.1,1,5) != substr(USCB_block_this.2,1,5),.(tot.c=.N),by=USCB_block_this.1]
	setnames(zz1,c("USCB_block_this.1"),c("USCB_block_this"))

	zz2 <- b.dt2[substr(USCB_block_this.1,1,5) != substr(USCB_block_this.2,1,5),.(tot.c=.N),by=USCB_block_this.2]
	setnames(zz2,c("USCB_block_this.2"),c("USCB_block_this"))

	block.dt2 <- rbindlist(list(zz1,zz2),use.names=TRUE,fill=TRUE)

	block_rel.dt <- merge(block.dt1,block.dt2,by="USCB_block_this",all.x=TRUE,all.y=TRUE)

	rm(zz1,zz2,block.dt1,block.dt2)

	block_rel.dt[,tot.c := ifelse(is.na(tot.c),0,tot.c)]

	block_rel.dt <- merge(block_rel.dt, api.data_cb[,c("USCB_block_this","USCB_pop"), with=FALSE], by="USCB_block_this", all.x=TRUE)

	###remove piers misassigned to another county (e.g., piers in Queens and Brooklyn misassigned to Manhattan)###
	piers.dt <- block_rel.dt[tot.c==tot & USCB_pop==0]

	b.dt2[,is.pier := ifelse(USCB_block_this.1 %in% piers.dt$USCB_block_this | USCB_block_this.2 %in% piers.dt$USCB_block_this,1,0)]

	###generate tracts from blocks###
	b.dt2[,USCB_tract_this.1 := substr(USCB_block_this.1,1,11)]
	b.dt2[,USCB_tract_this.2 := substr(USCB_block_this.2,1,11)]

	###aggregate by census tracts and exclude piers###
	neighbors.dt <- b.dt2[is.pier==0,.(edge_len=sum(as.numeric(edge_len))),by=list(USCB_tract_this.1,USCB_tract_this.2)]
	neighbors.dt <- neighbors.dt[USCB_tract_this.1 != USCB_tract_this.2]
	
	############################
	###generate a node tables###
	############################

	node.dt1 <- edges.dt[,c("TNIDF","TFIDL","TFIDR"),with=FALSE]
	setnames(node.dt1,c("TNIDF"),c("TNID"))

	node.dt2 <- edges.dt[,c("TNIDT","TFIDL","TFIDR"),with=FALSE]
	setnames(node.dt2,c("TNIDT"),c("TNID"))

	node.dt <- unique(rbindlist(list(node.dt1,node.dt2),use.names=TRUE,fill=TRUE))

	node.dt <- merge(node.dt,faces.dt[,c("TFID", "USCB_block_this", "USCB_tract_this","LWFLAG"),with=FALSE],by.x="TFIDL",by.y="TFID",all.x=TRUE)

	node.dt <- merge(node.dt,faces.dt[,c("TFID", "USCB_block_this", "USCB_tract_this","LWFLAG"),with=FALSE],by.x="TFIDR",by.y="TFID",all.x=TRUE)

	node.dt.m <- melt(node.dt, id.vars = c("TNID"), measure = list(c("USCB_tract_this.x", "USCB_tract_this.y"), c("LWFLAG.x", "LWFLAG.y")), value.name = c("USCB_tract_this", "LWFLAG"))

	node.dt.m[,variable:=NULL]
	node.dt.all <- unique(copy(node.dt.m))

	###remove water tracts###
	node.dt.m <- node.dt.m[LWFLAG != "P"]
	node.dt.m[,LWFLAG:=NULL]
	node.dt.m <- unique(node.dt.m)
	
	#########################################################
	###generate table for tracts that share a single point###
	#########################################################
	pt.dt <- melt(node.dt, id.vars = c("TNID"), measure = list(c("USCB_block_this.x", "USCB_block_this.y"), c("LWFLAG.x", "LWFLAG.y")), value.name = c("USCB_block_this", "LWFLAG"))

	###remove tables that are no longer needed###
	rm(node.dt1,node.dt2,node.dt)

	###remove water tracts###
	pt.dt <- pt.dt[LWFLAG != "P"]
	pt.dt[,c('variable','LWFLAG'):=NULL]
	pt.dt <- unique(pt.dt)

	pt.dt <- merge(pt.dt,pt.dt,by="TNID",allow.cartesian=TRUE)
	pt.dt <- pt.dt[USCB_block_this.x != USCB_block_this.y]

	###remove piers###
	pt.dt[,is.pier := ifelse((substr(USCB_block_this.x,1,5) != substr(USCB_block_this.y,1,5)) & ((USCB_block_this.x %in% piers.dt$USCB_block_this) | (USCB_block_this.y %in% piers.dt$USCB_block_this)), 1, 0)]

	pt.dt <- pt.dt[is.pier==0]

	pt.dt[,USCB_tract_this.1 := ifelse(as.numeric(substr(USCB_block_this.x,1,11)) < as.numeric(substr(USCB_block_this.y,1,11)), substr(USCB_block_this.x,1,11), substr(USCB_block_this.y,1,11))]

	pt.dt[,USCB_tract_this.2 := ifelse(as.numeric(substr(USCB_block_this.x,1,11)) < as.numeric(substr(USCB_block_this.y,1,11)), substr(USCB_block_this.y,1,11), substr(USCB_block_this.x,1,11))]

	pt.dt[,c("TNID","USCB_block_this.x","USCB_block_this.y","is.pier"):=NULL]
	pt.dt <- unique(pt.dt)[USCB_tract_this.1 != USCB_tract_this.2]

	pt.dt <- merge(pt.dt, neighbors.dt, by=c("USCB_tract_this.1", "USCB_tract_this.2"), all.x=TRUE)
	pt.dt <- pt.dt[is.na(edge_len) & !(is.na(USCB_tract_this.1)) & !(is.na(USCB_tract_this.2))]
	pt.dt[,edge_len := 0]

	if(nrow(pt.dt)>0){
		###combine point and line results###
		neighbors.dt <- rbindlist(list(neighbors.dt,pt.dt), use.names=TRUE, fill=TRUE)
	}

	neighbors.dt[,type := ifelse(as.numeric(edge_len)==0,"shared point","shared edge")]

	########################
	###next step: bridges###
	########################

	###if line is a bridge, capture all segments###
	###the rational being that bridges do not start immediatly at the shore but somewhat inland, instead###
	bridges_names <- unique(edges.dt[PASSFLG == "B" & grepl("Bri?d?ge?",FULLNAME,ignore.case=TRUE)]$FULLNAME)

	###return all bridge segments that are not railroad or boundary line types###
	bridges.dt <- edges.dt[((PASSFLG == "B") | (FULLNAME %in% bridges_names)) & !(MTFCC %in% c("P0001","P0004","R1011"))]
	
	#################################################
	###use igraph to form bridges from edges lines###
	#################################################

	e.dt <- unique(bridges.dt[,c("TNIDF","TNIDT"),with=FALSE])

	###generate graph object###
	net <- graph_from_data_frame(d=e.dt, directed=F) 
	 
	#dg <- decompose.graph(net)
	dg <- decompose(net)

	dg.dt <- rbindlist(lapply(1:length(dg),function(i){
		return(unique(data.table(TNID=as.numeric(trimws(V(dg[[i]])$name)), b_grp=i)))
	}),use.names=TRUE,fill=TRUE)
	
	
	#####################################################
	###capture all start and end nodes in bridge lines###
	#####################################################
	bridges.nodes.dt <- data.table(TNID=c(bridges.dt$TNIDF,bridges.dt$TNIDT))[,.(tot=.N),by=TNID]
	bridges.nodes.dt <- merge(bridges.nodes.dt,bridges.dt[,.(from_tot = .N),by=TNIDF],by.x="TNID",by.y="TNIDF",all.x=TRUE)
	bridges.nodes.dt <- merge(bridges.nodes.dt,bridges.dt[,.(to_tot = .N),by=TNIDT], by.x="TNID",by.y="TNIDT",all.x=TRUE)
	bridges.nodes.dt[,from_tot := ifelse(is.na(from_tot),0,from_tot)]
	bridges.nodes.dt[,to_tot := ifelse(is.na(to_tot),0,to_tot)]
	bridges.nodes.dt <- merge(bridges.nodes.dt, dg.dt, by="TNID", all.x=TRUE)

	bridges.nodes.dt_ends <- unique(merge(bridges.nodes.dt[to_tot==0 | from_tot==0], node.dt.m, by="TNID")[,c("b_grp","USCB_tract_this"),with=FALSE])

	bridges.nodes.dt_ends <- merge(bridges.nodes.dt_ends,bridges.nodes.dt_ends,by="b_grp", allow.cartesian=TRUE)
	bridges.nodes.dt_ends <- unique(bridges.nodes.dt_ends)[USCB_tract_this.x != USCB_tract_this.y]

	bridges.nodes.dt_ends[,USCB_tract_this.1 := ifelse(as.numeric(USCB_tract_this.x) < as.numeric(USCB_tract_this.y),USCB_tract_this.x,USCB_tract_this.y)]

	bridges.nodes.dt_ends[,USCB_tract_this.2 := ifelse(as.numeric(USCB_tract_this.x) < as.numeric(USCB_tract_this.y),USCB_tract_this.y,USCB_tract_this.x)]

	bridges.nodes.dt_ends <- unique(bridges.nodes.dt_ends[,c("b_grp","USCB_tract_this.1","USCB_tract_this.2"),with=FALSE])

	###get length by node###
	bridges.dt <- merge(bridges.dt, dg.dt, by.x = "TNIDF", by.y="TNID", all.x=TRUE)

	###create table connecting tracts on each side of bridge###
	bridges.dt_agg <- bridges.dt[,.(edge_len=sum(edge_len)),by=list(b_grp,FULLNAME)]
	bridges.dt_agg[,bridge_length := sum(edge_len), by=b_grp]
	bridges.dt_agg <- bridges.dt_agg[bridges.dt_agg[, .I[which.max(edge_len)], by=list(b_grp)]$V1]

	bridges.dt_agg[,edge_len := NULL]
	bridges.dt_agg <- merge(bridges.dt_agg, bridges.nodes.dt_ends, by="b_grp", all.x=TRUE)
	bridges.dt_agg <- bridges.dt_agg[bridges.dt_agg[, .I[which.min(bridge_length)], by=list(USCB_tract_this.1, USCB_tract_this.2)]$V1]


	###remove bridge relationships that are already in neighbor relationships table###
	keep.cols <- names(bridges.dt_agg)
	bridges.dt_agg <- merge(bridges.dt_agg,neighbors.dt,by=c("USCB_tract_this.1","USCB_tract_this.2"),all.x=TRUE)
	bridges.dt_agg <- bridges.dt_agg[is.na(edge_len),keep.cols,with=FALSE]
	bridges.dt_agg <- bridges.dt_agg[!is.na(USCB_tract_this.1) & !is.na(USCB_tract_this.2)]
	bridges.dt_agg[,b_grp := NULL]
	setorder(bridges.dt_agg,USCB_tract_this.1,USCB_tract_this.2)

	bridges.dt_agg[,type := "bridge connection"]
	
	
	#######################################
	###tracts connected by two addresses###
	#######################################

	if(!is.null(ADDR_dt)) {
	if("data.frame" %in% class(ADDR_dt)) {
	
	if(all(c("ADDR","CITY","STATE","ZIP","group_ID") %in% names(ADDR_dt)) & (nrow(ADDR_dt) > 1)) {
		
		addr.dt <- unique(copy(ADDR_dt)[,c("ADDR","CITY","STATE","ZIP","group_ID"), with=FALSE])
	
		###work around to replace removed censusxy package### 
		gc_URL <- 'https://geocoding.geo.census.gov/geocoder/geographies/onelineaddress?'
	
		gc.dt <- rbindlist(lapply(1:nrow(addr.dt), function(n) {
		
			temp.dt <- addr.dt[n]
		
			this.URL <- paste0(gc_URL,'address=',addr.dt[n]$ADDR,',',addr.dt[n]$CITY,',',addr.dt[n]$STATE,',',addr.dt[n]$ZIP,'&benchmark=2020&vintage=',geo.year,'&format=json')

			req <- jsonlite::fromJSON(URLencode(this.URL))$result$addressMatches
			
			if(length(req) > 0){
				if('geographies' %in% names(req)) {
					if("Census Tracts" %in% names(req$geographies)) {
						if(length(req$geographies[["Census Tracts"]]) > 0) {
							if(nrow(req$geographies[["Census Tracts"]][[1]])>0){
								ct <- req$geographies[["Census Tracts"]][[1]]$GEOID
							}
						}
					}
				}
			} else {
				ct <- NA
			}
			
			temp.dt[,USCB_tract_this := ct]
			
			return(temp.dt)
		
		}),use.names=TRUE)
	
		
		gc.dt <- gc.dt[!is.na(USCB_tract_this)]
	
		if(nrow(addr.dt) > nrow(gc.dt)){
			warning("Some addresses in your address table failed to geocode. Please check your address table.\n")
		}
		
		gc.dt[,tot := .N, by=group_ID]
		
		if(nrow(gc.dt) > nrow(gc.dt[tot > 1])){
			warning("Some addresses in your address table are without fellow group members. Please check your address table.\n")
		}
		
		gc.dt <- gc.dt[tot > 1]

		if(nrow(gc.dt) > 0){
			gc.dt <- merge(gc.dt[,c('group_ID','USCB_tract_this'),with=FALSE],gc.dt[,c('group_ID','USCB_tract_this'),with=FALSE],by="group_ID")[USCB_tract_this.x != USCB_tract_this.y]
			
			gc.dt[,USCB_tract_this.1 := ifelse(as.numeric(USCB_tract_this.x) < as.numeric(USCB_tract_this.y),USCB_tract_this.x,USCB_tract_this.y)]

			gc.dt[,USCB_tract_this.2 := ifelse(as.numeric(USCB_tract_this.x) < as.numeric(USCB_tract_this.y),USCB_tract_this.y,USCB_tract_this.x)]

			gc.dt[,c("USCB_tract_this.x","USCB_tract_this.y"):=NULL]

			gc.dt <- gc.dt[gc.dt[, .I[which.min(group_ID)], by=list(USCB_tract_this.1, USCB_tract_this.2)]$V1]

			###remove address pair relationships that are already in neighbor relationships table###
			keep.cols <- names(gc.dt)
			gc.dt <- merge(gc.dt,neighbors.dt,by=c("USCB_tract_this.1","USCB_tract_this.2"),all.x=TRUE)
			gc.dt <- gc.dt[is.na(edge_len),keep.cols,with=FALSE]
			setorder(gc.dt,USCB_tract_this.1,USCB_tract_this.2)

			if(isTRUE(use.bridges)){
				###remove address pair relationships that are already in bridges relationships table###
				gc.dt <- merge(gc.dt,bridges.dt_agg,by=c("USCB_tract_this.1","USCB_tract_this.2"),all.x=TRUE)
				gc.dt <- gc.dt[is.na(bridge_length),keep.cols,with=FALSE]
				setorder(gc.dt,USCB_tract_this.1,USCB_tract_this.2)
			} 
			
			###retain address connections###
			omit.dt <- omit.dt[!(USCB_tract_this %in% unique(c(gc.dt$USCB_tract_this.1,gc.dt$USCB_tract_this.2)))]
			
			gc.dt[,type := "manual"]
			
		} 
	} else{
			warning("Your address table is missing fields.  Please check your address table.\n")
		}
	}
	}	
	
	###merge to relationship tables###
	
	all_pairs.dt <- neighbors.dt[,c("USCB_tract_this.1","USCB_tract_this.2","type"),with=FALSE]
	
	if(exists("gc.dt")){
	if("data.table" %in% class(gc.dt)){
		if(all(c("group_ID","USCB_tract_this.1","USCB_tract_this.2","type" ) %in% names(gc.dt)) & (nrow(gc.dt) > 0)) {
			all_pairs.dt <- rbindlist(list(gc.dt[,c("USCB_tract_this.1","USCB_tract_this.2","type"),with=FALSE],all_pairs.dt), use.names=TRUE, fill=TRUE) 
		}
	}
	}
	
	
	
	if(isTRUE(use.bridges)){
		all_pairs.dt <- rbindlist(list(bridges.dt_agg[,c("USCB_tract_this.1","USCB_tract_this.2","type"),with=FALSE],all_pairs.dt), use.names=TRUE, fill=TRUE) 
		
		###retain bridge connections###
		omit.dt <- omit.dt[!(USCB_tract_this %in% unique(c(bridges.dt_agg$USCB_tract_this.1,bridges.dt_agg$USCB_tract_this.2)))]
		
	}
	
	
	if(isTRUE(omit.park_openspace)){
		omit.ct <- omit.dt[type=="park_openspace"]$USCB_tract_this
		all_pairs.dt <- all_pairs.dt[!(USCB_tract_this.1 %in% omit.ct) & !(USCB_tract_this.2 %in% omit.ct)]
	}
	
	if(isTRUE(omit.unpopulated)){
		omit.ct <- omit.dt[type=="unpopulated"]$USCB_tract_this
		all_pairs.dt <- all_pairs.dt[!(USCB_tract_this.1 %in% omit.ct) & !(USCB_tract_this.2 %in% omit.ct)]
	}

	

	###########################################
	###add rows for tracts without neighbors###
	###########################################
	
	unq.ct <- unique(substr(faces.dt[LWFLAG != "P"]$USCB_block_this,1,11))

	island.dt <- data.table(USCB_tract_this.1 = unq.ct[!(unq.ct %in% unique(c(all_pairs.dt$USCB_tract_this.1,all_pairs.dt$USCB_tract_this.2)))], USCB_tract_this.2 = NA)
	
	island.dt <- merge(island.dt, unique(omit.dt[,c("USCB_tract_this","type"),with=FALSE]), by.x="USCB_tract_this.1", by.y="USCB_tract_this", all.x=TRUE)
	
	island.dt[,type := trimws(paste("self",ifelse(is.na(type),"",type)))]
	
	all_pairs.dt <- rbindlist(list(all_pairs.dt,island.dt), use.names=TRUE, fill=TRUE)

	all_pairs.dt <- all_pairs.dt[(substr(USCB_tract_this.1,1,5) %in% paste0(FIPS.dt$state,FIPS.dt$county)) & ((substr(USCB_tract_this.2,1,5) %in% paste0(FIPS.dt$state,FIPS.dt$county)) | (is.na(USCB_tract_this.2)))]
	
	invisible(gc())
	
	setwd(old.wd)
	
	###save object as TXT file###
	if(!missing(output.path) & !missing(output.file_name) & !is.null(output.path) & !is.null(output.file_name)){
		if(dir.exists(output.path)){
			fwrite(all_pairs.dt[,c('USCB_tract_2010.1','USCB_tract_2010.2'), with=FALSE], file = file.path(output.path,output.file_name), sep = " ", col.names = FALSE, row.names = FALSE)
			cat(paste('\n Tract network file for',geo.year,'saved here:',file.path(output.path,output.file_name)))
		} else {
			warning("\nOutput file path does not exist. File will not be saved.\n")
		}
	} else {
		warning("\nOutput file path and/or file name missing. File will not be saved.\n")
	}
	
	
	if(as.character(geo.year)=="2010") {
		setnames(all_pairs.dt,c("USCB_tract_this.1","USCB_tract_this.2"),c("USCB_tract_2010.1","USCB_tract_2010.2"))
	} else {
		setnames(all_pairs.dt,c("USCB_tract_this.1","USCB_tract_this.2"),c("USCB_tract_2020.1","USCB_tract_2020.2"))
	}
	
	return(all_pairs.dt)
}
