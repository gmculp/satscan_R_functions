###disable scientific notation###
options(scipen = 999)

###load packages###
library(data.table)

data.table::setDTthreads(1)


download_USCB_TIGER_files <- function(FIPS.dt,USCB_TIGER.path){

	###pad state and county codes with leading zeros###
	FIPS.dt[,state := sprintf("%02d", as.numeric(state))]
	FIPS.dt[,county := sprintf("%03d", as.numeric(county))]

	base.URL <- "https://www2.census.gov/geo/tiger/TIGER2019"

	old.wd <- getwd()

	setwd(USCB_TIGER.path)
	
	###############################
	###download faces shapefiles###
	###############################

	main.URL <- file.path(base.URL,"FACES")
	main.path <- file.path(USCB_TIGER.path,"FACES")

	for (j in 1:nrow(FIPS.dt)){

		#j <- 1
		
		this.URL <- file.path(main.URL,paste0("tl_2019_",FIPS.dt[j]$state,FIPS.dt[j]$county,"_faces.zip"))

		file <- basename(this.URL)
		
		###download commpressed file###
		download.file(this.URL, file)
		
		###unzip compressed file###
		unzip(file, exdir = file.path(main.path,tools::file_path_sans_ext(file)))
		
		###removed compressed file###
		file.remove(file.path(USCB_TIGER.path,file))
		
	}
	
	cat("USCB TIGER FACES files downloaded.\n")

	invisible(gc())


	###############################
	###download edges shapefiles###
	###############################

	main.URL <- file.path(base.URL,"EDGES")
	main.path <- file.path(USCB_TIGER.path,"EDGES")

	for (j in 1:nrow(FIPS.dt)){

		#j <- 1
		
		this.URL <- file.path(main.URL,paste0("tl_2019_",FIPS.dt[j]$state,FIPS.dt[j]$county,"_edges.zip"))

		file <- basename(this.URL)
		
		###download commpressed file###
		download.file(this.URL, file)
		
		###unzip compressed file###
		unzip(file, exdir = file.path(main.path,tools::file_path_sans_ext(file)))
		
		###removed compressed file###
		file.remove(file.path(USCB_TIGER.path,file))
		
	}
	
	cat("USCB TIGER EDGES files downloaded.\n")

	invisible(gc())


	################################################################
	###download Topological Faces-Area Landmark Relationship File###
	################################################################

	main.URL <- file.path(base.URL,"FACESAL")
	main.path <- file.path(USCB_TIGER.path,"FACESAL")

	for (j in unique(FIPS.dt$state)){

		#j <- 1
		
		this.URL <- file.path(main.URL,paste0("tl_2019_",j,"_facesal.zip"))

		file <- basename(this.URL)
		
		###download commpressed file###
		download.file(this.URL, file)
		
		###unzip compressed file###
		unzip(file, exdir = file.path(main.path,tools::file_path_sans_ext(file)))
		
		###removed compressed file###
		file.remove(file.path(USCB_TIGER.path,file))
		
	}
	
	cat("USCB TIGER Topological Faces-Area Landmark relationship files downloaded.\n")

	invisible(gc())

	##############################################
	###download Area Landmark Relationship File###
	##############################################

	main.URL <- file.path(base.URL,"AREALM")
	main.path <- file.path(USCB_TIGER.path,"AREALM")

	for (j in unique(FIPS.dt$state)){

		#j <- 1
		
		this.URL <- file.path(main.URL,paste0("tl_2019_",j,"_arealm.zip"))

		file <- basename(this.URL)
		
		###download commpressed file###
		download.file(this.URL, file)
		
		###unzip compressed file###
		unzip(file, exdir = file.path(main.path,tools::file_path_sans_ext(file)))
		
		###removed compressed file###
		file.remove(file.path(USCB_TIGER.path,file))
		
	}
	
	cat("USCB TIGER Area Landmark relationship files downloaded.\n")

	invisible(gc())
	
	setwd(old.wd)

	cat("USCB TIGER file download complete.\n")
}



