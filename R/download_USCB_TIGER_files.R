###disable scientific notation###
options(scipen = 999)
options(timeout=240)

###load packages###
library(data.table)

data.table::setDTthreads(1)

d_f <- function(this.URL,main.path,USCB_TIGER.path) {
	file <- basename(this.URL)
				
	###download commpressed file###
	download.file(this.URL, file)
				
	###unzip compressed file###
	unzip(file, exdir = file.path(main.path,tools::file_path_sans_ext(file)))
				
	###removed compressed file###
	file.remove(file.path(USCB_TIGER.path,file))
}

download_USCB_TIGER_files <- function(FIPS_dt,USCB_TIGER.path,geo.year="2010"){

	if(as.character(geo.year)=="2010") {
		f_year <- "2019"
	} else {
		f_year <- "2022"
	}
	
	FIPS.dt <- copy(as.data.table(FIPS_dt))

	###pad state and county codes with leading zeros###
	FIPS.dt[,state := sprintf("%02d", as.numeric(state))]
	FIPS.dt[,county := sprintf("%03d", as.numeric(county))]
	
	FIPS.dt <- unique(FIPS.dt[,c("state","county"),with=FALSE])

	base.URL <- paste0("https://www2.census.gov/geo/tiger/TIGER",f_year)

	old.wd <- getwd()

	setwd(USCB_TIGER.path)
	
	file.dt <- data.table(f_name=c("FACES","EDGES","FACESAL","AREALM"),f_desc=c("FACES","EDGES","Topological Faces-Area Landmark relationship","Area Landmark relationship"),f_type=c("county","county","state","state"))
	
	
	for (k in 1:nrow(file.dt)){
	
		this.f_name <- file.dt[k]$f_name
		this.f_desc <- file.dt[k]$f_desc
		this.f_type <- file.dt[k]$f_type
	
		main.URL <- file.path(base.URL,this.f_name)
		main.path <- file.path(USCB_TIGER.path,this.f_name)
		
		if (this.f_type=="county") {
			for (j in 1:nrow(FIPS.dt)){
			
				this.URL <- file.path(main.URL,paste0("tl_",f_year,"_",FIPS.dt[j]$state,FIPS.dt[j]$county,"_",tolower(this.f_name),".zip"))

				d_f(this.URL,main.path,USCB_TIGER.path)
			}
		} else{
			for (j in unique(FIPS.dt$state)){
				
				this.URL <- file.path(main.URL,paste0("tl_",f_year,"_",j,"_",tolower(this.f_name),".zip"))

				d_f(this.URL,main.path,USCB_TIGER.path)
			}
		}
		
		cat(paste0("\n\nUSCB TIGER ",this.f_desc," files downloaded.\n\n"))
		
		invisible(gc())
	
	}
	
}	
	
