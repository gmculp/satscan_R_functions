# satscan_R_functions
This repository contains functions for generating SaTScan related items.  It is a work in progress (in other words, there are likely some bugs) and will be updated reglarly. 

At this point in time, the following functions are available:

* ```download_USCB_TIGER_files.R```: Function to download necessary files from USCB TIGER website. *Note that if you are using this code from within an organization's firewall, the USCB TIGER website (https://www2.census.gov/geo/tiger) may need to be whitelisted.*

* ```generate_USCB_tract_network_file.R```: Function to generate census tract relationship file in CSV format. Includes optional logical arguments to enable bridge connectivity, omit connectivity with parks and other open spaces and omit connectivity with unpopulated tracts. In addition, there is the option to input a data frame of to and from addresses.

Required packages that must be installed to run this code:

* ```data.table```: for handling large data.frames more efficiently

* ```foreign```: for reading in dbf files

* ```sf```: for reading in shapefiles

* ```censusapi```: for reading in population data

* ```censusxy```: for geocoding addresses to census block level

* ```igraph```: for collapsing directional multipart polylines in edges files into single part polylines

          
Here is a code sample for generating a network file for tracts within NYC...
```
source("R/download_USCB_TIGER_files.R")
source("R/generate_USCB_tract_network_file.R")

###specify place to store USCB TIGER files###
USCB_TIGER.path <- "C:/SaTScan_resources/census_files"

###specify data table containing state and county FIPS codes###
FIPS_dt <- data.table(state=rep("36",5),county=c("061","005","047","081","085"))

###option to bring in CSV of addresses with column for connection ID###
###each group must consist of at least two addresses###
#Staten Island Ferry connection#
ADDR_dt <- data.table(ADDR=c("1 Bay St","4 South Street"), CITY = c("Staten Island","New York"), STATE = c("NY","NY"), ZIP=c("10301","10004"),group_ID=c(1,1))

###automatically download all necessary files from USCB TIGER website###
###you will only have to do this once###
download_USCB_TIGER_files(FIPS.dt,USCB_TIGER.path)

###generate census tract relationship file###
all_pairs.dt <- generate_USCB_tract_network_file(FIPS_dt, USCB_TIGER.path, omit.park_openspace=TRUE, omit.unpopulated=TRUE, use.bridges=TRUE, ADDR_dt)



