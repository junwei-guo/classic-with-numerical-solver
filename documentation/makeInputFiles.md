# Creation and formatting of model input files {#makeInputFiles}

1. @ref makeMet
2. @ref makeInit
3. @ref ghgfiles
4. @ref makeOther
5. @ref inputFileForm

---

# Preparation of files for meteorological inputs {#makeMet}

CLASSIC requires the meteorological inputs as described here. The input files are netCDF format. A netcdf dump of the files will look like this:

        netcdf pres_v1.1.5_T63_chunked_1700_2017 {
        dimensions:
            lon = 128 ;
            lat = 64 ;
            time = UNLIMITED ; // (464280 currently)
        variables:
            float lon(lon) ;
                lon:standard_name = "longitude" ;
                lon:long_name = "longitude" ;
                lon:units = "degrees_east" ;
                lon:axis = "X" ;
                lon:_Storage = "contiguous" ;
                lon:_Endianness = "little" ;
            float lat(lat) ;
                lat:standard_name = "latitude" ;
                lat:long_name = "latitude" ;
                lat:units = "degrees_north" ;
                lat:axis = "Y" ;
                lat:_Storage = "contiguous" ;
                lat:_Endianness = "little" ;
            double time(time) ;
                time:standard_name = "time" ;
                time:units = "**day as %Y%m%d.%f**" ;
                time:calendar = "proleptic_gregorian" ;
                time:axis = "T" ;
                time:_Storage = "chunked" ;
                time:_ChunkSizes = 1 ;
                time:_Endianness = "little" ;
            float pres(time, lat, lon) ;
                pres:long_name = "Pressure" ;
                pres:units = "Pa" ;
                pres:_FillValue = 9.96921e+36f ;
                pres:missing_value = 9.96921e+36f ;
                pres:_Storage = "chunked" ;
                pres:_ChunkSizes = 464280, 8, 16 ;
                pres:_Endianness = "little" ;

        }
Importantly,

- Only one variable per file besides lon, lat, and time.
- Note the time units.
- The file is chunked (which depends on the grid being used. What you see here is optimal for T63 global runs). Chunking is not needed for site-level runs. More on chunking in @ref inputFileForm.
- CLASSIC expects the first time step to be 0 hour 0 minute (in hh:mm:ss format 00:00:00).

If you have existing ACSCII met files from a site, see [ASCII to NetCDF met file loader](@ref asciiMet) to use the provided tool to convert them to the appropriate netCDF format.

# Preparation of the model initialization file {#makeInit}

If you have the old format .INI (and CTEM's .CTM) initialization files, there is a [tool to convert them to netCDF format](@ref initTool) for use in CLASSIC. The tool itself is located in tools/initFileConverter.

If you have a netCDF format initialization/restart file that you wish to edit, you can use the [script created for that purpose](@ref modifyRS). It is located in tools/modifyRestartFile. This is designed for regional/global scale edits. For site-level files (single point) it is likely easier to use a combination of ncdump/ncgen.

An ncdump of a properly formatted, global-scale, initialization file is included below. This initialization file is setup for a biogeochemistry (CTEM on) run with peatlands and PFT competition variables included. Note that for a physics-only run (CTEM off), many of the variables will not be read in/required, similarly for a CTEM on run some CLASS-only variables are not required. See the [manual's mainpage](@ref main) for links to sections describing the variables.

          netcdf initFile {
            dimensions:
            	tile = 1 ;
            	lat = 64 ;
            	lon = 128 ;
            	icp1 = 5 ;
            	layer = 20 ;
            	ic = 4 ;
            	icc = 9 ;
            	iccp1 = 10 ;
            	iccp2 = 11 ;
            	months = 12 ;
            	slope = 8 ;
            variables:
            	float ALBS(tile, lat, lon) ;
            		ALBS:_FillValue = -999.f ;
            		ALBS:units = "-" ;
            		ALBS:long_name = "Snow albedo" ;
            	float ALIC(tile, icp1, lat, lon) ;
            		ALIC:_FillValue = -999.f ;
            		ALIC:units = "-" ;
            		ALIC:long_name = "Average near-IR albedo of vegetation category when fully-leafed" ;
            	float ALVC(tile, icp1, lat, lon) ;
            		ALVC:_FillValue = -999.f ;
            		ALVC:units = "-" ;
            		ALVC:long_name = "Average visible albedo of vegetation category when fully-leafed" ;
            	float CLAY(tile, layer, lat, lon) ;
            		CLAY:_FillValue = -999.f ;
            		CLAY:units = "%" ;
            		CLAY:long_name = "Percentage clay content" ;
            		CLAY:provenance = "GSDE soil texture" ;
            	float CMAS(tile, ic, lat, lon) ;
            		CMAS:_FillValue = -999.f ;
            		CMAS:units = "$[kg m^{-2} ]$" ;
            		CMAS:long_name = " Annual maximum canopy mass for vegetation category" ;
            	float Cmossmas(tile, lat, lon) ;
            		Cmossmas:_FillValue = -999.f ;
            		Cmossmas:units = "kgC/m2" ;
            		Cmossmas:long_name = "C in moss biomass" ;
            	float DELZ(layer) ;
            		DELZ:_FillValue = -999.f ;
            		DELZ:units = "m" ;
            		DELZ:long_name = "Ground layer thickness" ;
            	float DRN(tile, lat, lon) ;
            		DRN:_FillValue = -999.f ;
            		DRN:units = "-" ;
            		DRN:long_name = "Soil drainage index" ;
            	float FARE(tile, lat, lon) ;
            		FARE:_FillValue = -999.f ;
            		FARE:units = "fraction" ;
            		FARE:long_name = "Tile fractional area of gridcell" ;
            	float FCAN(tile, icp1, lat, lon) ;
            		FCAN:_FillValue = -999.f ;
            		FCAN:units = "-" ;
            		FCAN:long_name = "Annual maximum fractional coverage of modelled area (read in for CLASS only runs)" ;
            	int GC(lat, lon) ;
            		GC:_FillValue = -999 ;
            		GC:units = "-" ;
            		GC:long_name = "GCM surface descriptor - land surfaces (inc. inland water) is -1" ;
            	float GRO(tile, lat, lon) ;
            		GRO:_FillValue = -999.f ;
            		GRO:units = "-" ;
            		GRO:long_name = "Vegetation growth index" ;
            	float LNZ0(tile, icp1, lat, lon) ;
            		LNZ0:_FillValue = -999.f ;
            		LNZ0:units = "-" ;
            		LNZ0:long_name = "Natural logarithm of maximum vegetation roughness length" ;
            	int MID(tile, lat, lon) ;
            		MID:_FillValue = -999 ;
            		MID:units = "-" ;
            		MID:long_name = "Mosaic tile type identifier (1 for land surface, 0 for inland lake)" ;
            	float ORGM(tile, layer, lat, lon) ;
            		ORGM:_FillValue = -999.f ;
            		ORGM:units = "%" ;
            		ORGM:long_name = "Percentage organic matter content" ;
            		ORGM:provenance = "Zobler soil texture" ;
            	float PAMN(tile, ic, lat, lon) ;
            		PAMN:_FillValue = -999.f ;
            		PAMN:units = "-" ;
            		PAMN:long_name = "Annual minimum plant area index of vegetation category" ;
            	float PAMX(tile, ic, lat, lon) ;
            		PAMX:_FillValue = -999.f ;
            		PAMX:units = "-" ;
            		PAMX:long_name = "Annual maximum plant area index of vegetation category" ;
            	float RCAN(tile, lat, lon) ;
            		RCAN:_FillValue = -999.f ;
            		RCAN:units = "kg/m2" ;
            		RCAN:long_name = "Intercepted liquid water stored on canopy" ;
            	float RHOS(tile, lat, lon) ;
            		RHOS:_FillValue = -999.f ;
            		RHOS:units = "kg/m3" ;
            		RHOS:long_name = "Density of snow" ;
            	float ROOT(tile, ic, lat, lon) ;
            		ROOT:_FillValue = -999.f ;
            		ROOT:units = "m" ;
            		ROOT:long_name = "Annual maximum rooting depth of vegetation category" ;
            	float SAND(tile, layer, lat, lon) ;
            		SAND:_FillValue = -999.f ;
            		SAND:units = "%" ;
            		SAND:long_name = "Percentage sand content" ;
            		SAND:provenance = "Zobler soil texture" ;
            	float SCAN(tile, lat, lon) ;
            		SCAN:_FillValue = -999.f ;
            		SCAN:units = "kg/m2" ;
            		SCAN:long_name = "Intercepted frozen water stored on canopy" ;
            	float SDEP(tile, lat, lon) ;
            		SDEP:_FillValue = -999.f ;
            		SDEP:units = "m" ;
            		SDEP:long_name = "Soil permeable depth" ;
            		SDEP:provenance = "Shangguan et al. (2017) soil depths" ;
            	float SNO(tile, lat, lon) ;
            		SNO:_FillValue = -999.f ;
            		SNO:units = "kg/m2" ;
            		SNO:long_name = "Mass of snow pack" ;
            	int SOCI(tile, lat, lon) ;
            		SOCI:_FillValue = -999 ;
            		SOCI:units = "index" ;
            		SOCI:long_name = "Soil colour index" ;
            	float TBAR(tile, layer, lat, lon) ;
            		TBAR:_FillValue = -999.f ;
            		TBAR:units = "C" ;
            		TBAR:long_name = "Temperature of soil layers" ;
            	float TCAN(tile, lat, lon) ;
            		TCAN:_FillValue = -999.f ;
            		TCAN:units = "C" ;
            		TCAN:long_name = "Vegetation canopy temperature" ;
            	float THIC(tile, layer, lat, lon) ;
            		THIC:_FillValue = -999.f ;
            		THIC:units = "m3/m3" ;
            		THIC:long_name = "Volumetric frozen water content of soil layers" ;
            	float THLQ(tile, layer, lat, lon) ;
            		THLQ:_FillValue = -999.f ;
            		THLQ:units = "m3/m3" ;
            		THLQ:long_name = "Volumetric liquid water content of soil layers" ;
            	float TPND(tile, lat, lon) ;
            		TPND:_FillValue = -999.f ;
            		TPND:units = "C" ;
            		TPND:long_name = "Temperature of ponded water" ;
            	float TSNO(tile, lat, lon) ;
            		TSNO:_FillValue = -999.f ;
            		TSNO:units = "C" ;
            		TSNO:long_name = "Snowpack temperature" ;
            	float ZPND(tile, lat, lon) ;
            		ZPND:_FillValue = -999.f ;
            		ZPND:units = "m" ;
            		ZPND:long_name = "Depth of ponded water on surface" ;
            	float anndefct(lat, lon) ;
            		anndefct:_FillValue = -999.f ;
            		anndefct:units = "mm" ;
            		anndefct:long_name = "Annual water deficit (PFTCompetition variable)" ;
            	float annpcp(lat, lon) ;
            		annpcp:_FillValue = -999.f ;
            		annpcp:units = "mm" ;
            		annpcp:long_name = "Annual precipitation (PFTCompetition variable)" ;
            	float annsrpls(lat, lon) ;
            		annsrpls:_FillValue = -999.f ;
            		annsrpls:units = "mm" ;
            		annsrpls:long_name = "Annual water surplus (PFTCompetition variable)" ;
            	float aridity(lat, lon) ;
            		aridity:_FillValue = -999.f ;
            		aridity:units = "-" ;
            		aridity:long_name = "Aridity index, ratio of potential evaporation to precipitation (PFTCompetition variable)" ;
            	float bleafmas(tile, icc, lat, lon) ;
            		bleafmas:_FillValue = -999.f ;
            		bleafmas:units = "kgC/m2" ;
            		bleafmas:long_name = "Brown leaf mass" ;
            	float defctmon(lat, lon) ;
            		defctmon:_FillValue = -999.f ;
            		defctmon:units = "months" ;
            		defctmon:long_name = "Number of months in a year with water deficit i.e. precipitation less than potential evaporation (PFTCompetition variable)" ;
            	float dmoss(tile, lat, lon) ;
            		dmoss:_FillValue = -999.f ;
            		dmoss:units = "m" ;
            		dmoss:long_name = "Depth of living moss" ;
            	float dry_season_length(lat, lon) ;
            		dry_season_length:_FillValue = -999.f ;
            		dry_season_length:units = "months" ;
            		dry_season_length:long_name = "Length of dry season (PFTCompetition variable)" ;
            	float fcancmx(tile, icc, lat, lon) ;
            		fcancmx:_FillValue = -999.f ;
            		fcancmx:units = "-" ;
            		fcancmx:long_name = "PFT fractional coverage per grid cell" ;
            	float gdd5(lat, lon) ;
            		gdd5:_FillValue = -999.f ;
            		gdd5:long_name = "Growing degree days above 5 C (PFTCompetition variable)" ;
            		gdd5:units = "days" ;
            	float gleafmas(tile, icc, lat, lon) ;
            		gleafmas:_FillValue = -999.f ;
            		gleafmas:units = "kgC/m2" ;
            		gleafmas:long_name = "Green leaf mass" ;
            	float grclarea(lat, lon) ;
            		grclarea:_FillValue = -999.f ;
            		grclarea:units = "km2" ;
            		grclarea:long_name = "Area of grid cell" ;
            	double ic(ic) ;
            		ic:_FillValue = NaN ;
            		ic:long_name = "CLASS PFTs (needleleaved tree, broadleaved tree, crops, grass)" ;
            	double icc(icc) ;
            		icc:_FillValue = NaN ;
            		icc:long_name = "CTEM PFTs (tree:NDL-EVG, NDL-DCD, BDL-EVG, BDL-COLD, BDL-DRY; crops C3,C4; grasses C3,C4)" ;
            	double iccp1(iccp1) ;
            		iccp1:_FillValue = NaN ;
            		iccp1:long_name = "CTEM PFTs + bareground" ;
            	double icp1(icp1) ;
            		icp1:_FillValue = NaN ;
            		icp1:long_name = "CLASS PFTs + bareground" ;
            	float ipeatland(tile, lat, lon) ;
            		ipeatland:_FillValue = -999.f ;
            		ipeatland:units = "-" ;
            		ipeatland:long_name = "Peatland flag: 0 = not a peatland, 1= bog, 2 = fen" ;
            	double lat(lat) ;
            		lat:_FillValue = NaN ;
            		lat:actual_range = "-87.8637987364, -87.8637987364" ;
            		lat:units = "degrees_north" ;
            		lat:long_name = "latitude" ;
            		lat:standard_name = "latitude" ;
            		lat:axis = "Y" ;
            	double layer(layer) ;
            		layer:_FillValue = NaN ;
            		layer:long_name = "ground column layers" ;
            	float lfstatus(tile, icc, lat, lon) ;
            		lfstatus:_FillValue = -999.f ;
            		lfstatus:units = "-" ;
            		lfstatus:long_name = "Leaf status, see Phenology" ;
            	double litrmass(tile, iccp2, lat, lon) ;
            		litrmass:_FillValue = NaN ;
            		litrmass:units = "kgC/m2" ;
            		litrmass:long_name = "Litter mass per soil layer" ;
            	float litrmsmoss(tile, lat, lon) ;
            		litrmsmoss:_FillValue = -999.f ;
            		litrmsmoss:units = "kgC/m2" ;
            		litrmsmoss:long_name = "Moss litter mass" ;
            	double lon(lon) ;
            		lon:_FillValue = NaN ;
            		lon:actual_range = "0.0, 0.0" ;
            		lon:units = "degrees_east" ;
            		lon:long_name = "longitude" ;
            		lon:standard_name = "longitude" ;
            		lon:axis = "X" ;
            	double months(months) ;
            		months:_FillValue = NaN ;
            		months:long_name = "Months" ;
            	int nmtest(lat, lon) ;
            		nmtest:_FillValue = -999 ;
            		nmtest:long_name = "Number of tiles in each grid cell" ;
            	float pandays(tile, icc, lat, lon) ;
            		pandays:_FillValue = -999.f ;
            		pandays:units = "-" ;
            		pandays:long_name = "Days with +ve new photosynthesis, see Phenology" ;
            	float rootmass(tile, icc, lat, lon) ;
            		rootmass:_FillValue = -999.f ;
            		rootmass:units = "kgC/m2" ;
            		rootmass:long_name = "Root mass" ;
                double grwtheff(tile, icctem, lat, lon) ;
                	grwtheff:units = "(kg C/m^2)/(m2/m2)" ;
                	grwtheff:_FillValue = -999. ;
                	grwtheff:long_name = "Growth efficiency. Change in biomass per year per unit max. LAI,for use in mortality subroutine" ;
            	double slope(slope) ;
            		slope:_FillValue = NaN ;
            		slope:long_name = "wetland slope fractions for 0.025, 0.05, 0.1, 0.15, 0.20, 0.25, 0.3 and 0.35 percent slope threshold" ;
            	double soilcmas(tile, iccp2, lat, lon) ;
            		soilcmas:_FillValue = NaN ;
            		soilcmas:units = "kgC/m2" ;
            		soilcmas:long_name = "Soil C mass per soil layer" ;
            	float srplsmon(lat, lon) ;
            		srplsmon:_FillValue = -999.f ;
            		srplsmon:units = "months" ;
            		srplsmon:long_name = "Number of months in a year with surplus water i.e. precipitation more than potential evaporation (PFTCompetition variable)" ;
            	float stemmass(tile, icc, lat, lon) ;
            		stemmass:_FillValue = -999.f ;
            		stemmass:units = "kgC/m2" ;
            		stemmass:long_name = "Stem mass" ;
            	float tcoldm(lat, lon) ;
            		tcoldm:_FillValue = -999.f ;
            		tcoldm:units = "C" ;
            		tcoldm:long_name = "Temperature of the coldest month (PFTCompetition variable)" ;
            	double tile(tile) ;
            		tile:_FillValue = NaN ;
            		tile:long_name = "tiles" ;
            		tile:axis = "Z" ;
            	float twarmm(lat, lon) ;
            		twarmm:_FillValue = -999.f ;
            		twarmm:units = "C" ;
            		twarmm:long_name = "Temperature of the warmest month (PFTCompetition variable)" ;
            	float slopefrac(tile, slope, lat, lon) ;
            		slopefrac:_FillValue = -999.f ;
            		slopefrac:units = "-" ;
            		slopefrac:long_name = "Slope-based fraction for dynamic wetlands" ;
            	double iccp2(iccp2) ;
            		iccp2:_FillValue = NaN ;
            		iccp2:long_name = "CTEM PFTs, bareground, LUC product pools" ;
            	double maxAnnualActLyr(tile, lat, lon) ;
            		maxAnnualActLyr:long_name = "!< Active layer thickness maximum over the e-folding period specified by parameter eftime" ;
            		maxAnnualActLyr:units = "m" ;
            		maxAnnualActLyr:_FillValue = -999. ;
            		maxAnnualActLyr:missing_value = -999. ;
          }

# Greenhouse gas inputs files {#ghgfiles}

If you have existing ACSCII GHG files, see [our tool to create GHG input files](@ref makeGHGfiles) to convert them to the appropriate netCDF format.

# Making input files for other input variables {#makeOther}

Most other CLASSIC input files are not desired or needed for point runs of the model. If you require input files for regional or global simulations we may be able to provide you with versions we use in our runs. The possible inputs are listed in Additional inputs depending on model configuration. Because these files are not required of most users we have not set up tools to help generate the files.

# Some notes on input file format {#inputFileForm}

These points below are most relevant for those running globally or at high resolution regionally. They will have minimal impact upon site-level runs. 

**Input File Format**

Because CLASSIC utilizes the parallel version of netCDF-4 libraries, all input files must be in netCDF-4 format. Many data files continue to be distributed in netCDF "classic" format, sometimes referred to as netCDF-3. There are various tools available to convert from one format to the other. One such utilitiy is "ncks" which is part of the [nco package](http://nco.sourceforge.net).  For example:
```
ncks -4 input_file output_file
```
Further information regarding netCDF can be found on Unidata's [FAQ](https://www.unidata.ucar.edu/software/netcdf/docs/faq.html). 

> **Simply converting to netCDF-4 format is not usually all that is required to obtain good I/O performance when running the parallel version of CLASSIC. Please note carefully the following section.**

**Chunking Input Files**

One of the requirements of the netCDF-4 format is that if a record dimension is present, the file must be chunked. So when a netCDF classic file with an unlimited time dimension is converted to netCDF-4 format, it is automatically "chunked". Chunking is meant to enhance I/O performance, but this is only the case if the chunks are chosen in accordance to how the file is to be accessed. See this link for details of why [chunking matters](https://www.unidata.ucar.edu/blogs/developer/en/entry/chunking_data_why_it_matters). `It is important to note that a poor choice of chunk sizes can, in fact, lead to very poor I/O performance. Popular utilities such as "ncks" and "cdo" typically do not, by default, chunk netCDF-4 files in a way that works well with CLASSIC.`

The parallel version of CLASSIC runs the entire simulation at each grid cell in a domain on an individual processor. This means that all time steps in the input file must be read in for that grid cell at the very start of the simulation. In order to ensure good read performance, the input files should be chunked so that all the time steps for a given grid cell are included in the same chunk. However, it is not recommended to chunk a file so that each chunk only contains the time steps at a single grid point as it would mean very poor read performance when accessing the file as 2D grids in latitude/longitude. The best compromise is to chunk the files so that each chunk contains all time steps for multiple grid points. For example, using ncks on a netCDF file of T63 global grids with a "time" dimension consisting of 39420 time steps, a "lat" dimension of 64 latitudes and a "lon" dimension of 128 longitudes, a suitable command line would be: 
```
ncks -h -4 --cnk_plc=g3d --cnk_dmn=time,39420 --cnk_dmn=lat,8 --cnk_dmn=lon,16 input_file output_file
```
And for the high resolution Canada domain (310x160 rotated lat-lon grid):
```
ncks -h -4 --cnk_plc=g3d --cnk_dmn=time,39420 --cnk_dmn=lat,10 --cnk_dmn=lon,10 input_file output_file
```
It is recommended that:

1.  All processing of the input file with other utilities such as "cdo" or other "nco" tools be completed before this final step since those utilities may rechunk the file using different chunk sizes.
2.  The input file is in netCDF classic format so that the chunking and conversion to netCDF-4 format (via the "-4" switch) is done at the same time.
3.  The entire number of time steps in the input file is used for the chunk size in the time dimension. 
4.  The chunk sizes for the lat and lon dimensions in the above examples be used for those domains. However, it may be that for files with longer time series, the chunk sizes for the lat and lon dimensions need to be adjusted (i.e. reduced).
5.  The output file is checked to ensure that the chunking was done properly. This can be done using "ncdump -hs output_file". The output from "ncdump" should show the variable with the chunk sizes used by the "ncks" utility and the time variable should also be chunked with the total number of time steps. The lat and lon variables will not be chunked (i.e. will be shown as "contiguous" in the ncdump output).

> **Caveats for ECCC users**

> **Local implementations of nco/ncks may not support all the features needed to properly process and chunk files suitable for use on the ECCC Science network. On each ppp, /usr/bin/ncks appears to work well and it is recommended to use it rather than a local implementation of nco to do the chunking if there are issues encountered locally.**

> **Chunking (via ncks) requires that the entire file be placed into memory. If there is insufficient RAM, the chunking will fail. Note that there is over 200 GB of RAM available on each node of the ppp.**
