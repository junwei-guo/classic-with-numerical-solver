# Basic model inputs to run CLASSIC {#basicInputs}

1. @subpage modelParams
2. @subpage forcingData
3. @subpage vegetationData
  - @subpage vegCLASSonly
  - @subpage vegCTEMtoo
4. @subpage soilData
5. @subpage initProgVar
  - @subpage initPhysProgVar
  - @subpage initCTEMProgVar
6. @subpage exModSets 

----

# Model parameters {#modelParams}

Model parameters are read in from an external fortran namelist file. Three options are supplied in the configurationFiles folder (PFTs for each are listed in @ref PFTsCLASSIC)

- *template_run_parameters.txt* is the default model setup with 9 PFTs as published in Melton and Arora (2016) \cite Melton2016-zx 
- *template_run_parameters_peatlands.txt* additionally includes the peatland PFTs as described in Wu et al. (2016) \cite Wu2016-zt. 
- *template_run_parameters_shrubs.txt* includes shrubs at both the CLASS PFT level and CTEM. (**BETA version**)

In all cases it is important to match up the PFTs in the parameter namelist file with those in the initialization file. E.g. if you want to simulate shrubs you need to ensure you have shrub parameters in your parameter namelist as well as some non-zero fractional cover for shrub PFTs in your initialization file (or LUC file).

The parameter namelist file has two distinct sections. The upper section *&classicPFTbasic* contains the information needed to setup the arrays for variables that vary by PFT and must correspond to the information in the *&classicparams* section of the namelist file (number of PFTs, etc.) 

          &classicPFTbasic

          ! In this header, the basic numbers of PFTs is listed. Below this initial namelist is the parameters for each PFT.

          ! Number of PFTs considered by the physics subroutines. NOTE: The number specified here must match the data in your init netcdf file.
          ican = 4 ,    

          ! Number of CTEM level PFTs. NOTE: The number specified here must match the data in your init netcdf file.
          icc = 9 ,           

          ! Maximum number of level 2 CTEM PFTs. This is the maximum number of CTEM PFTs associated with a single CLASS PFT.
          l2max = 3 ,         

          /

Based upon this information, the variables that have CLASS or CTEM number of PFTs (*ican* and *icc*, respectively) are allocated in allocateParamsCLASSIC in classic_params.

The rest of the parameters are then read in and stored in the classic_params module.

# Atmospheric Forcing Data {#forcingData}

At each physics time step, for each grid cell or modelled area, the following atmospheric forcing data are,

**Required**:
- *FDL* Downwelling longwave sky radiation [\f$ W m^{-2}\f$ ]
- *PRE* Surface precipitation rate [\f$kg m^{-2} s^{-1}\f$ ]
  - CLASSIC is able to run with total incoming precipitation by partitioning it into rainfall and snowfall on the basis of empirically-derived equations. If the rainfall rate (*RPRE*) and snowfall rate (*SPRE*) are available, they could be used instead. The @ref model_state_driver.getMet and @ref model_state_drivers.updateMet subroutines should be modified accordingly, and the job options switch *IPCP* should be set to 4 (more on this in @ref setupJobOpts).
- *PRES* Surface air pressure [\f$P_a\f$ ]
- *QA* Specific humidity at reference height [\f$kg kg^{-1}\f$ ]
- *TA* Air temperature at reference height [degree C]
  - For atmospheric models, the air temperature supplied to CLASS should be the lowest level air temperature extrapolated using the dry adiabatic lapse rate to the bottom of the atmosphere, i.e. to where the wind speed is zero and the pressure is equal to the surface pressure, Pa.
  - For field data, the actual measured air temperature at the reference height should be used, since in this case the adiabatic extrapolation is performed within CLASS.
  - **Note: In CLASSIC, *TA* is in Kelvin, however the driver expects the air temperature in units of degrees Celsius.**
- *VMOD* Wind speed at reference height [\f$m s^-1\f$ ]
  - Atmospheric models provide the zonal and meridional components of the wind velocity, but CLASS does not actually require information on wind direction. Thus, if only the scalar wind speed is available, either *UL* or *VL* can be set to it, and the other to zero. (Both of these terms, plus the scalar wind speed *VMOD*, must be supplied to CLASS.)
- *FSS* Downwelling shortwave radiation incident on a horizontal surface [\f$W m^{-2}\f$ ]
  -CLASS ideally recieves the forcing incoming shortwave radiation partitioned into the visible and near-infrared components. Since these are rarely available, they can each be roughly estimated as approximately half of the total incoming solar radiation (as is done in the model code in main.f90).  Note: if the ISNOALB switch (see below) is set to 1, the incoming shortwave radiation is required in four wavelength bands, and the direct and diffuse components are required as well. At present this model option is not operational (the test version does not use four band inputs)


**Optional**:
- *FCLO* Fractional cloud cover [ ]
  - The fractional cloud cover is used in the calculation of the transmissivity of the vegetation canopy. If it is not available it is estimated on the basis of the solar zenith angle and the occurrence of precipitation. If you do have the fractional cloud cover, you will need to properly input it to the model and then override the estimation done by @ref generalUtils.findCloudiness
- *FSIH* Near infrared shortwave radiation incident on a horizontal surface [\f$W m^{-2}\f$ ] (see above for *FSS*, set in @ref main.f90)
- *FSVH* Visible shortwave radiation incident on a horizontal surface [\f$W m^{-2}\f$ ] (see above for *FSS*, set in @ref main.f90)
- *UL* Zonal component of wind velocity [\f$m s^{-1}\f$ ] (see above for *VMOD*, set in @ref main.f90)
- *VL* Meridional component of wind velocity [\f$m s^{-1}\f$ ] (see above for *VMOD*, set in @ref main.f90)


**Specified**:
- *ZBLD* Atmospheric blending height for surface roughness length averaging [m]
  - If the surface being modelled is a very heterogeneous one, care must be taken to ensure that the reference heights are greater than the “blending height”, the distance above the surface at which the atmospheric variables are not dominated by any one surface type. In principle this height depends on the length scale of the roughness elements; it may be as large as 50-100 m. In CLASSIC the blending height is used in averaging the roughness lengths over the modelled area, and is read in separately from *ZRFM* and *ZRFH* as *ZBLD*.
- *ZRFH* Reference height associated with forcing air temperature and humidity [m]
- *ZRFM* Reference height associated with forcing wind speed [m]
  - In atmospheric models the forcing wind speed, air temperature and specific humidity are obtained from the lowest modelled atmospheric layer, and thus the reference height will be the height above the “surface” (i.e. the location where the wind speed is zero and the pressure is equal to the surface pressure, Pa) corresponding to that lowest layer. Some atmospheric models use a vertical co-ordinate system in which the momentum and thermodynamic levels are staggered, and if so, *ZFRM* and *ZRFH* will have different values. If that is the case, the switch *ISLFD* in the job options file should be set to 2, so that the subroutines @ref FLXSURFZ.f and @ref DIASURFZ.f are called, since the other options do not support different reference heights.
  - In the case of field data, the reference height is the height above the ground surface at which the variables are measured. If the measurement height for wind speed is different from that for the air temperature and specific humidity, again the *ISLFD* switch in the job options file should be set to 2.
  - **Note** neither *ZRFH* nor *ZRFM* may be smaller than the vegetation canopy height, as this will cause the model run to crash.
- *GC* GCM surface descriptor
  - For land surfaces (inc. inland water) set it to -1

  **Note** If you are using gridded meteorology and it is all in 'local' time you will need to set allLocalTime = .true. in metModule.f90.  This will prevent it from adjusting the timezone of your meteorology. For site-level simulations as long as you provide meteorology on the same timestep as the model physics (typically 30 minutes) your meteorology will be used as is. To determine if your meteorology is relative to Greenwich or in local time, look at your shortwave radiation. As you move through time do you see the sun move across longitudes (allLocalTime = .false.) or across latitudes (allLocalTime = .true.). It seems reanalysis will generally be false while climate model outputs are generally true.
  
  **Note** You will need to start on timestep 0 minutes, 0 hours, day 1 for your meteorology. In main.f90 it looks for this timestep to assign various needed variables. You will get weird behaviour without this step!
  
## Advisement regarding the physics timestep

The length of the time step should be carefully considered in assembling the forcing data. CLASS has been designed to run at a time step of 30 minutes or less, and the explicit prognostic time stepping scheme used for the soil, snow and vegetation variables is based on this assumption. Longer time steps may lead to the emergence of numerical instabilities in the modelled prognostic variables. The physics timestep can be changed in the run parameters namelist file.


# Input Vegetation Data {#vegetationData}

CLASSIC can be run with either dynamic vegetation (CTEM+CLASS) or a physics only simulation (CLASS). The model inputs differ between the two configurations with some input vegetation data for a CLASS-only run ignored when CTEM is turned on. As well CTEM requires some additional inputs not needed for CLASS-only runs as described below.

## Required vegetation data {#vegCLASSonly}

The typical main vegetation categories for the model physics (CLASS) are needleleaf trees, broadleaf trees, broadleaf shrubs, crops and grass (e.g. ican = 5). However this is adaptable and can include further PFTs depending on model configuration. There are checks in the model code that look for known PFTs (e.g. 'NdlTr', 'BdlTr', 'BdlSh', 'Crops', 'Grass'). If a PFT is introduced that is not known, the model will bail and let the user know where the PFT check failed. This is designed to prevent poorly considered additions and ensure developers are aware of where the different PFTs branch in the code. Urban areas are also treated as “vegetation” in the CLASS code, and have associated values for *FCAN*, *ALVC*, *ALIC* and *LNZ0* (see below). Thus these arrays have a larger dimension of ican + 1 rather than ican.  

For each of the CLASS PFTs the following data are required for each mosaic tile over each grid cell or modelled area (**NOTE**: When CTEM is turned on, i.e. dynamic vegetation is desired, the variables indicated in bold font are overwritten by CTEM during the model run (specifically in @ref calcLandSurfParams.f90). *FCAN* may be overwritten if land use change or competition between PFTs is turned on.)

 1. *ALIC* Average near-IR albedo of vegetation category when fully-leafed [ ]
 2. *ALVC* Average visible albedo of vegetation category when fully-leafed [ ]
 3. **CMAS** Annual maximum canopy mass for vegetation category [\f$kg m^{-2}\f$ ]
 4. **FCAN** Annual maximum fractional coverage of modelled area [ ]
   - This variable is only read-in when CTEM is off.
 5. **LNZ0** Natural logarithm of maximum vegetation roughness length [ ]
 6. **PAMN** Annual minimum plant area index of vegetation category [ ]
 7. **PAMX** Annual maximum plant area index of vegetation category [ ]
 8. **ROOT** Annual maximum rooting depth of vegetation category [m]


### Variables not presently read in

These variables are not presently read-in by CLASSIC but could be if desired.

9. *PSGA* Soil moisture suction coefficient (used in stomatal resistance calculation) [ ]
10. *PSGB* Soil moisture suction coefficient (used in stomatal resistance calculation) [ ]
11. *QA50* Reference value of incoming shortwave radiation (used in stomatal resistance calculation) [\f$W m^{-2}\f$ ]
12. *RSMN* Minimum stomatal resistance of vegetation category [\f$s m^{-1}\f$ ]
13. *VPDA* Vapour pressure deficit coefficient (used in stomatal resistance calculation) [ ]
14. *VPDB* Vapour pressure deficit coefficient (used in stomatal resistance calculation) [ ]


In **physics only runs (CLASS only)**, the vegetation is prescribed as follows (For full details of these calculations, see the documentation for subroutine @ref calcLandSurfParams.f90):

- CLASS models the physiological characteristics of trees as remaining constant throughout the year except for the leaf area index and plant area index, which vary seasonally between the limits defined by *PAMX* and *PAMN*.
- The areal coverage of crops varies from zero in the winter to *FCAN* at the height of the growing season, and their physiological characteristics undergo a corresponding cycle.
- Grasses remain constant year-round.

Ideally the vegetation parameters should be measured at the modelled location. Of course this is not always possible, especially when running over a large modelling domain. As a guide, the table below provides generic values from the literature for the 20 categories of globally significant vegetation types. If more than one type of vegetation in a given category is present on the modelled area, the parameters for the category should be areally averaged over the vegetation types present.

\image html "landcovercat_table.png" ""
\image latex "landcovercat_table.png" ""

For the non-required stomatal resistance parameters, typical values for the four principal vegetation types are given below:

\f[
\begin{array}{ | l | c | c | c | c | c | c | }
 & \text{RSMN} & \text{QA50} & \text{VPDA} & \text{VPDB} & \text{PSGA} & \text{PSGB} \\
\text{Needleleaf trees} & 200.0 & 30.0 & 0.65 & 1.05 & 100.0 & 5.0 \\
\text{Broadleaf trees} & 125.0 & 40.0 & 0.50 & 0.60 & 100.0 & 5.0 \\
\text{Crops} & 85.0 & 30.0 & 0.50 & 1.00 & 100.0 & 5.0 \\
\text{Grass} & 100.0 & 30.0 & 0.50 & 1.00 & 100.0 & 5.0 \\
\end{array}
\f]

## Required vegetation data for a biogeochemical simulation (CLASS+CTEM) {#vegCTEMtoo}

In addition to the CLASS variables described above, CTEM requires the following further information about the vegetation:

- *fcancmx* Fractional coverage of CTEM PFTs per grid cell [ ]

If land use is being simulated this value will come from a land use change file (see @ref inputLUC) otherwise it is taken from the model initialization file and kept constant thoughout a run (provided competition between PFTs is not turned on).

# Input Soil Data {#soilData}

The following specifications are required for each modelled soil layer:

- *DELZ* Layer thickness [m]

CLASSIC supports any number and thickness of ground layers. However, because the temperature stepping scheme used in CLASS is of an explicit formulation, care must be taken not to make the layers too thin, since this may lead to numerical instability problems. As a rule of thumb, the thicknesses of layers should be limited to \f$\geq\f$ 0.10 m. *DELZ* is the same across all grid cells (the soil permeable depth, *SDEP*, can vary, see below).

For each of the modelled soil layers on each of the mosaic tiles, the following texture data are required:

- *CLAY* Percentage clay content
- *ORGM* Percentage organic matter content
- *SAND* Percentage sand content

\image html "percentSand.png" "Percent Sand"
\image latex "percentSand.png" "Percent Sand"

- For mineral soils, the percentages of sand, clay and organic matter content need not add up to 100%, since the residual is assigned to silt content. If the exact sand, clay and organic matter contents are not known, estimates can be made for the general soil type on the basis of the standard USDA texture triangle shown above. Organic matter contents in mineral soils are typically not more than a few percent.
- If the layer consists of bedrock, *SAND* is assigned a flag value of -3. If it is part of a continental ice sheet, it is assigned a flag value of -4. In both cases, *CLAYROT* and *ORGMROT* are not used and can be set to zero.
- Highly organic soils have different behaviour if the area is being modelled as a peatland:
  - If the peatland flag is 0 (ipeatland = 0 in the model initialization file), that is the tile is not being modelled as a peatland, and the soil layer is a fully organic one (generally takes as \f$\geq\f$ 30% organic matter), *SANDROT*, *CLAYROT* and *ORGMROT* are used differently. The sand content should be assigned a flag value of -2. The peat texture is then assigned to be fibric, hemic or sapric (see Letts et al. (2000) \cite Letts2000-pg). The current default is for the first layer to be assumed as fibric, the second and third as hemic and any lower layers as sapric until bedrock or a sand value > 0 (i.e. mineral soil) is reached. *CLAYROT* is not used and can be set to zero.
  - If the tile is being treated as a peatland then the first soil layer is considered moss following Wu et al. (2016) \cite Wu2016-zt. The lower soil layers are treated such that the lower layers are assigned fibric, hemic or sapric characteristics (see run parameters namelist file)

  \f[
  \begin{array}{ | l | c | c | c |  }
   \text{Soil layer type}  & \text{SAND} & \text{CLAY} & \text{ORGM}  \\
   \text{Mineral}   & >= 0 & >= 0           & >= 0 \\
   \text{Peat}      & -2   & \text{ignored} & \text{ignored} \\
   \text{Bedrock}   & -3   & \text{ignored} & \text{ignored}  \\
   \text{Ice sheet} & -4   & \text{ignored} & \text{ignored}  \\
  \end{array}
  \f]


*SAND*, *CLAY* and *ORGM* are utilized in the calculation of the soil layer thermal and hydraulic properties in @ref soilProperties.f90. If measured values of these properties are available, they could be used instead (with modificiations to the code).

For each of the mosaic tiles over the modelled area, the following surface parameters must be specified:

- *DRN* Soil drainage index
  - The drainage index is usually set to 0.005 except in cases of deep soils where it is desired to suppress drainage from the bottom of the soil profile (e.g. in bogs, or in deep soils with a high water table). In such cases it is set to 0. A value of 1 is freely draining.
- *FARE* Fractional coverage of mosaic tile on the modelled area (also discussed [here](@ref compvsmosaic))
- *MID* Mosaic tile type identifier (1 for land surface, 0 for inland lake) (also discussed [here](@ref compvsmosaic))
- *SDEP* Soil permeable depth [m]
  - The soil permeable depth, i.e. the depth to bedrock, may be less than the modelled thermal depth of the soil profile. If the depth to bedrock occurs within a soil layer, CLASS assigns the specified mineral or organic soil characteristics to the part of the layer above bedrock, and values corresponding to rock to the portion below. All layers fully below *SDEP* are treated as rock. These layers will be given a sand value of -3.
- *SOCI* Soil colour index
  - The soil colour index is used to assign the soil albedo.  It ranges from 1 to 20; low values indicate bright soils and high values indicate dark (see Lawrence and Chase (2007) \cite Lawrence2007-bc).  The wet and dry visible and near-infrared albedos for the given index are obtained from lookup tables, which can be found in @ref soilProperties.f90 .

CLASS provides a means of accounting for the possibility of the depth to bedrock falling within a layer, and therefore of phase changes of water taking place in only the upper part of the layer, by introducing the variable *TBAS*, which refers to the temperature of the lower part of the layer containing the bedrock. At the beginning of the time step the temperature of the upper part of the layer is disaggregated from the overall average layer temperature using the saved value of *TBAS*. The heat flow between the upper part of the soil layer and the lower part is diagnosed from the heat flux at the top of the layer. The upper layer temperature and *TBAS* are stepped ahead separately, and the net heat flux in the upper part of the layer is used in the phase change of water if appropriate. The upper layer temperature and *TBAS* are re-aggregated at the end of the time step to yield once again the overall average layer temperature.

Two variables, assumed to be constant over the grid cell, are provided if required for atmospheric model runs:

- *GGEO* Geothermal heat flux [\f$W m^{-2} \f$]
  - Unless the soil depth is very large and/or the run is very long, the geothermal heat flux can be set to zero. Since this is rarely used, this is not presently read in from the initialization file.
- *Z0OR* Orographic roughness length [m]
  - *Z0OR* is the surface roughness length representing the contribution of orography or other terrain effects to the overall roughness, which becomes important when the modelled grid cell is very large (e.g. in a GCM). For field studies it can be set to zero. It is presently not required as input for a model run and is set to zero in @ref model_state_drivers.read_initialstate

Four parameters are required for modelling lateral movement of soil water: *GRKF*, *WFCI*, *WFSF* and *XSLP*. However, the routines for interflow and streamflow modelling are not implemented in this version of CLASSIC, so they are not read in.

- *grclarea* Area of grid cell [\f$km^{2}\f$ ]

Area of the grid cell is required if CTEM is on. It is used for calculations relating to fire (@ref disturbance_scheme.disturb) and land use change (@ref landuse_change.luc). Both of those subroutines are not generally used (or meaningful) at the point scale so for point scale runs *grclarea* can be set to an arbitrary number like 100 \f$km^{2}\f$.

# Initialization of Prognostic Variables {#initProgVar}

## Initialization of Physics Prognostic Variables {#initPhysProgVar}

CLASSIC requires initial values of the land surface prognostic variables, either from the most recent atmospheric model integration or from field measurements. These are listed below, with guidelines for specifying values for each if measured values are not available.

- *TBAR* Temperature of soil layers [K]
- *THIC* Volumetric frozen water content of soil layers [\f$m^3 m^{-3}\f$ ]
- *THLQ* Volumetric liquid water content of soil layers [\f$m^3 m^{-3} \f$]

*TBAR*, *THLQ* and *THIC* are required for each of the modelled soil layers. Thin soil layers near the surface equilibrate quickly, but thicker, deeper layers respond more slowly, and long-term biases can be introduced into the simulation if their temperatures and moisture contents are not initialized to reasonable values or the model is not spun up for a suitably long period. If measured values are not available, for the moisture contents, it may be better to err on the low side, since soil moisture recharge typically takes place on shorter time scales than soil moisture loss. Field capacity is commonly used as an initial value. If the soil layer temperature is above freezing, the liquid moisture content would be set to the field capacity and the frozen moisture content to zero; if the layer temperature is below zero, the liquid moisture content would be set to the minimum value and the frozen moisture content to the field capacity minus the minimum value (*THLMIN*; see @ref soilProperties.f90). Very deep soil temperatures do not have a large effect on surface fluxes, but errors in their initial values can adversely affect hydrological simulations. For rock or ice sheet layers, *THLQ* and *THIC* should both be set to zero.

- *TBAS* Temperature of bedrock in the last permeable soil layer [K]

If CLASSIC is run with more than 3 soil layers, a separate calculation of *TBAS* is deemed unnecessary (see @ref **TMCALC.f**). If the model is run with only 3 layers then *TBAS* is used. Generally this can be initialized to the temperature of the third soil layer (as is presently done in @ref model_state_drivers.read_initialstate) and it is not required to be read in from the initialization file.

- *SNO* Mass of snow pack [\f$kg m^{-2}\f$ ]
- *TSNO* Snowpack temperature [K]
- *WSNO* Liquid water content of snow pack [\f$kg m^{-2} \f$]
- *ALBS* Snow albedo [ ]
- *RHOS* Density of snow [\f$kg m^{-3}\f$ ]

It is best to begin a simulation in snow-free conditions, so that the snow simulation can start from the simplest possible state where *SNO*, *TSNO*, *ALBS*, *RHOS* and *WSNO* are all initialized to zero. If erroneous values of the snow variables are specified as initial conditions, this can lead to a persistent bias in the land surface simulation. It can also lead to instability and resulting model crashes. These variables are, however, written to the restart file allowing one to restart the model from the same state (for example, after a spinup to begin a transient simulation). *WSNO* is not written or read-in from the restart file as its influence is small.

- *TAC* Temperature of air within vegetation canopy [K]
- *TCAN* Vegetation canopy temperature [K]
- *RCAN* Intercepted liquid water stored on canopy [\f$kg m^{-2} \f$]
- *SCAN* Intercepted frozen water stored on canopy [\f$kg m^{-2} \f$]
- *QAC* Specific humidity of air within vegetation canopy space [\f$kg kg^{-1} \f$]
- *CMAI* Aggregated mass of vegetation canopy [\f$kg m^{-2} \f$]

The vegetation canopy has a relatively small heat capacity and water storage capacity compared with the soil, so its temperature and intercepted water stores equilibrate quite quickly. *TCAN* and *TAC* can be initialized to the air temperature and *QAC* to the air specific humidity. *RCAN* and *SCAN* can be initialized to zero. *CMAI*, which is used only in the diagnostic energy balance check during the time step, can also be set to zero. At present the initialization file does not contain *TAC*, *QAC* or *CMAI*. All are initialized in @ref model_state_drivers.read_initialstate

- *GRO* Vegetation growth index [ ]
  - GROROT should be initialized to 1 during the growing season and to 0 otherwise.

- *TPND* Temperature of ponded water [K]
- *TSFS* Ground surface temperature over subarea [K]
- *ZPND* Depth of ponded water on surface [m]

Surface ponded water is a small term and is ephemeral in nature, so *ZPND* and *TPND* can both be initialized to zero. *TSFS* is included simply to provide a first guess for the surface temperature iteration in the next time step, so it can be initialized to an arbitrary value (presently it is not read in but rather set in @ref model_state_drivers.read_initialstate for the snow-covered subareas of the surface to the freezing point of water; for the snow-free subareas to the temperature of the first soil layer).

## Initialization of Biogeochemical Prognostic Variables {#initCTEMProgVar}

CTEM's initialization variables are generally PFT dependent. Several model options require extra initialization variables which are described in @ref CTEMaddInputs.

- *gleafmas* Green leaf mass [\f$kg C m^{-2} \f$]
- *bleafmas* Brown leaf mass [\f$kg C m^{-2} \f$]
  - Brown leaf mass is only for grasses, all of other PFTs have a value of zero.
- *rootmass* Root mass [\f$kg C m^{-2} \f$]
  - Root mass is assumed to include both the fine and coarse root masses of the plants. No distinction is made within the model.
- *stemmass* Stem mass [\f$kg C m^{-2} \f$]
  - Grasses and crops have no stem mass so are given values of zero.
- *soilcmas* Soil C mass per soil layer [\f$kg C m^{-2} \f$]
- *litrmass* Litter mass per soil layer [\f$kg C m^{-2} \f$]
  - Both *litrmass* and *soilcmas* are per PFT but also the arrays contain two additional values:
  - **bareground (icc + 1)**, where icc is the number of CTEM PFTs, and the **land use change (LUC) product pool (icc + 2)**. Soil C can be within the bare ground when land use change, competition or fire create bare ground that used to have vegetation. That C then continues to respire even though the vegetation have left that area. The LUC product pool represent the C that is converted into LUC products (soil C is 'furniture'- effectively a long lived C storage pool for the C removed from the landscape due to LUC whereas litter is 'paper', a short lived pool). The LUC C respires in response to environmental cues similar to litter and soil C in a manner similar to that material being in a landfill, for e.g.

Specify initialization amounts for green leaf, brown leaf, root and stem biomass, litter and soil carbon mass for each CTEM PFT. When growing vegetation from bare ground set these to zero or if no values are known.

- *lfstatus* Leaf pheological state [ ]
  - The phenology subroutine of CTEM tracks leaf status according to four plants states, 1) leaf onset or maximum growth, 2) normal growth, 3) leaf fall, and 4) no leaves state. When leaf status is set to 4 the model thinks there are no leaves. When initializing from bare ground with no vegetation set this to 4.
- *pandays* Days with positive net photosynthesis [Days]
  - The phenology subroutine of CTEM tracks days with positive net photosynthesis to initiate leaf onset. When net photosynthesis is positive for seven days, favourable weather is assumed to arrive and, leaf onset begins. If this variable is set to 7 then the model will think that favourable weather has arrived. When initializing from bare ground with no vegetation set this variable to 0.

# Example model setups {#exModSets}

  Here are a few example situations for physics only runs (CLASS alone).

  1. sum(*FCAN*) = 1 and *FARE* > 0 - Simulate all vegetated ground with no bare ground. If *FCAN(ican+1)* > 0 this includes some urban area.
  2. sum(*FCAN*) < 1 and *FARE* > 0  - Simulate vegetated ground and some bare ground (1 - sum(*FCAN*)). If *FCAN(ican+1)* > 0 this includes some urban area.
  3. sum(*FCAN*) = 0  and *FARE* > 0 - Only simulate bare ground.
  4. sum(*FCAN*) = 0  and *FARE* = 0 - Don't simulate this cell, it is a lake or something.

If biogeochemistry is turned on (CTEM on) then we have a subtle difference in that *fcancmx* is used and determines the *fcan(1:ican)* values. Assumedly if you have some vegetation then at least one ground layer should be soil (*SAND* > 0).
