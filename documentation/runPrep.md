
# Preparing a CLASSIC run {#runPrep}

1. @ref Environ
  1. @ref Containers
2. @ref compilingMod
3. @ref setupJobOpts
4. @ref xmlSystem

----

# Setting up the runtime environment {#Environ}

To run CLASSIC you can either use your own immediate environment or use our Singularity container (**RECOMMENDED** and described below). If you use your own immediate environment, the following Linux packages are required at a minimum:

- make
- libnetcdff-dev
- git
- gfortran
- netcdf-bin
- zlib1g

These packages will allow serial compiling and running of the model. To run in parallel add mpich. To run the documentation tool add doxygen.

Note that the above packages are available for Ubuntu/Debian and may have different names under other Linux distros. 

# Running CLASSIC in a Singularity Container {#Containers}

## What are containers?

Containers are a tool to package and distribute all the elements and dependencies of a Linux-based application. Within a container it is possible to store the computing environment along with applications such as model code. Containers bring an ease of transportation, installation, and execution across operating systems such as Linux (local or cloud), Mac, and Windows.
[Docker is a type of container](https://www.docker.com/). [Singularity](https://www.sylabs.io/) is a platform that will allow us to build, run or shell into one of these Docker containers.

In order to run CLASSIC, we need several specific software tools such as compilers (e.g. GNU, Intel, Cray, etc.) and libraries (MPI, NetCDF etc.) that have to be present on our machine to be able to run the model. If we use a container then we can simply access the model through the container, a process that eliminates the cumbersome and time consuming process of installing all of these compilers and libraries locally.

Another significant advantage is that the versions of each library is "frozen" in the container. This has several advantages for scientific reproducibility, mobility, and model development. For example, we could have a container with the library versions at the time of a major model release version. So, even several years later, we can run the model with the original intended library versions recreating the software environment exactly. A version of the WRF weather prediction model has been containerized to aid in it use in teaching and research (Hacker et al. 2016) \cite Hacker2016-qg .

## Benefits of {Singularity} containers

From Kurtzer et al. (2017) \cite Kurtzer2017-xc :

> Singularity offers mobility of compute by enabling environments to be completely portable via a single image file,  and is designed with the features necessary to allow seamless integration with any scientific computational resources. ... Mobility of compute is defined as the ability to define, create, and maintain a workflow locally while remaining confident that the workflow can be executed on different hosts, Linux operating systems, and/or cloud service providers. In essence, mobility of compute means being able to contain the entire software stack, from data files up through the library stack, and reliability move it from system to system. Mobility of compute is an essential building block for reproducible science, and consistent and continuous deployment of applications. ... Many of the same features that facilitate mobility also facilitate reproducibility. Once a contained workflow has been defined, the container image can be snapshotted, archived, and locked down such that it can be used later and the user can be confident that the code within the container has not changed. The container is not subject to any external influence from the host operating system (aside from the kernel which is ubiquitous of any OS level virtualization solution).... Singularity can give the user the freedom they need to install the applications, versions, and dependencies for their workflows without impacting the system in any way. Users can define their own working environment and literally copy that environment image (a single file) to a shared resource, and run their workflow inside that image.

A nice benefit of containers is that they are designed to be easy to use (Kurtzer et al., 2017) \cite Kurtzer2017-xc :

> The goal of Singularity is to support existing and traditional HPC resources as easily as installing a single package onto the host operating system. For the administrators of the hosts, some configuration may be required via a single configuration file, however the default values are tuned to be generally applicable for shared environments.

System administrators will be comforted to know (Kurtzer et al., 2017) \cite Kurtzer2017-xc :

> Singularity does not provide a pathway for privilege escalation (which makes it truly applicable for multi-tenant shared scientific compute resources). This means that in the runtime environment, a user inside a Singularity container is the same user as outside the container. If a user wants to be root inside the container, they must first become root outside the container. Considering on most shared resources the user will not have root access means they will not have root access within their containers either. This simple concept thus defines the Singularity usage workflow.

## How to use Singularity containers?

In order to use Singularity containers, one must first make certain that a local installation of Singularity is available.

For most up-to-date instructions on installing Singularity on Linux, Mac, or Windows see the [Singularity documentation](https://www.sylabs.io/docs/).

Generally, on a Linux machine (Ubuntu in our particular case), one may use aptitude with the following command to install singularity (providing the user has administrative privelidges).

`sudo apt install singularity-container`

## Obtaining the CLASSIC Singularity container

Download our container from our CLASSIC community Zenodo page: https://zenodo.org/communities/classic. It is recommended you follow the [Quick Start Tutorial](https://cccma.gitlab.io/classic_pages/info/get_started/) if possible.

You can either use the `exec` command form or shell into the container (example assumes you are in the same folder as the .simg file):

`singularity shell CLASSIC_container.simg`

E.g.

        acrnrjm@cccsing: ~> singularity shell /user/nphome1/rjm/CLASSIC_container.simg
        Singularity: Invoking an interactive shell within container...

        Singularity CLASSIC_container.simg:~>

If that is successful, you are now in the CLASSIC container environment. This environment contains all the libraries needed to run the model (Note this is a bare-bones installation with only the run-time environment. It does not presently contain a workflow or compiled model code - see the [Quick Start Tutorial](https://cccma.gitlab.io/classic_pages/info/get_started/).

E.g. test if gfortran is installed:

        Singularity CLASSIC_container.simg:~> gfortran
        gfortran: fatal error: no input files
        compilation terminated.
        Singularity CLASSIC_container.simg:~>

And test for something that we know is not installed:

        Singularity CLASSIC_container.simg:~> okular
        bash: okular: command not found
        Singularity CLASSIC_container.simg:~>

You can now navigate to the location of CLASSIC code, compile and run the model.

E.g.

        Singularity CLASSIC_container.simg:~/Documents/CLASSIC> bin/CLASSIC
         Usage is as follows

         bin/CLASSIC joboptions_file longitude/{longitude}/latitude/{latitude}

         - joboptions_file - an example is
           configurationFiles/template_job_options_file.txt.

         - longitude/latitude
           e.g. 105.23/40.91

          *OR*
          if you wish to run a region then you give
          the corners of the box you wish to run

         - longitude/longitude/latitude/latitude
           e.g. 90/105/30/45

## Dealing with access to file systems while in a container

From the [Singularity manual](https://www.sylabs.io/guides/2.6/user-guide/bind_paths_and_mounts.html?highlight=bind%20mount):

> When Singularity ‘swaps’ the host operating system for the one inside your container, the host file systems becomes inaccessible. But you may want to read and write files on the host system from within the container. To enable this functionality, Singularity will bind directories back in via two primary methods: system-defined bind points and conditional user-defined bind points.

This can demonstrated by an example in which a remote server is visible on the host operating system but not within the Singularity container. Using the -B option binds that path allowing it to now be accessible:

        acrnrjm@cccsing: ~/Documents/CLASSIC> singularity shell ../../CLASSIC_container.simg 
        Singularity: Invoking an interactive shell within container...

        Singularity CLASSIC_container.simg:~/Documents/CLASSIC> ls /raid/ra40/data/rjm/meteorologicalDatasets/CRU_JRA_v1.0.5_1901_2017
        ls: cannot access '/raid/ra40/data/rjm/meteorologicalDatasets/CRU_JRA_v1.0.5_1901_2017': No such file or directory

        Singularity CLASSIC_container.simg:~/Documents/CLASSIC> exit
        exit

        acrnrjm@cccsing: ~/Documents/CLASSIC> singularity shell -B /raid/ra40/data/rjm/meteorologicalDatasets/CRU_JRA_v1.0.5_1901_2017 ../../CLASSIC_container.simg 
        Singularity: Invoking an interactive shell within container...

        Singularity CLASSIC_container.simg:~/Documents/CLASSIC> ls /raid/ra40/data/rjm/meteorologicalDatasets/CRU_JRA_v1.0.5_1901_2017
        dlwrf_T63_chunked_1700_2017.nc dswrf_v1.1.5_T63_chunked_1700_2017.nc spfh_T63_chunked_1700_2017.nc ...
        
This -B flag can be very useful when using a Vagrant box to run CLASSIC on a Windows, Mac or even Linux machine. The -B flag can point to the [synced folder](https://www.vagrantup.com/docs/synced-folders/basic_usage.html).

# Compiling CLASSIC for serial and parallel simulations {#compilingMod}

CLASSIC's Makefile facilitates the compilation of the model code for running in serial mode or in parallel.

The basic command "make" (without options) compiles the code for running in serial mode. The executable placed in the bin directory is then 'CLASSIC_serial'

The "mode" option can be specified on the "make" command line to target different platforms and/or compilation modes as follows:

make mode=<mode>

where <mode> is:
  
        serial        - compiles code for serial runs (default)
        parallel      - compiles code for parallel runs

The serial compilation is presently set to use the GNU compiler (gfortran). 

Parallel compilation uses the MPI wrapper "mpif90". Use "mpif90 -show" to display the pre-configured compiler command for the system. The Makefile is currently set to use the GNU compiler with "mpif90". The settings for the compiler as well as the include/library flags are system dependent and will likely require modification on each system.

There are two additional modes that are particularly useful for users with access to the ECCC supercomputer environment. These are:

        ppp           - compiles code for parallel runs using the Intel compiler on the ECCC pre/post-processing clusters
        supercomputer - compiles code for parallel runs using the Intel compiler on the ECCC supercomputers

These modes are unlikely to be useful for non-ECCC users.

Upon compilation, all object (.o) and module (.mod) files are placed into an "objectFiles" directory. These directories are labelled with an additional suffix according to the compilation mode used. Executables are placed in the "bin" directory, again labelled according to the compilation mode.

A useful command is 'make mode={serial,parallel,etc.} clean', which removes all *.o *.mod files (for a specified mode). This can allow a fresh compilation which can be handy if some parameters are changed that aren't being refreshed on a make. Note that a re-compilation will overwrite an existing executable of the same name (i.e. mode) in the "bin" directory.

# Setting up the joboptions file {#setupJobOpts}

The joboptions file controls the model configuration, inputs files used, and model outputs. The template joboptions file is located in the configurationFiles folder. Use this as your starting point.

If, for example, we wanted to run CLASSIC at a site with observed meteorology, with leap years, from 1991 to 2017 then we can set up the meteorological options as follows:

        &joboptions

        ! Meteorological options:
            readMetStartYear = 1991  ! First year of meteorological forcing to read in from the met file
            readMetEndYear = 2017    ! Last year of meteorological forcing to read in from the met file
            metLoop = 1 ,            ! no. of times to cycle over the read-in meteorology
            leap = .true. ,         ! True if your meteorological forcing includes leap years

If we are interested in spinning up the C pools, we could set metLoop to run the model over the readMetStartYear to readMetEndYear a metLoop number of times. We should also point to our meteorological forcing files like,

        ! Meteorological forcing files:
            metFileFss = pathToFile/dswrf.nc',        ! location of the incoming shortwave radiation meteorology file
            metFileFdl = 'pathToFile/dlwrf.nc',        ! location of the incoming longwave radiation meteorology file
            metFilePre = 'pathToFile/pre.nc',        ! location of the precipitation meteorology file
            metFileTa = 'pathToFile/tmp.nc',         ! location of the air temperature meteorology file
            metFileQa = 'pathToFile/spfh.nc',         ! location of the specific humidity meteorology file
            metFileUv = 'pathToFile/wind.nc',         ! location of the wind speed meteorology file
            metFilePres = 'pathToFile/pres.nc',       ! location of the atmospheric pressure meteorology file

The model initialization and restart files need to be pointed to. Note the rs_file_to_overwrite is **overwritten**. The simplest thing to do at a start of a run is to make a copy of the init_file to be the rs_file_to_overwrite. CLASSIC only overwrites the prognostic variables leaving any others unchanged (see model_state_drivers.write_restart). The model parameters file also needs to be pointed to.

        ! Initialization and restart files
            init_file = 'inputFiles/CLASSCTEM_initialization_bulkdetrital_nosnow.nc' ,     ! location of the model initialization file
            rs_file_to_overwrite = 'rsFile.nc' ,       ! location of the existing netcdf file that will be **overwritten** for the restart file
                                                       !! typically here you will just copy the init_file and rename it so it can be overwritten.
        ! Namelist of model parameters
            runparams_file = 'configurationFiles/template_run_parameters.txt' ,    ! location of the namelist file containing the model parameters

The next series of switches relate to CTEM. If you wish to run a physics only run with prescribed vegetation parameters, ctem_on is set to .false. and the remainder of this section is ignored. The physics switches are farther down.

        ! CTEM (biogeochemistry) switches:
         ctem_on = .true. ,     ! set this to true for using ctem simulated dynamic lai and canopy mass, else class simulated specified
                                ! lai and canopy mass are used. with this switch on, all the main ctem subroutines are run.
            spinfast = 3 ,      ! Set this to a higher number up to 10 to spin up soil carbon pool faster. Set to 1 for final round of spin up and transient runs.

The CO2 and CH4 switches behave similarly. If a constant [CO2] is desired, transientCO2 is set to .false. and the year of observed CO2 is specified in fixedYearCO2. The [CO2] of the corresponding year is selected from the CO2File and used for the run. Simlarly for CH4. If are **not** interested in simulating methane related variables then simply copy your CO2File to be your CH4file (e.g. 'cp CO2file.nc fakeCH4file.nc', and then specify the fakeCH4file.nc for the CH4File since the model expects a file there). If for some reason you wish to run with non-historical [CO2], you could edit your CO2File (using a combination of the NCO tools ncdump and ncgen, for example) to include a future year with associated CO2 value (like 2050 and 500ppm for example)

        !CO2 switches:
            transientCO2 = .true. , ! Read in time varying CO2 concentration from CO2File or if set to false then use only the year of fixedYearCO2 value
            CO2File = 'inputFiles/mole_fraction_of_carbon_dioxide_in_air_input4MIPs_GHGConcentrations_CMIP_UoM-CMIP-1-2-0_gr1-GMNHSH_1850-2014_absTimeIndex.nc' ,
            fixedYearCO2 = 1850 ,   ! If transientCO2 = .true., this is ignored.

        !CH4 switches:
            ! If you don't care about running with methane, you can just copy your CO2 file (cp CO2file.nc fakeCH4file.nc) and point to it here.
            ! It won't affect the rest of your run. However you do need to have a file specified if CTEM is on.
            transientCH4 = .true. , ! Read in time varying CH4 concentration from CH4File or if set to false then use only the year of fixedYearCH4 value
            CH4File = 'inputFiles/mole_fraction_of_methane_in_air_input4MIPs_GHGConcentrations_CMIP_UoM-CMIP-1-2-0_gr1-GMNHSH_1850-2014_absTimeAxis.nc',
            fixedYearCH4 = 1850 ,   !If transientCH4 = .true., this is ignored.

Disturbance in the form of fire is optional in CLASSIC. If fire is turned on then the model requires a population density input file (see @ref initPopd) which is used in a manner similar to CO2 and CH4, i.e. it can use a fixed or transient value. Lightning strikes (see @ref initLightFire) are also required and can be specified similarly to population density.

        !Fire switches
            dofire = .true. ,               ! If true the fire disturbance parameterization is turned on.

                transientPOPD = .true. ,    ! Read in time varying population density from POPDFile or if set to false then use only the year of fixedYearPOPD.
                POPDFile = 'inputFiles/POPD_annual_1850_2016_T63_Sep_2017_chunked.nc' ,
                fixedYearPOPD = 1850 ,      ! If transientPOPD = .true., this is ignored.

                transientLGHT= .true.      ! use lightning strike time series, otherwise use fixedYearLGHT
                LGHTFile = 'inputFiles/lisotd_1995_2014_climtlgl_lghtng_as_ts_1850_2050_chunked.nc' , ! Location of the netcdf file containing lightning strike values
                fixedYearLGHT = 1850 ,    ! set the year to use for lightning strikes if transientLGHT is false.

Competition for space between plant functional types is parameterized in CLASSIC. If PFTCompetition is set to false, the PFT fractional coverage follows rules as will be outlined next. If PFTCompetition is true then CLASSIC has two more switches of interest. Competition uses bioclimatic indices to determine whether PFTs should be allowed to compete within a gridcell (see competition_scheme.existence and @ref initClimComp). The bioclimatic indices are either read-in from the init_file or are calculated anew for the run underway. If the values are in the init_file from a spinup run then set inibioclim to true, otherwise set inibioclim to false. It is also possible to start the model from bare ground (rather than the PFT configuration found in your init_file) by setting start_bare to true.

        ! Competition switches:
            PFTCompetition = .false. ,      ! If true, competition between PFTs for space on a grid cell is implimented
                inibioclim = .false. ,      ! set this to true if competition between pfts is to be implimented and you have the mean climate values
                                            ! in the init netcdf file.
                start_bare = .false.,       ! Set this to true if competition is true, and if you wish to start from bare ground. if this is set to false, the
                                            ! init netcdf file info will be used to set up the run. NOTE: This still keeps the crop fractions
                                            ! (while setting all pools to zero)

Land use change is possible via a LUCFile (see @ref inputLUC) that has the fractional coverage for each PFT annually. This switch has an additional option as described in the comment below. Pay attention here.

        ! Land Use switches:
            !** If you wish to use your own PFT fractional covers (specified in the init_file), set fixedYearLUC to -9999, otherwise set it
            ! to the year of land cover you want to use. If you wish to have transient land cover changes, set
            ! lnduseon to true, it will update the fractional coverages from LUCFile. When lnduseon is false it is
            ! not updated beyond the initial read in of landcover for fixedYearLUC, or if -9999 then the LUCFile is
            ! not used at all.**
            lnduseon = .true. ,
            LUCFile = 'inputFiles/LUH_HYDE_based_crop_area_adjusted_land_cover_CTEM_fractions_1850_2017_T63_chunked.nc' ,
            fixedYearLUC = 1901 ,

CLASSIC can determine dynamics wetland locations for wetland methane emissions. Alternatively CLASSIC can read in time evolving wetland fractions from an external file.

        ! Wetland switches:
            ! If you wish to read in and use observed wetland fractions, there are two options. If you wish time
            ! evolving wetland fractions set transientOBSWETF to true and give a OBSWETFFile. If you wish to use
            ! a single year of that file set transientOBSWETF to false, give a OBSWETFFile, and set fixedYearOBSWETF
            ! to some valid year. If you wish to use only dynamically determined wetland fractions set transientOBSWETF
            ! to false and set fixedYearOBSWETF to -9999. The slope fractions in the init_file will then be used to
            ! dynamically determine wetland extent.
            transientOBSWETF = .false. ,  ! use observed wetland fraction time series, otherwise use fixedYearOBSWETF
            OBSWETFFile = '',             ! Location of the netcdf file containing observed wetland fraction
            fixedYearOBSWETF = -9999 ,    ! set the year to use for observed wetland fraction if transientOBSWETF is false.

CLASS switches determine the configuration of the physics only as well as CLASS+CTEM (physics and biogeochemistry) runs.

        ! Physics switches:

            IDISP = 0 ,    !  This switch controls the calculation of the vegetation        displacement height. In most
                          !atmospheric models a “terrain-following” coordinate system is used, in which the vegetation
                          !displacement height is considered to be part of the “terrain”, and is therefore neglected. 
                          ! For such applications IDISP is set to 0. For studies making use of field data, it should !be set to 1.
            IZREF = 2 ,    !This switch indicates where the bottom of the atmosphere is conceptually located. In
                          !most atmospheric models the bottom is assumed to lie at the local surface roughness length, 
                          ! i.e. where the horizontal wind speed is zero; for such simulations IZREF is set to 2. For 
                          ! all other cases, including sites using field data, it should be set to 1.
            ZRFH = 40.0,  ! The reference heights at which the energy variables (air temperature and specific humidity) are provided. When 
                          ! using field data this is the measurement height. 
            ZRFM = 40.0,  ! The reference heights at which the momentum variables (wind speed) are provided. When using field data this is the measurement height. 
            ZBLD = 50.0,  ! The atmospheric blending height.  Technically this variable depends on the length scale of the
                          ! patches of roughness elements on the land surface, but this is difficult to ascertain.  Usually it is assigned a value of 50 m.
                            
            ISLFD = 0 ,   ! This switch indicates which surface layer flux stability
                          ! correction and screen-level diagnostic subroutines are to be
                          ! used. If ISLFD=0, the CCCma stability correction subroutine
                          ! DRCOEF is used with simple similarity-based screen-level
                          ! diagnostic calculations. If ISLFD=1, DRCOEF is used for the
                          ! stability corrections and the RPN subroutine SLDIAG is used
                          ! for the screen level diagnostic calculations. If ISLFD=2, the
                          ! RPN stability correction subroutine FLXSURFZ is used with the
                          ! companion RPN subroutine DIASURFZ for the screen level
                          ! diagnostics. (When running CLASS coupled to an atmospheric 
                          ! model with staggered vertical thermodynamic and momentum
                          ! levels, ISLFD must be set to 2, since FLXSURFZ allows inputs
                          ! on staggered levels and DRCOEF does not.)

! The implications of the ISLFD switch are discussed more in CLASST.f

            IPCP = 1 ,     ! if ipcp=1, the rainfall-snowfall cutoff is taken to lie at 0 C. if ipcp=2, a linear partitioning of precipitation between
                            ! rainfall and snowfall is done between 0 C and 2 C. if ipcp=3, rainfall and snowfall are partitioned according to
                            ! a polynomial curve between 0 C and 6 C.
                            
This switch is used to specify which approach is to be used in subroutine atmosphericVarsCalc.f90 for the
partitioning of precipitation between rainfall and snowfall. If IPCP=1, precipitation is assumed to be all rainfall at air temperatures above 0 C and all snowfall at air temperatures ≤ 0 C. If IPCP=2,the partitioning between rainfall and snowfall varies linearly between all rainfall at temperatures above 2 C, and all snowfall at temperatures below 0 C. If IPCP=3, rainfall and snowfall are partitioned according to a polynomial curve, from all rainfall at temperatures above 6 C to all snowfall at temperatures below 0 C. If IPCP=4, the rainfall and snowfall rates, *RPRE* and *SPRE*, are read in directly at each time step, and the total precipitation is calculated as their sum. (Note that if the latter option is used, the model will need to be adapted to read the extra variables in)
                            
            IWF = 0 ,      ! if iwf=0, only overland flow and baseflow are modelled, and the ground surface slope is not modelled. if iwf=n (0<n<4) ,
                           ! the watflood calculations of overland flow and interflow are performed; interflow is drawn from the top n soil layers.
                           
**Note** the interflow and streamflow parts of the code are still under development, so unless the user is engaged in this development, this switch should be set to 0.
                           
            isnoalb = 0 ,  ! if isnoalb is set to 0, the original two-band snow albedo algorithms are used. if it is set to 1, the new four-band routines are used.
            
**Note** the four-band albedo, while implemented within the code, is still being tested. Users are advised to use isnoalb = 0.

The iteration scheme for canopy or ground surface temperatures can be either a bisection or Newton-Raphson method. **Note:** Recently problems have been discovered with the Newton-Raphson scheme, involving instabilities and occasional failure to converge, so currently users are advised not to select this option.

         ! Iteration scheme
            ! ITC, ITCG and ITG are switches to choose the iteration scheme to be used in calculating the canopy or ground surface temperature
            ! respectively.  if the switch is set to 1, a bisection method is used; if to 2, the newton-raphson method is used.
            ITC = 1 ,   ! Canopy
            ITCG = 1 ,  ! Ground under canopy
            ITG = 1 ,   ! Ground

User supplied values can be used for plant area index, vegetation height, and canopy, soil, or snow albedos. If any of these inputs are supplied the model_state_drivers.f90 needs to be adapted to read-in the user-supplied values.

         ! User-supplied values:
            ! if ipai, ihgt, ialc, ials and ialg are zero, the values of plant area index, vegetation height, canopy albedo, snow albedo
            ! and soil albedo respectively calculated by class are used. if any of these switches is set to 1, the value of the
            ! corresponding parameter calculated by class is overridden by a user-supplied input value.
            IPAI = 0 ,  ! Plant area index
            IHGT = 0 ,  ! Vegetation height
            IALC = 0 ,  ! Canopy albedo
            IALS = 0 ,  ! Snow albedo
            IALG = 0 ,  ! Soil albedo

Model outputs are in netcdf format. The outputs metadata is read in from an xmlFile (see @ref xmlSystem). Outputs can be grid-cell average (default), per PFT, or per tile. In all cases the output has to be specified in the xmlFile and also properly handled in prepareOutputs.f90. Temporal resolution of output files are half-hourly (for physics variables as well as photosynthesis and canopy conductance only), daily, monthly and annually. For half-hourly and daily outputs the start and end days as well as start and end years can be specified. Monthly file can specify the year to start writing the outputs.

        ! Output options:

            output_directory = 'outputFiles' ,        ! Directory where the output netcdfs will be placed
            xmlFile = 'configurationFiles/outputVariableDescriptors_v1.2.xml' ,  ! location of the xml file that outlines the possible netcdf output files

            doperpftoutput = .true. ,   ! Switch for making extra output files that are at the per PFT level
            dopertileoutput = .false. , ! Switch for making extra output files that are at the per tile level

            dohhoutput = .false. ,      ! Switch for making half hourly output files (annual are always outputted)
            JHHSTD = 166 ,                ! day of the year to start writing the half-hourly output
            JHHENDD = 185 ,             ! day of the year to stop writing the half-hourly output
            JHHSTY = 1901 ,             ! simulation year (iyear) to start writing the half-hourly output
            JHHENDY = 1901 ,            ! simulation year (iyear) to stop writing the half-hourly output

            dodayoutput = .false. ,     ! Switch for making daily output files (annual are always outputted)
            JDSTD = 20 ,                 ! day of the year to start writing the daily output
            JDENDD = 30 ,              ! day of the year to stop writing the daily output
            JDSTY = 1902 ,              ! simulation year (iyear) to start writing the daily output
            JDENDY = 1903 ,             ! simulation year (iyear) to stop writing the daily output

            domonthoutput = .true. ,    ! Switch for making monthly output files (annual are always outputted)
            JMOSTY = 1901 ,             ! Year to start writing out the monthly output files.
            
            doAnnualOutput = .true. ,   ! Switch for making annual output files 

CLASSIC has the capability to produce checksums to ensure model outputs haven't changed. The checksum module takes the model restart variables and calculates a sum based upon the binary representation of the variables' values. This 'checksum' can then be compared to previous runs to see if the checksum values have changed. While this test is not infalliable, the probability of a false positive decreases exponentially with the number of outputs being tested. For further info see https://en.wikipedia.org/wiki/Checksum. See @ref regTest for framework and checksum.f90 for specifics.

            doChecksums = .false.,      ! checksums can be generated if you wish to ensure you are making no
                                        ! functional changes to the model results. See the CLASSIC documentation and
                                        ! the regression_testing tool.

Comments can be added to output files using the Comment field below. Also comments can be left in the joboptions file after the backslash.

            Comment = ' test '          ! Comment about the run that will be written to the output netcdfs

         /

        This area can be used for comments about the run.
