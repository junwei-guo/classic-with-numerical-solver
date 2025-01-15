# How do I ... / Something has gone wrong! {#howDoI}

1. @ref chgSoil
2. @ref failStart
3. @ref makeClean
4. @ref soilColourIndex
5. @ref runInWindows
6. @ref windowsEndings
7. @ref izrefsite


----

# Change the number/depth/etc of soil layers? {#chgSoil}

Information about the soil layers is taken in from the initialization netcdf file. model_state_drivers::read_modelsetup reads the netcdf for the number of soil layers (*ignd*). The *DELZ* variable in the intialization file is the thickness of each layer. So change the initialization file for the number and thicknesses of the soil layers as required and CLASSIC will simply read them in.

# My run starts and I get an immediate fail! {#failStart}

If you start a new run and it immediately fails with an output something like:

        Singularity jormelton-containerCLASSIC-master-latest.simg:~/Documents/CLASSIC> bin/CLASSIC configurationFiles/template_job_options_file.txt 254.85/53.98
         in process           0           1           1   254.84999999999999        53.979999999999997                1           1           1           1
          CANOPY ENERGY BALANCE         1       8092.77620036       8077.49643236
                 0.000000       1.875648      -2.045030      -0.189964       0.000000    8073.385790
                 0.000000       1.339229     257.925000
         died on   254.84999999999999        53.979999999999997 
         
This is a general indication of a problem in your initialization file. Often it is fixed by setting snow/liquid in the canopy to zero (see @ref initPhysProgVar). You may also get snow energy/water balance failing, again look to @ref initPhysProgVar for advice with setup.

# My run won't compile but I didn't do anything I can think of. {#makeClean}

Sometimes you need to run a '`make mode=???? clean`' (where mode is the one you used to compile the code, for site-level it is typically `serial`) to clean out old .mod and .o files that might be causing issues. Do that and then recompile.

# initFileConverter tells me I need a soil colour index {#soilColourIndex}

Soil colour index is described in @ref soilData and used in soilProperties.f90. A global file of soil colour index produced by Peter Lawrence (NCAR) can be obtained from ftp://ftp.cccma.ec.gc.ca/pub/jmelton/mksrf_soilcol_global_c090324.nc

# Run on a Windows Machine {#runInWindows}

Courtesy of E. Humphreys. 

1. Install Oracle VirtualBox
2. Install Ubuntu as a virtual machine with ~ 100GB drive space
3. Install Dropbox on Ubuntu (to share files between machines but there are other ways to do this)
4. Install Singularity.  
  1. In Terminal, type: >sudo apt install singularity-container
5. Download our container into a working directory that is shared via Dropbox with the windows machine, e.g. ~/Dropbox/CLASSIC_working/:
  1. In Terminal, cd to your working directory and type: >singularity pull shub://jormelton/containerCLASSIC
6. Shell into the container:  >singularity shell jormelton-containerCLASSIC-master-latest.simg
7. In this container, get CLASSIC code and folders:  

**Note the text files need to have unix line endings if created in Windows (see @ref windowsEndings), make sure to convert to Unix (LF) from Windows (CR LF).**

# Gotcha for Windows users {#windowsEndings}

Text files created on DOS/Windows machines have different line endings than files created on Unix/Linux. DOS uses carriage return and line feed ("\r\n") as a line ending, which Unix uses just line feed ("\n"). You need to be careful about transferring files between Windows machines and Unix machines to make sure the line endings are translated properly.(from http://www.cs.toronto.edu/~krueger/csc209h/tut/line-endings.html)

Any CLASSIC tool that reads in ASCII expects the Linux/Unix line endings.

# My site-level run fails with IZREF = 1 but is okay with IZREF = 2 {#izrefsite}

Site-level runs should use IZREF = 1 with appropriate values of ZRFM and ZRFH for your site (described in @ref setupJobOpts). However it is possible, especially with CTEM (biogeochemistry) on that the vegetation will grow taller than ZRFM/ZRFH causing the model to crash. The model may simulate vegetation at the site taller than reality due to prior land use history at the site, model bias, etc. In this case you will need to raise the ZRFM/ZRFH heights to remain above the canopy. 