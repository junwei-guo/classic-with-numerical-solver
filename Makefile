# Makefile for CLASSIC.

# Usage:

# make mode=<mode>

#   or, for the default mode, just:

# make

# where <mode> is:
#
# serial        - compiles code for serial runs using the GNU compiler on the local platform (default)
# parallel      - compiles code for parallel runs using the GNU compiler on the local platform
# supercomputer - compiles code for parallel runs using the Intel compiler on the ECCC supercomputers
# ppp           - compiles code for parallel runs using the Intel compiler on the ECCC front-ends
#                 *** NB: Before running make, support for the parallel netCDF library must first be enabled 
#                         using this command:
#                         . ssmuse-sh -x hpco/exp/hdf5-netcdf4/parallel/openmpi-3.1.2/static/intel-19.0.3.199/01

# The "cray" mode uses the native Cray compiler on the supercomputers. It still requires testing.

# ==================================================================================================================

# Modifications: Ed Chan, Dec 2019.
#   - Revised supercomputer/ppp compiler/linker options, which are now the same on both platforms.
#   - On ppp, support for the parallel netCDF library must now be enabled PRIOR to running make.
# Modifications: Ed Chan, Sep 2018.
#   - Updated various compiler options.
#   - Added "cray" and "ppp" modes.
#   - All objects go into directories specific to the mode used.
#   - All executable names are labelled with the mode appended.

# Object files
OBJ = brent.o numSolver.o classicParams.o ctemStateVars.o classStateVars.o peatlandsMod.o \
	generalUtils.o fileIOModule.o ctemUtilities.o calcLandSurfParams.o groundAlbedo.o mvidx.o fourBandAlbedo.o \
	SNOW_ALBVAL.o SNOW_TRANVAL.o snowAlbedoTransmiss.o soilHeatFluxPrep.o waterInfiltrateUnsat.o \
	canopyInterception.o waterFlowNonInfiltrate.o snowAging.o  energyBudgetPrep.o waterInfiltrateSat.o canopyAlbedoTransmiss.o \
	classGatherScatter.o waterFlowInfiltrate.o snowSublimation.o energBalVegSolve.o waterCalcPrep.o methaneProcesses.o \
	ctemGatherScatter.o applyAllometry.o photosynCanopyConduct.o canopyWaterUpdate.o \
	atmosphericVarsCalc.o canopyPhaseChange.o iceSheetBalance.o waterUnderCanopy.o energBalNoVegSolve.o errorHandler.o classGrowthIndex.o \
	DIASURFZ.o screenRelativeHumidity.o pondedWaterFreeze.o snowTempUpdate.o checkWaterBudget.o energyBudgetDriver.o DRCOEF.o SLDIAG.o waterUpdates.o snowHeatCond.o \
	radiationDriver.o waterBudgetDriver.o FLXSURFZ.o snowInfiltrateRipen.o snowMelt.o soilWaterPhaseChg.o soilProperties.o energyWaterBalanceCheck.o snowAddNew.o \
	soilHeatFluxCleanup.o waterBaseflow.o  balanceCarbon.o autotrophicRespiration.o phenolgy.o  \
	turnover.o mortality.o disturb.o competitionMod.o landuseChangeMod.o \
	allocateCarbon.o heterotrophicRespirationMod.o tracer.o ctemDriver.o outputManager.o \
	prepareOutputs.o  modelStateDrivers.o readFromJobOptions.o metModule.o main.o xmlParser.o \
	xmlManager.o CLASSIC.o
	# COMBAK PERLAY
	#allocateCarbon.o heterotrophicRespiration.o soilCProcesses.o tracer.o ctemDriver.o outputManager.o \
	prepareOutputs.o  modelStateDrivers.o readFromJobOptions.o metModule.o main.o xmlParser.o \
	xmlManager.o CLASSIC.o
	# COMBAK PERLAY

# Object files (.o and .mod) are placed into directories named according to the mode. In the serial
# compiler case, these are placed into the "objectFiles" directory.
ODIR = objectFiles_$(mode)

# Set variables based on the "mode" specified.

# *** NB: These settings are particular to the ECCC platforms. Other users will need to find their own recipes.
#         The use most general option is "serial".

ifeq ($(mode), supercomputer)
	# Wrapper to the default Fortran compiler loaded via a module. The following is specific to the Intel compiler.
	# Include/library flags should NOT be specified as the appropriate ones should be available via the loaded module.
	COMPILER = ftn
	# Fortran Flags.
        FFLAGS = -DPARALLEL -r8 -g -O2 -mp1 -xCORE-AVX2 -align array64byte -init=arrays -init=zero -traceback -module $(ODIR)
else ifeq ($(mode), cray)
	# Wrapper to the default Fortran compiler loaded via a module. The following is specific to the Cray compiler.
	# Include/library flags should NOT be specified as the appropriate ones should be available via the loaded module.
	COMPILER = ftn
	# Fortran Flags. Note: -g prevents most optimizations, -rm gives listings.
	FFLAGS = -DPARALLEL -s real64 -e0 -ez -O2 -g
else ifeq ($(mode), ppp)
        # NB: Support for parallel netCDF library must first be enabled prior to running make via:
        #     . ssmuse-sh -x hpco/exp/hdf5-netcdf4/parallel/openmpi-3.1.2/static/intel-19.0.3.199/01
        # Open MPI wrapper to the default Fortran compiler. The following is specific to the Intel compiler.
        # Include/library flags for Open MPI should NOT be specified as the appropriate ones should be available via the wrapper.
        COMPILER = mpif90
        # Fortran Flags.
        FFLAGS = -DPARALLEL -r8 -g -O2 -mp1 -xCORE-AVX2 -align array64byte -init=arrays -init=zero -traceback -module $(ODIR)
        # For debugging:
        #FFLAGS = -DPARALLEL -r8 -g -O1 -mp1 -xCORE-AVX2 -align array64byte -init=arrays -init=zero -traceback -fpe0 -module $(ODIR)
        # Library Flags for netCDF.
        LFLAGS = -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz -lcurl
else ifeq ($(mode), parallel)
        # Parallel compiler.
        COMPILER = mpif90
        # Fortran Flags. The following is specific to the gfortran compiler.
        FFLAGS = -DPARALLEL -O3 -g -fdefault-real-8 -ffree-line-length-none -fbacktrace -ffpe-trap=invalid,zero,overflow -fbounds-check -J$(ODIR)
        # Include Flags.
        IFLAGS =  -I/para_netcdf_hdf-4.6.3/MPI/include
        # Library Flags.
        LFLAGS = -L/para_netcdf_hdf-4.6.3/MPI/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5
else
	# Serial compiler.
	COMPILER = gfortran
	mode = serial
	ODIR = objectFiles
	# Fortran Flags.
	#FFLAGS = -O0 -g -fdefault-real-8 -ffree-line-length-none -fbacktrace -ffpe-trap=invalid,zero,overflow -fbounds-check -J$(ODIR) #-Wall -Wextra
	FFLAGS = -O3 -fdefault-real-8 -ffree-line-length-none -fbacktrace -ffpe-trap=invalid,zero,overflow -fbounds-check -J$(ODIR) #-Wall -Wextra
	# Include Flags.
	IFLAGS = -I/usr/include -I/usr/lib64/gfortran/modules
	# Library Flags
	LFLAGS = -lnetcdff -ldl -lz -lm
endif

# Create required directory/.gitignore file, if missing.
VOID := $(shell mkdir -p $(ODIR))
#VOID := $(shell [ ! -f $(ODIR)/.gitignore ] && cp objectFiles/.gitignore $(ODIR))

# RECIPES
# Compile object files from .F90 sources
$(ODIR)/%.o: src/%.F90
	$(COMPILER) $(FFLAGS) $(IFLAGS) -c $< -o $@

# Compile object files from .f90 sources
$(ODIR)/%.o: src/%.f90
	$(COMPILER) $(FFLAGS) $(IFLAGS) -c $< -o $@

# Compile object files from .f (Fortran 77) sources
$(ODIR)/%.o: src/%.f
	$(COMPILER) $(FFLAGS) $(IFLAGS) -c $< -o $@

# Properly reference the ODIR for the linking
OBJD = $(patsubst %,$(ODIR)/%,$(OBJ))

# Link objects together and put executable in the bin/ directory
CLASSIC: $(OBJD)
	$(COMPILER) $(FFLAGS) $(IFLAGS) -o bin/CLASSIC_$(mode) $(OBJD) $(LFLAGS)

# "make clean mode=supercomputer" removes all object files in objectFiles_supercomputer
clean:
	rm -f $(ODIR)/*.o $(ODIR)/*.mod bin/CLASSIC_$(mode)
