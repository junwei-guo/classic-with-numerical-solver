# Overview of the Canadian Land Surface Scheme (CLASS) {#overviewCLASS}

The Canadian Land Surface Scheme, CLASS, was originally developed for use with the Canadian Global Climate Model (CanGCM) (Verseghy, 1991 \cite Verseghy1991-635 ; Verseghy et al., 1993 \cite Verseghy1993-1ee ). The table at the end of this overview summarizes the development of CLASS from the late 1980’s onward.

The basic function of CLASS is to integrate the energy and water balances of the land surface forward in time from an initial starting point, making use of atmospheric forcing data to drive the simulation. When CLASS is run in coupled mode with a global or regional atmospheric model, the required forcing data are passed to it at each time step over each modeled grid cell from the atmospheric driver. CLASS then performs its internal calculations, evaluating a suite of prognostic and diagnostic variables such as albedo and surface radiative and turbulent fluxes, which are in turn passed back to the driver. CLASS can also be run in uncoupled or offline mode, using forcing data derived from field measurements, and the output values of its prognostic and diagnostic variables can then be validated against observations.

CLASS models separately the energy and water balances of the soil, snow, and vegetation canopy (see the diagram below). The basic prognostic variables consist of the temperatures and the liquid and frozen moisture contents of the soil layers; the mass, temperature, density, albedo and liquid water content of the snow pack; the temperature of the vegetation canopy and the mass of intercepted rain and snow present on it; the temperature and depth of ponded water on the soil surface; and an empirical vegetation growth index (which is not used when CTEM is turned on). These variables must be initialized, and a set of physical parameters describing the soil and vegetation existing on the modelled area must be assigned background values, at the beginning of the simulation (see @ref forcingData).

At each time step, CLASS calculates the bulk characteristics of the vegetation canopy on the basis of the vegetation types present over the modelled area. In a pre-processing step, each vegetation type is assigned representative values of parameters such as albedo, roughness length, annual maximum and minimum plant area index, rooting depth and so on (see @ref initProgVar). These values are then aggregated over vegetation categories identified by CLASS (commonly includes: needleleaf trees, broadleaf trees, shrubs, crops, and grass (i.e. short vegetation). The physiological characteristics of the vegetation in each category are determined at the current time step using the aggregated background parameters and assumed annual or diurnal variation functions. These physiological characteristics are then aggregated to produce the bulk canopy characteristics for the current time step.

\image html "schematicDiagramOfClass.png" "Schematic Diagram Of CLASS"
\image latex "schematicDiagramOfClass.png" "Schematic Diagram Of CLASS"

In performing the surface flux calculations the modeled area is divided into up to four subareas: bare soil, vegetation over soil, snow over bare soil, and vegetation over snow. The fluxes of these sub-regions are determined each CLASS timestep, the average value over the modelled area is found and the sub-region fluxes are initialized with the average value at the start of the next time step. The fractional snow coverage is determined using the concept of a threshold snow depth. If the calculated snow depth is less than this value, the snow depth is set to the threshold value and the fractional snow cover is calculated on the basis of conservation of snow mass. The fluxes are calculated for each of the four subareas, and these and the prognostic variables are then areally averaged before being passed back to the atmospheric model.

Originally CLASS performed only one set of these calculations for each grid cell of the model domain. In more recent versions, a “mosaic” option has been added to handle sub-grid scale heterogeneity more effectively. When this option is utilized, each grid cell is divided into a user-specified number of mosaic “tiles”, and the CLASS calculations are performed in turn over each. The surface fluxes are averaged, but the prognostic variables are kept separate for each of the tiles of the mosaic between time steps (see [here for more](@ref compvsmosaic)).

In the CLASSIC offline driver, a gather-scatter operation is included in the driver, mimicking the practice in atmospheric models of “gathering” land surface points on latitude circles onto long vectors prior to the calculations (e.g. src/ctemGatherScatter.f90 or src/classGatherScatter.f90), for improved computational efficiency on vector supercomputers. For CLASS, the mosaic tiles on each of the modelled grid cells are “gathered” onto long arrays prior to calling the CLASS subroutines (thus collapsing the first two dimensions of the arrays into one), and subsequently “scattered” back onto the grid cells before performing the diagnostic averaging calculations. Future code developments will work to make this necessity easier to work with, however developments will follow those adopted in the CanESM framework.

## Development history of CLASS {#devHistory}

\f[
\begin{array}{ | c | c | l | }
1.0 & \text{April 1989} & \text{Basic thermal and hydrological model of snow and soil.} \\
2.0 & \text{August 1991} & \text{Addition of vegetation thermal and hydrological model.} \\
2.1 & \text{May 1993} & \text{Full vectorization of code to enable efficienr running on vector supercomputers.} \\
2.2 & \text{April 1994} & \text{Augmentation of diagnostic calculations; incorporation of in-line comments throughout;} \\
    & & \text{development of a parallel stand-alone version of the model for use with field data.} \\
2.3 & \text{December 1994} & \text{Revisions to diagnostic calculations; new near-surface atmospheric stability functions.} \\
2.4 & \text{August 1995} & \text{Complete set of water budget diagnostic calculations; parametrizations of organic soils} \\
    & & \text{and rock soils; allowance for inhomegeneity between soil layers; incorporation of variable} \\
    & & \text{surface detention capacity.} \\
2.5 & \text{January 1996} & \text{Completion of energy budget diagnostic calculations.} \\
2.6 & \text{August 1997} & \text{Revisions to surface stability function calculations.} \\
2.7 & \text{December 1997} & \text{Incorporation of variable soil permeable depth; calculation of soil thermal and hydraulic} \\
    & & \text{properties based on textural composition; modified surface temperature iteration scheme.} \\
3.0 & \text{December 2002} & \text{Improved treatment of soil evaporation; complete treatment of organic soils; new canopy} \\
    & & \text{conductance formulation; preliminary routines for lateral movement of soil water; enhanced} \\
    & & \text{snow density and snow interception; improved turbulent transfer from vegetation; mosaic} \\
    & & \text{formulation.} \\
3.1 & \text{April 2005} & \text{Faster surface temperature iteration scheme; refinements to leaf boundary resistance} \\
    & & \text{formulation; improved treatment of snow sublimation and interception; transition to} \\
    & & \text{Fortran 90 and single precision variables.} \\
3.2 & \text{May 2006} & \text{Option for multiple soil layers at depth; additional liquid water content of snow pack;} \\
    & & \text{revised radiation transmission in vegetation.} \\
3.3 & \text{December 2006} & \text{Separate temperature profile curve fit for snow and soil; multiple-layer option for ice} \\
    & & \text{sheets; water and energy balance checks for each time step; modifications to soil hydraulic} \\
    & & \text{conductivity calculations.} \\
3.4 & \text{April 2008} & \text{Streamline and clean up code; updated soil thermal conductivity calculations; revisions to} \\
    & & \text{handling of water stored on vegetation.} \\
3.5 & \text{December 2010} & \text{Updated field capacity calculation; revised treatment of water on canopy; reworked} \\
    & & \text{calculation of baseflow.} \\
3.6 & \text{December 2011} & \text{Revised ponding depth over organic soils; revised snow albedo refreshment threshold; new} \\
    & & \text{snow thermal conductivity algorithm; interface with Canadian Terrestrial Ecosystem Model} \\
    & & \text{(CTEM).} \\
3.6.1 & \text{December 2016} & \text{New treatment of bare soil albedo; new optional four-band snow albedo formulation;} \\
    & & \text{fixes to guard against overshoots in water drawdown by evapotranspiration; upper limit on} \\
    & & \text{snow depth.} \\
3.6.2 & \text{July 2019} & \text{CLASS is formerly incorporated into the Canadian Land Surface Scheme including} \\
    & & \text{Biogeochemical Cycles (CLASSIC).} \\
\end{array}
\f]
