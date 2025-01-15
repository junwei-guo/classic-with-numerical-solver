# Running CLASSIC {#runCLASSIC}

1. @ref runStandAloneMode "At a point"
2. @ref runGrid "Over a grid"


---


# Running CLASSIC for a point location {#runStandAloneMode}

Once all input files are prepared CLASSIC can be run over a point location by ([first make sure you compiled using the correct option!](@ref compilingMod)):

        bin/CLASSIC_serial pathToJobOptionsFile longitude/latitude

E.g.

        bin/CLASSIC_serial configurationFiles/template_job_options_file.txt 105.23/40.91

or, we can use the shorthand *0/0* to specify the whole domain since there is only one gridcell:

        bin/CLASSIC_serial configurationFiles/template_job_options_file.txt 0/0

## Special note about fire and mortality for runs at a point location

The fire parameterization for CLASSIC should be turned off for runs at a point location. At a single location fire is a stochastic process that cannot be simulated reasonably by the model.

Just like fire, mortality is also a spatial process. For example, in CTEM approximately 1-2% of trees are killed to account for age related mortality. Clearly, trees in a small plot do not die 1% every year, but on a landscape scale on average 1% of trees may die. So when comparing model simulated biomass to point scale observations it may be desirable to switch off mortality. To turn off mortality you can insert 'return' after loop 140 in mortality.f90


# Running CLASSIC over a grid (Global or Regional) {#runGrid}

Once all input files are prepared CLASSIC can be run over a grid by ([make sure you compiled using the correct option!](@ref compilingMod)):

        mpirun -np # bin/CLASSIC_parallel pathToJobOptionsFile Wlongitude/Elongitude/Slatitude/Nlatitude

Where the # is replaced by the number of cores you wish to run the model on. E.g.

        mpirun -np 2 bin/CLASSIC_parallel configurationFiles/template_job_options_file.txt 90.5/105.5/30.4/45.5

or, we can use the shorthand *0/0/0/0* to specify the whole domain if you want to just run your whole region:

        mpirun -np 2 bin/CLASSIC_parallel configurationFiles/template_job_options_file.txt 0/0/0/0
