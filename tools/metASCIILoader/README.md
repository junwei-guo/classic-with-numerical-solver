ASCII to NetCDF met file loader {#asciiMet}
========

The metASCIILoader program takes in ASCII met files and converts them to NetCDF format. It is located in tools/metASCIILoader.

# Compile

To compile, use the `make` command from the console.
This will generate the binary executable in the *bin* folder.

# Execute

To run the program, use the following command:

`bin/metASCIILoader [file.MET/csv] [longitude] [latitude] [containsLeaps]`

Where *longitude* and *latitude* are the desired real values for the local file and the *file.MET or file.csv* is your local met file. *containsLeaps* is a logical switch. True if your input file contains leap years, i.e. years that contain Feb 29th, false if not. It is expecting either 'true' or 'false'.

The provided longitude and latitude values will be written as a property in the netcdf file.

See @ref forcingData for expected units and variables.

# Structure

The original, legacy met file (*.MET) looks like this:

        Minutes, Hour, Day, Year, Shortwave, Longwave, Precipitation, Temperature, Humidity, Wind, Pressure
        [...]
          6  0    1  1901     0.00   264.07    0.0000E+00     2.78   4.286E-03    3.02    99308.05
          6 30    1  1901     0.00   259.21    0.0000E+00     4.06   4.607E-03    1.33    99228.51
          7  0    1  1901     0.00   259.21    0.0000E+00     4.43   4.724E-03    1.29    99230.83
          7 30    1  1901     0.00   259.21    0.0000E+00     4.79   4.841E-03    1.26    99233.15
          8  0    1  1901    52.61   259.21    0.0000E+00     5.16   4.958E-03    1.23    99235.47
          8 30    1  1901   106.40   259.21    0.0000E+00     5.53   5.075E-03    1.19    99237.79
          9  0    1  1901   148.44   259.21    0.0000E+00     5.90   5.192E-03    1.16    99240.11
          9 30    1  1901   177.86   259.21    0.0000E+00     6.27   5.309E-03    1.13    99242.43
         10  0    1  1901   195.90   259.21    0.0000E+00     6.64   5.426E-03    1.09    99244.75

         [...]

It is important to ensure there are no extra columns or trailing white spaces after the pressure field. Also if your filename contains *MET* anywhere in it (case sensitive) the script assumes you have a legacy format file. Otherwise it assumes CSV.

You can also give a CSV file like:

          00,30,1,2010,0.0,215.456,1.388888888888889e-05, 24.787,0.0003890691366542001,1.9680000000000002,99497.00000000001
          01,00,1,2010,0.0,215.456,1.388888888888889e-05,-24.663,0.0003958684331327244,1.935,99500.0
          01,30,1,2010,0.0,215.456,0.0,-24.539,0.00040272515229766096,1.902,99503.99999999999
          02,00,1,2010,0.0,215.456,0.0,-24.415,0.0004089842829792506,1.869,99507.0
          02,30,1,2010,0.0,215.456,0.0,-24.291,0.00041596956902755834,1.8359999999999999,99510.0
          03,00,1,2010,0.0,215.456,0.0,-24.166999999999998,0.0004230137137325519,1.804,99514.0
          03,30,1,2010,0.0,199.73,0.0,-24.043000000000003,0.00042946217462082526,1.771,99517.00000000001
          04,00,1,2010,0.0,199.73,0.0,-24.104,0.00042630990863214304,1.859,99509.0
          04,30,1,2010,0.0,199.73,0.0,-24.166,0.000423116599139658,1.9480000000000002,99502.00000000001
          05,00,1,2010,0.0,199.73,0.0,-24.226999999999997,0.0004193268030774631,2.036,99495.0
