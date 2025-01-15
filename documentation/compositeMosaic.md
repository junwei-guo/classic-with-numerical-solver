# Composite vs. mosaic representation {#compvsmosaic}

CLASSIC can be run in two different configurations: composite or mosaic (see figure below)

In the composite mode, the structural vegetation attributes (including leaf area index, vegetation height, rooting depth) of PFTs that exist in a grid cell are averaged in proportion to their fractional coverages and then used in the grid-averaged energy and water balance calculations. As a result the entire grid cell is characterized by land surface physical environment (including soil temperature, soil moisture, fractional snow cover, and net radiation) that is common to all PFTs.

In contrast, in the mosaic mode, a grid box is split into multiple tiles representing individual PFTs (or land unit, such as soil texture e.g. Melton et al. (2017) \cite Melton2017-gp) for each of which energy, water and carbon balance calculations are performed separately. In principle, the mosaic mode may be used to represent tiles that are characterized by any chosen distinction such a lowlands vs. uplands, soil texture, vegetation, soil depth, etc.

As a result of the differences between the mosaic and composite configurations, the simulated carbon balance evolves somewhat differently in the two configurations despite being driven with identical climate forcing (see Melton and Arora (2014) \cite Melton2014-xk). For previous publications using this model capability see: Li and Arora (2012) \cite Li2012-f7f, Melton and Arora (2014) \cite Melton2014-xk, Shrestha et al. (2016) \cite Shrestha2016-do, and Melton et al. (2017) \cite Melton2017-gp.

Although pretty clever and powerful, running the mosaic version and interpreting the model results can be a logical nightmare, especially with competition on where the fractions of different tiles/mosaic change with time. So if you are new at this, please consider running the model in the composite mode.


Some initialization variables that impact upon composite vs. mosaic model runs:

- **FAREROT** Fractional coverage of mosaic tile on the modelled area
  - If you run with multiple tiles you need to specify what fraction of the grid cell the tile occupies.
- **MIDROT** Mosaic tile type identifier (1 for land surface, 0 for inland lake)
  - CLASSIC runs in the coupled models with a sub-grid lakes scheme (The Canadian Small Lake Model (CSLM); Verseghy and MacKay, 2017 \cite Verseghy2017-ys). The lakes can then be represented as a tile.


\image html "compVsMosaic_MeltonArora_BG_2014.png" "Schematic representation of the composite and mosaic approaches for the coupling of CLASS and CTEM models in a stand-alone mode. (From Melton and Arora, 2014)"
\image latex "compVsMosaic_MeltonArora_BG_2014.png" "Schematic representation of the composite and mosaic approaches for the coupling of CLASS and CTEM models in a stand-alone mode. (From Melton and Arora, 2014)"
