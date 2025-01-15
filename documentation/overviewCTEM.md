# Overview of the Canadian Terrestrial Ecosystem Model (CTEM) {#overviewCTEM}

Version 1 of the CTEM is the terrestrial carbon cycle component of the second generation Canadian Earth System Model (CanESM2) (Arora et al., 2011)\cite Arora2011-79f where it was coupled to version 2.7 of the Canadian Land Surface Scheme (CLASS). CTEM v. 2.0 (Melton and Arora, 2016) \cite Melton2016-zx has been coupled to CLASS v. 3.6 (Verseghy, 2012) \cite Verseghy2012-c0e. Together CLASS and CTEM form the CLASSIC model which is capable of being run online in the CanESM model or offline, driven by observation-based meteorological forcings. CLASSIC models terrestrial ecosystem processes by tracking the flow of carbon through three living vegetation components (leaves, stem and roots) and two dead carbon pools (litter and soil).

<!-- \f[
\begin{table}[]
\caption{CTEM and peatland PFTs and their mapping to the CLASS PFTs}
\label{my-label}
\begin{array}{|l|l|l|l|l|l|}
\hline
\multicolumn{1}{|c|}{CLASS PFTs} & \multicolumn{3}{c|}{CTEM PFTs}                    & \multicolumn{2}{c|}{Peatland PFTs}  \\ \hline
Needleleaf tree                  & Evergreen & Deciduous      &                      &                  &                  \\ \hline
Broadleaf tree                   & Evergreen & Cold Deciduous & Drought/Dry Decidous & Evergreen Shrubs & Deciduous Shrubs \\ \hline
Crop                             & C$_3$     & C$_4$          &                      &                  &                  \\ \hline
Grass                            & C$_3$     & C$_4$          &                      & Sedges           &                  \\ \hline
\end{array}
\end{table}
\f] -->

The amount of carbon in these five carbon pools is simulated prognostically (see below). In the CLASSIC framework, CLASS uses structural vegetation attributes (including LAI, vegetation height, canopy mass and rooting depth) simulated by CTEM, and CTEM uses soil moisture, soil temperature and net radiation calculated by CLASS. Combined, CLASS and CTEM simulate the atmosphere--land fluxes of energy, water and \f$CO_2\f$.

Version 1.0 of CTEM is described in a collection of papers detailing parametrization of photosynthesis, autotrophic and heterotrophic respiration (Arora, 2003) \cite Arora2003-3b7; phenology, carbon allocation, biomass turnover and conversion of biomass to structural attributes (Arora and Boer, 2005) \cite Arora2005-6b1; dynamic root distribution (Arora and Boer, 2003) \cite Arora2003838; and disturbance (fire) (Arora and Boer, 2005) \cite Arora20052ac. These processes are modelled over prescribed fractional coverage of (at the time, typically) nine PFTs (Wang et al., 2006) \cite Wang2006-he and determine the structural vegetation dynamics including vegetation biomass, LAI, vegetation height, fraction of roots in each of the three soil layers, leaf onset and offset times and primary \f$CO_2\f$ fluxes of gross primary productivity (GPP) and NPP.

CTEM v. 2.0 is described in Melton and Arora, 2016 \cite Melton2016-zx. CTEM v. 2.0 can be run in two different modes, either (i) using specified fractional coverage of its PFTs, or (ii) allowing the fractional coverage of its non-crop PFTs to be dynamically determined based on competition between PFTs. The parametrization for simulating competition between PFTs is found in @ref competition_mod.f90 and described in Arora and Boer (2006a,b), \cite Arora2006-ax \cite Arora2006-pp,  Melton and Arora, 2016 \cite Melton2016-zx, and Shrestha et al. (2016) \cite Shrestha2016-do. The fire parametrization has also been refined in the new model version as described in @ref disturb.f90 and Melton and Arora, 2016 \cite Melton2016-zx.

# Plant Functional Types (PFTs) in CLASSIC {#PFTsCLASSIC}

The original PFT scheme of CLASSIC is as follows (parameters for this configuration are found in the configurationFiles/template_run_parameters.txt file):

| CLASS PFTs | CTEM PFTs --- | ---| ---|
|:---------------:|:---------:|----------------|------------|
| Needleleaf tree (NdlTr) | Evergreen ('NdlEvgTr') | Deciduous ('NdlDcdTr') |  |
| Broadleaf tree (BdlTr) | Evergreen ('BdlEvgTr') | Cold Deciduous ('BdlDCoTr') | Drought/Dry Decidous ('BdlDDrTr') |
| Crop (Crops) | \f$C_3\f$ ('CropC3  ') | \f$C_4\f$ ('CropC4  ')|  |
| Grass (Grass) | \f$C_3\f$ ('GrassC3 ') | \f$C_4\f$ ('GrassC4 ')|  |


The peatland module of CLASSIC (Wu et al. 2016) \cite Wu2016-zt introduced three peatland-specific PFTs. These new PFTs impacted upon the biogeochemistry but are treated like existing CLASS PFTs for the physics calculations (parameters for this configuration are found in the configurationFiles/template_run_parameters_peatlands.txt file):

| CLASS PFTs | CTEM PFTs --- | ---| ---| Peatland PFTs | ---|
|:---------------:|:---------:|----------------|----------------------|:----------------:|------------------|
| Needleleaf tree (NdlTr) | Evergreen ('NdlEvgTr') | Deciduous ('NdlDcdTr') |  |  |  |
| Broadleaf tree (BdlTr) | Evergreen ('BdlEvgTr') | Cold Deciduous ('BdlDCoTr') | Drought/Dry Decidous ('BdlDDrTr') | Evergreen Shrubs ('BdlEvgSh') | Deciduous Shrubs ('BdlDCoSh')|
| Crop | \f$C_3\f$ | \f$C_4\f$ |  |  |  |
| Grass | \f$C_3\f$ | \f$C_4\f$ |  Sedges ('Sedge   ') |  | |

Shrubs have also been implemented as a fifth CLASS PFT allowing for the model physics to distinguish shrubs from trees. Parameters for this configuration are found in the configurationFiles/template_run_parameters_shrubs.txt file):

| CLASS PFTs | CTEM PFTs --- | ---| ---|
|:---------------:|:---------:|----------------|------------|
| Needleleaf tree (NdlTr) | Evergreen ('NdlEvgTr') | Deciduous ('NdlDcdTr') |  |
| Broadleaf tree (BdlTr) | Evergreen ('BdlEvgTr') | Cold Deciduous ('BdlDCoTr') | Drought/Dry Decidous ('BdlDDrTr') |
| Crop (Crops) | \f$C_3\f$ ('CropC3  ') | \f$C_4\f$ ('CropC4  ')|  |
| Grass (Grass) | \f$C_3\f$ ('GrassC3 ') | \f$C_4\f$ ('GrassC4 ')| Sedges ('Sedge   ')  |
| Broadleaf shrub (BdlSh) | Evergreen Shrubs ('BdlEvgSh') | Deciduous Shrubs ('BdlDCoSh')|  |

The model is not limited to these PFTs provided PFT-specific parameters can be determined. As well it is important to thoughtfully integrate new PFTs. Within the model code there are checks for unknown PFTs based upon the short names of the known PFTs in the tables above.

# Rate change equations for carbon pools {#CTEMRateChgEqns}


From the gross canopy photosynthesis rate (\f$G_{\text{canopy}}\f$, @ref PHTSYN3.f), maintenance and growth respirations (\f$R_\mathrm{m}\f$ and \f$R_\mathrm{g}\f$, @ref mainres.f), and
heterotrophic respiration components (\f$R_{\text{h,H}}\f$ and \f$R_{\text{h,D}}\f$, @ref hetres_mod.f90), it is possible to estimate the change in carbon amount of the model's five pools.

When the daily NPP (\f$ G_{canopy} - R_\mathrm{m} - R_\mathrm{g}\f$) is positive, carbon is allocated to the plant's live carbon pools and the rate of change is given by

  \f[
  \frac{\mathrm{d}C_i}{\mathrm{d}t} = a_{fi} \left(G_{canopy}-R_\mathrm{m}-R_\mathrm{g} \right) - D_i - H_i - M_i \\ \quad i = {L, S, R} \qquad (Eqn 1)
   \f]
   <!-- {#rate_change_eqns_live_pools} -->

where \f$a_{fi}\f$ is the corresponding allocation fractions for each pool (stem, root and leaves) and \f$D_i\f$ is the litter produced from these components as explained in @ref phenolgy.f90. \f$H_i\f$ is the loss associated with fire that releases \f$CO_2\f$ and other trace gases to the atmosphere and \f$M_i\f$ is the mortality associated with fire that contributes to the litter pool as explained in @ref disturb.f90.

If the daily NPP is negative (\f$G_{canopy} < R_\mathrm{m}\f$, \f$R_\mathrm{g} = 0\f$), the rate of change is given by

\f[
 \frac{\mathrm{d}C_i}{\mathrm{d}t} = a_{fi}G_{canopy} - R_{m,i}  - D_i  - H_i - M_i, \\ \quad i = {L, S, R} \qquad (Eqn 2)
 \f]
<!-- \label{rate_change_eqns_live_pools2} -->

Negative NPP causes the plant to lose carbon from its live carbon pools due to respiratory costs in addition to the losses due to litter production (\f$D_i\f$) and disturbance (\f$H_i\f$, \f$M_i\f$).

The rate change equations for the litter and soil carbon pools are given by
\f{eqnarray*}{
\frac{\mathrm{d}C_\mathrm{D}}{\mathrm{d}t} &=& D_\mathrm{L} + D_\mathrm{S} +
D_\mathrm{R} + M_\mathrm{L} + M_\mathrm{R} + M_\mathrm{S} - H_\mathrm{D} -C_{\mathrm{D} \rightarrow \mathrm{H}} - R_{h,D} \\
\frac{\mathrm{d}C_\mathrm{H}}{\mathrm{d}t} &=& C_{\mathrm{D} \rightarrow
\mathrm{H}} - R_{h,H}
\qquad (Eqn 3)
\f}
<!-- \label{rate_change_eqns_dead_pools}, -->

where \f$C_{\mathrm{D} \rightarrow \mathrm{H}}\f$ represents the transfer of
humified litter to the soil carbon pool and \f$H_\mathrm{D}\f$
is loss associated with burning of litter associated with fire that releases
\f$CO_2\f$ and other trace gases to the atmosphere.
