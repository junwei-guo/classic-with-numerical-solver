# Legacy Filename Reference Table {#legacyFileNames}

These tables link CLASS and CTEM filenames prior to the code refactoring to the new filesnames (if changed).

| CLASS Legacy Filename | New Filename  |
|--------------|------------------------|
| APREP.f      | calcLandSurfParams.f90 |
| CANADD.f     | canopyInterception.f90 |
| CANALB.f     | canopyAlbedoTransmiss.f90|
| CANVAP.f     | canopyWaterUpdate.f90  |
| CGROW.f      | classGrowthIndex.f90   |
| CHKWAT.f     | checkWaterBudget.f90   |
| CLASSA.f     | radiationDriver.f90    |
| CLASSB.f     | soilProperties.f90     |
| CLASSG.f     | classGather - classGatherScatter.f90|
| CLASSI.f     | atmosphericVarsCalc.f90|
| CLASSS.f     | classScatter - classGatherScatter.f90|
| CLASST.f     | energyBudgetDriver.f90 |
| CLASSW.f     | waterBudgetDriver.f90  |
| CLASSZ.f     | energyWaterBalanceCheck.f90|
| CWCALC.f     | canopyPhaseChange.f90  |
| GATPREP.f    | classGatherPrep - classGatherScatter.f90|
| GRALB.f      | groundAlbedo.f90       |
| GRDRAN.f     | waterFlowNonInfiltrate.f90|
| GRINFL.f     | waterFlowInfiltrate.f90|
| ICEBAL.f     | iceSheetBalance.f90    |
| PHTSYN3.f    | photosynCanopyConduct.f90|
| SCREENRH.f   | screenRelativeHumidity.f90|
| SNINFL.f     | snowInfiltrateRipen.f90|
| SNOADD.f     | snowAddNew.f90         |
| SNOALBA.f    | snowAlbedoTransmiss.f90|
| SNOALBW.f    | snowAging.f90          |
| SNOVAP.f     | snowSublimation.f90    |
| SUBCAN.f     | waterUnderCanopy.f90   |
| TFREEZ.f     | pondedWaterFreeze.f90  |
| TMELT.f      | snowMelt.f90           |
| TMCALC.f     | waterUpdates.f90       |
| TNPOST.f     | soilHeatFluxCleanup.f90|
| TNPREP.f     | soilHeatFluxPrep.f90   |
| TPREP.f      | energyBudgetPrep.f90   |
| TSOLVC.f     | energBalVegSolve.f90   |
| TSOLVE.f     | energBalNoVegSolve.f90 |
| TSPOST.f     | snowTempUpdate.f90     |
| TSPREP.f     | snowHeatCond.f90       |
| TWCALC.f     | soilWaterPhaseChg.f90  |
| WEND.f       | waterBaseflow.f90      |
| WFILL.f      | waterInfiltrateUnsat.f90|
| WFLOW.f      | waterInfiltrateSat.f90 |
| WPREP.f      | waterCalcPrep.f90      |
| XIT.f        | errorHandler.f90       |



| CTEM Legacy Filename | New Filename  |
|--------------|------------------------|
| allocate.f  | allocateCarbon.f90  |
| balcar.f    | balanceCarbon.f90   |
| bio2str.f   | applyAllometry.f90  |
| ctem.f      | ctemDriver.f90      |
| ctemg1.f    | ctemg1 in ctemGatherScatter.f90  |
| ctems2      | ctems2 in ctemGatherScatter.f90  |
| ctemg2      | ctemg2 in ctemGatherScatter.f90  |
| hetres_mod.f90   | heterotrophicRespiration.f90  |
| mainres.f   | autotrophicRespiration.f90   |
| mortality.f | mortality.f90  |
| phenolgy.f   | phenolgy.f90  |
| soil_ch4uptake | soilCH4Uptake - methaneProcesses.f90|
| wetland_methane.f90 | wetlandMethane - methaneProcesses.f90|
