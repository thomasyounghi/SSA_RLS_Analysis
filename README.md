# SSA_RLS_Analysis
Analysis scripts for measuring the effect of replicative aging on single strand annealing in yeast cells.

## Order to run analyses scripts
1_CombiningData.R

2_CleaningDataforFl.R
2_IdentifyingSSACellsToCircle.R

3_CombiningManualBG_autocellfl.R
3_CombiningManualBG_circledcellfl.R

4_MergingCircledAndAutoCellFl_bgcorrected.R

5_brightneighboringcells.R
6_adjustmentforbrightneighbors.R

7_ClassifyingYFP.R  (the classifications here were not ultimately used)
7_ClassifyingYFP_alltimes.R
7_ControlStrain_RFPandYFP.R

9_AssessingRepair_5hafterdoxremoval.R
9_AssessingRepairatLaterTimePoints.R
9_PlottingAgeatDox.R

10_AssessingmCherryLoss_mCherryDegron_AbsoluteCutoff.R
10_ComparingYFPAppearanceTimes.R
10_PlottingTrajectories_YFPoffatdox_Alive5hafterDoxRemoval_allstrains_cowplot.R

11_Plotting_repairfrac_byexp.R
11_ComparingYoungandOld_ContingencyTables.R
11_Ratio_oldrepair_over_youngrepair.R

11_BTin5hafterdoxremoval_roundingtonearestbudbefore.R
11_BT_Comparison_PrePostDox_notbyyfpclass.R




#YFP classification
7_Revised_YFPrepair_assessment_cellstomanuallycheck
7_Revised_YFPrepair_assessment.R
<-these two files were combined.

#This was originally used,but the inputs changed. Remove this
8_CombiningManualAutoYFPclass

#This was ultimately used
8_CombiningManualAutoYFPclass_allstrains.R



