%% Step 1
% Clean MAT files
% RUN
% convertOLD_NEW_AOmat.m
% Function
% Inputs: 
% 1. Location for old .mat files
% 2. Location for new .mat files
% Example 
% oldmat = 'I:\01_Coding_Datasets\SLEEP VIDEO\UNMC_12\oldmat';
% newmat = 'I:\01_Coding_Datasets\SLEEP VIDEO\UNMC_12\newmat';
% convertOLD_NEW_AOmat(oldmat, newmat)
% cleanAOstruct
%% Step 2
% Create EEG_LFP TimeTable
% RUN
% SleepDataCreateTT_EEGoLFP_v5.mlapp




%% Step 3 - Consensus setup
% download individual folders from Box
% loadF = 'J:\01_Coding_Datasets\SLEEP_Score_Consensus\3_UNMC_FinalCon\Final'
% saveF = 'J:\01_Coding_Datasets\SLEEP_Score_Consensus\3_UNMC_FinalCon\FinalCC'
folderLOC = 'C:\Users\Admin\Downloads\LOad_UNMC_21';
saveLOC = 'C:\Users\Admin\Downloads\SAve_UNMC_21';
createConsensusTTFirst_newFormat(folderLOC , saveLOC)


%% Step 4 -- Final Score Consensus

