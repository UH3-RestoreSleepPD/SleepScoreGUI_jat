%% STEPS for ACCEL PROCESSING

% Process RAW data
% use blockTrialExtract_as
close all
rLC = 'C:\Users\Admin\Documents\Github\SleepScoreGUI_jat\RCpS_Accel_OOBed\Thompson_0622\Trial_1';

blockTrialExtract_as(rLC,2,0,1)

%%

close all
rLC2 = 'C:\Users\johna\OneDrive\Documents\Github\SleepScoreGUI_jat\RCpS_Accel_OOBed\Platt_0622';

blockTrialExtract_as(rLC2,2,0,1)

%%

close all
rLC2 = 'C:\Users\Admin\Documents\Github\SleepScoreGUI_jat\RCpS_Accel_OOBed\QuickTESTjat';

cd(rLC2)
matDirall = dir('*.mat');
matDirnames = {matDirall.name};
dataFileName = matDirnames{contains(matDirnames,'Data')};
timeFileName = matDirnames{contains(matDirnames,'Times')};
load(timeFileName,'Acctable');
load(dataFileName,'GeneratedData');
accelTABLE = GeneratedData(:,8:10);
accelTABLEn = accelTABLE(~isnan(table2array(accelTABLE(:,1))),:);

stackedplot(accelTABLEn)

start = transpose([411 , 1083, 1211, 1345 , 1500, 1608, 1702, 1881, 1988, 2058, 2203]);
stop = transpose([1031 , 1199 , 1318 , 1470 , 1599, 1698, 1785, 1965, 2053,2125, 2436]);

allIndex = [start , stop];
fulltrialID = transpose({'M1_T1','M2_T1','M2_T2','M2_T3','M3_T1','M3_T2','M3_T3',...
               'M4_T1','M4_T2','M4_T3','M5_T1'});
MoveID = cellfun(@(x) x(1:2) , fulltrialID,'UniformOutput',false);
trialID = cellfun(@(x) str2double(x(5)), fulltrialID, 'UniformOutput',true);

useTable = table(allIndex(:,1),allIndex(:,2),fulltrialID,MoveID,trialID,...
    'VariableNames',{'StartIndraw','StopTIndraw','FullTrialID','MoveID','TrialID'});

blockTrialExtract_as(rLC2,2,0,1,useTable)


%%

close all
rLC2 = 'C:\Users\johna\OneDrive\Documents\Github\SleepScoreGUI_jat\RCpS_Accel_OOBed\Machnika';

blockTrialExtract_as(rLC2,2,0,1)


%%


close all
rLC2 = 'C:\Users\Admin\Documents\Github\SleepScoreGUI_jat\RCpS_Accel_OOBed\subjectUN25';

cd(rLC2)
matDirall = dir('*.mat');
matDirnames = {matDirall.name};
dataFileName = matDirnames{contains(matDirnames,'Data')};
timeFileName = matDirnames{contains(matDirnames,'Times')};
load(timeFileName,'Acctable');
load(dataFileName,'GeneratedData');
accelTABLE = GeneratedData(:,8:10);
accelTABLEn = accelTABLE(~isnan(table2array(accelTABLE(:,1))),:);

figure;
stackedplot(accelTABLEn)

start = transpose([1426 , 4021, 4571, 4896, 5930, 6429, 6763, 9974, 11588, 13053, 15127]);
stop =  transpose([2799 , 4451, 4842, 5267, 6371, 6734, 7117, 11474, 12851,14256, 16335]);

allIndex = [start , stop];
fulltrialID = transpose({'M1_T1','M2_T1','M2_T2','M2_T3','M3_T1','M3_T2','M3_T3',...
               'M4_T1','M4_T2','M4_T3','M5_T1'});
MoveID = cellfun(@(x) x(1:2) , fulltrialID,'UniformOutput',false);
trialID = cellfun(@(x) str2double(x(5)), fulltrialID, 'UniformOutput',true);

useTable = table(allIndex(:,1),allIndex(:,2),fulltrialID,MoveID,trialID,...
    'VariableNames',{'StartIndraw','StopTIndraw','FullTrialID','MoveID','TrialID'});

blockTrialExtract_as(rLC2,2,0,1,useTable)
