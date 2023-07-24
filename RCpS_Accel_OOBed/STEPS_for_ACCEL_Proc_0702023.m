%% STEPS for ACCEL PROCESSING

% Process RAW data
% use blockTrialExtract_as
close all
rLC = 'C:\Users\Admin\Documents\Github\SleepScoreGUI_jat\RCpS_Accel_OOBed\Thompson_0622\Trial_1';

blockTrialExtract_as(rLC,2,0,1)

%%

close all
rLC2 = 'C:\Users\Admin\Documents\Github\SleepScoreGUI_jat\RCpS_Accel_OOBed\Platt_0622';

blockTrialExtract_as(rLC2,2,0,1)