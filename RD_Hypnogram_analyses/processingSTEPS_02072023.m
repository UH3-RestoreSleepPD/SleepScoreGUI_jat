% 1. actigraphyProcess_RD_pre
% 2. processACTpyData_RD
% 3. ?

%%

dirI = 'D:\VRAWLFP\ActigraphyProcessALL\INTRA_ALL';

actigraphyProcess_RD_pre(dirI)

%%

actmd = 'D:\VRAWLFP\ActigraphyProcessALL\INTRA_All_Mat';
pymd = 'D:\VRAWLFP\IntraAct_Proc';
sD = 'D:\VRAWLFP\ActigraphyProcessALL\INTRA_ACT_ALL';

processACTpyData_RD(actmd,pymd,sD)