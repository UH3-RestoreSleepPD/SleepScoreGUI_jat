function [] = convertOLD_NEW_AOmat(oldLOC , newLOC)

% cd('I:\01_Coding_Datasets\SLEEP VIDEO\UNMC_11\oldMat')
% 
% raw = 'I:\01_Coding_Datasets\SLEEP VIDEO\UNMC_11\oldMat';
% new = 'I:\01_Coding_Datasets\SLEEP VIDEO\UNMC_11\newMat';

cd(oldLOC)

raw = oldLOC;
new = newLOC;

mdir = dir('*.mat');
mdir2 = {mdir.name};
fileLIST = mdir2;


for fi = 1:length(fileLIST)
    cd(raw)
    matob = matfile(fileLIST{fi});
    varlist = who(matob);
    chanlist = varlist(~contains(varlist,{'SF','Ports','CADD','Channel','CANALOG','new','raw','fileLIST'}));
    
    cd(raw)
    load(fileLIST{fi},chanlist{:});
    cd(new)
    save(fileLIST{fi},chanlist{:});
    
    
end

disp('All AO mat files have been processed')
