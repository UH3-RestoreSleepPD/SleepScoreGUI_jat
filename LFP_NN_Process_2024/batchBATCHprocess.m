cd('D:\VRAWLFP\LFP')
getLISTi = dir();
getLISTii = {getLISTi.name};
getLISTiii = getLISTii(~ismember(getLISTii,{'.','..'}));

% 1 and 2 are done
for gi = 8:10
    tmpName = getLISTiii{gi};
    batchProcess_EEG_LFP_BuildTT(tmpName,'psg')
end