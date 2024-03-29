cd('D:\VRAWLFP\LFP')
getLISTi = dir();
getLISTii = {getLISTi.name};
getLISTiii = getLISTii(~ismember(getLISTii,{'.','..'}));

% 1 and 2 are done
for gi = 1:20
    tmpName = getLISTiii{gi};
    disp(tmpName)


    
    if matches(tmpName,{'UNMC_1','UNMC_26'})
        batchProcess_EEG_LFP_BuildTT(tmpName,'lfpNF')
    else
        continue

    end
end