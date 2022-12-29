
% SigNal = 'lfpNF';
caseES = {'UNMC_21','UNMC_20','UNMC_22'};

% RE-PROCESS F221016-0001 - UNMC_21 Night 3, first mat file

for ci = 1:length(caseES)

    batchProcessLFPconvert(caseES{ci},'psg')
    batchProcessLFPconvert(caseES{ci},'lfpNF')

end


%%

