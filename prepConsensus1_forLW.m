

maindir = 'C:\Users\Admin\Downloads\FinalConsensus-selected';
cd(maindir)
saveDIR = 'C:\Users\Admin\Downloads\FinalConsensus-selected\sumLOC';

matDir = dir('*.mat');
matDir2 = {matDir.name};

for mi = 1:length(matDir2)
    cd(maindir)
    load(matDir2{mi},'finalCS');

    fildNme = fieldnames(finalCS);
    for fi = 1:length(fildNme)

        tmpFN = finalCS.(fildNme{fi});

        fparTS = split(matDir2{mi},'_');

        saveName = [fparTS{1},'_',fparTS{2},'_',num2str(fi),'.mat'];

        TT = table(tmpFN,'VariableNames',{'FINALSCORE'});

        cd(saveDIR)
        save(saveName,"TT");
    end
end