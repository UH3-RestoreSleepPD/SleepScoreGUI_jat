mainLIST = dir('*.mat');
secLIST = {mainLIST.name};

%%
start = 24;

for si = start:length(secLIST)

    tmpName = secLIST{si};
    tmpElems = split(tmpName,'_');
    subID = [tmpElems{1},'_',tmpElems{2}];
    nightID = tmpElems{3};

    fscoresM = ['G:\SleepStudy2_Aim1\finalSCORES\',subID,'_Fscores.mat'];

    load(fscoresM,'finalCS')

    tmpScores = finalCS.(['N',nightID]);

    load(tmpName,'LFPTT');

    lfpttmp = removevars(LFPTT,'ActTime');

    lfpttmp.FSScore = tmpScores;

    LFPTT = lfpttmp;

    % Save in JZ folder 
    saveLOC = ['G:\SleepStudy2_Aim1\JZ_ANN\',tmpName];
    save(saveLOC,"LFPTT")

end