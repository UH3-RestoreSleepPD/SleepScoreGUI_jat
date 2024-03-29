function [] = LW_consensusAnalysis_v1(stagE , dirID)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
%
%
% Folder structure
% Consensus Folder
% Separate cases Folder
% ---- CK Folder
%      ---- NIGHT 1 data
%      ---- NIGHT 2 data
%      ---- NIGHT 3 data
% ---- ST Folder
% ---- LW Folder
% ---- MS Folder
% ---- NIGHT 1 data (any case)
% ---- NIGHT 2 data (any case)
% ---- NIGHT 3 data (any case)

mainLOC = [dirID , ':\01_Coding_Datasets\LW_ConsensusStudy']; % NeuroWork
saveDIR = [dirID , ':\01_Coding_Datasets\LW_ConsensusStudy\Summary']; % NeuroWork

switch stagE
    case 1 % extract initial scores
        stageDir = [mainLOC , filesep , 'InitialReview'];
        cd(stageDir)
        allLIST = getFList(1);

        tmpDat = cell(1,length(allLIST));
        for ci = 1:length(allLIST)

            load(allLIST{ci},'TT')
            tmpTab = TT(:,{'LW','ST','MS','CK'});
            tmpDat{ci} = tmpTab;

        end

        initialDat = tmpDat;

        cd(saveDIR)
        save('InitialReview.mat','initialDat','allLIST');

        disp('Initial data sort done!')

        return

    case 2 % extract consensus 1 scores and final scores
        stageDir = [mainLOC , filesep , 'Consensus1'];
        cd(stageDir)
        allLIST = getFList(1);

        tmpDat1 = cell(1,length(allLIST));
        tmpDat2 = cell(1,length(allLIST));
        for ci = 1:length(allLIST)

            load(allLIST{ci},'TT')
            tmpTab1 = TT(:,{'LW','ST','MS','CK'});
            tmpTab2 = TT(:,"FINALSCORE");
            tmpDat1{ci} = tmpTab1;
            tmpDat2{ci} = tmpTab2;

        end

        Con1Dat = tmpDat1;
        ConFDat = tmpDat2;

        cd(saveDIR)
        save('Consensus1.mat','Con1Dat','allLIST');
        save('ConsensusF.mat','ConFDat','allLIST');

        disp('Consensus data sort done!')

        return

    case 3 % Analysis

        % Go to summary folder
        cd(saveDIR)

        % Load in each
        load('InitialReview.mat','initialDat','allLIST')
        initLIST = allLIST;
        load('Consensus1.mat','Con1Dat','allLIST')
        consen1LIST = allLIST;
        load('ConsensusF.mat','ConFDat','allLIST')
        consenFLIST = allLIST;

        for di = 1:3

            switch di
                case 1
                    initialScores = processAlldat(initialDat);
                    initialprocess = epochSummarize(initialScores,1);
                case 2
                    consen1Scores = processAlldat(Con1Dat);
                    consen1process = epochSummarize(consen1Scores,2);
                case 3
                    consenFprocess = epochSummarize(ConFDat,3);

            end
        end

        % For each % Agree per stage

        % Loop through and mash together in new table with DATE name in
        % column
        finTable = table;
        stageCol = {};
        caseCol = {};
        for sti = 1:3

            tmpTable = table;
            switch sti
                case 1

                    for ei1 = 1:length(initLIST)

                        tmpTable = [tmpTable ; initialprocess{ei1}];
                        caseCol = [caseCol ; repmat(initLIST(ei1),6,1)];

                    end

                    stageCol = [stageCol ; repmat({'Initial'},height(tmpTable),1)];

                case 2

                    for ei1 = 1:length(consen1LIST)

                        tmpTable = [tmpTable ; consen1process{ei1}];
                        caseCol = [caseCol ; repmat(consen1LIST(ei1),6,1)];

                    end

                    stageCol = [stageCol ; repmat({'Consensus1'},height(tmpTable),1)];

                case 3

                    for ei1 = 1:length(consenFLIST)

                        tmpTable = [tmpTable ; consenFprocess{ei1}];
                        caseCol = [caseCol ; repmat(consenFLIST(ei1),6,1)];

                    end

                    stageCol = [stageCol ; repmat({'ConsensusF'},height(tmpTable),1)];


            end
            finTable = [finTable ; tmpTable];
        end

        % SAVE file
        finalTABLE = finTable;
        finalTABLE.CaseID = caseCol;
        finalTABLE.StageID = stageCol;
        writetable(finalTABLE, 'ConsensusScoreSummary.xlsx')



end

end % END of function



function [allScoreDat] = processAlldat(allcell)


allScoreDat = cell(size(allcell));
for ii = 1:length(allcell)
    tmpDay = allcell{ii};
    finScore = cell(height(tmpDay),1);
    perAgree = zeros(height(tmpDay),1);
    for ei = 1:height(tmpDay)
        tmpE = table2cell(tmpDay(ei,1:4));
        tabU = tabulate(tmpE);
        if tabU{1,3} >= 75
            finScore{ei} = tabU{1,1};
            perAgree(ei) = tabU{1,3};
        else
            finScore{ei} = 'U';
            perAgree(ei) = tabU{1,3};
        end
    end
    allScoreDat{ii} = table(finScore,perAgree,'VariableNames',{'FinalScore','PercentAgree'});
end

end



function [stageProcess] = epochSummarize(inScores, stageEE)

stageIDS = {'N1','N2','N3','R','W','U'};

stageProcess = cell(size(inScores));
for is = 1:length(inScores)
    tmpDay = inScores{is};

    if ismember(stageEE,[1,2])
        tmpTAB = tabulate(tmpDay.FinalScore);
    else
        tmpDay(cellfun(@(x) isempty(x), tmpDay.FINALSCORE),:) = [];
        tmpTAB = tabulate(tmpDay.FINALSCORE);
    end

    tmpTABo = cell2table(tmpTAB,'VariableNames',{'Stage','Count','Percent'});

    if length(tmpTABo.Stage) < length(stageIDS)

        % Determine which are missing
        st2add = stageIDS(~ismember(stageIDS,tmpTABo.Stage));

        for ssi = 1:length(st2add)
            newRow = height(tmpTABo) + 1;
            tmpTABo.Stage(newRow) = st2add(ssi);
            tmpTABo.Count(newRow) = 0;
            tmpTABo.Percent(newRow) = 0;
        end
    end
    stageProcess{is} = tmpTABo;
end


end















function [fLIST] = getFList(type)

if type == 1
    mdir1 = dir('*.mat');
    fLIST = {mdir1.name};

else
    dir1 = dir();
    dir2 = {dir1.name};
    fLIST = dir2(~ismember(dir2,{'.','..'}) & [dir1.isdir]);
end


end