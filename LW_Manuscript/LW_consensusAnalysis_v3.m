function [] = LW_consensusAnalysis_v3(stagE , dirID)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%

mainLOC = [dirID , ':\01_Coding_Datasets\LW_ConsensusStudy']; % NeuroWork
saveDIR = [mainLOC, '\SummaryAnalysis_12_2022\AllData']; % NeuroWork
saveDirFc = [dirID   , ':\01_Coding_Datasets\LW_ConsensusStudy\SummaryF'];

switch stagE
    
    case 1 % Analysis

        % Go to summary folder
        cd(saveDIR)

        % Load in each
        load('InitialReview.mat','initialDat','allLIST')
        initLIST = allLIST;
        load('finalReview.mat','finalConfL','finalNlist')
        finalLIST = finalNlist;
        finalDat = finalConfL;

        for di = 1:2

            switch di
                case 1
                    initialScores = processAlldat(initialDat);
                    initialprocess = epochSummarize(initialScores,1);
                case 2
                    finalScores = processAlldat(finalDat);
                    finalprocess = epochSummarize(finalScores,1);
            end
        end

        % For each % Agree per stage

        % Loop through and mash together in new table with DATE name in
        % column
        finTable = table;
        stageCol = {};
        caseCol = {};
        for sti = 1:2

            tmpTable = table;
            switch sti
                case 1

                    for ei1 = 1:length(initLIST)

                        tmpTable = [tmpTable ; initialprocess{ei1}];
                        caseCol = [caseCol ; repmat(initLIST(ei1),6,1)];

                    end

                    stageCol = [stageCol ; repmat({'Initial'},height(tmpTable),1)];


                case 2

                    for ei1 = 1:length(finalLIST)

                        tmpTable = [tmpTable ; finalprocess{ei1}];
                        caseCol = [caseCol ; repmat(finalLIST(ei1),6,1)];

                    end

                    stageCol = [stageCol ; repmat({'Final'},height(tmpTable),1)];


            end
            finTable = [finTable ; tmpTable];
        end

        % SAVE file
        finalTABLE = finTable;
        finalTABLE.CaseID = caseCol;
        finalTABLE.StageID = stageCol;
        writetable(finalTABLE, 'ConsensusScoreSummary2.xlsx')
        writetable(finalTABLE, 'ConsensusScoreSummary2.csv')
        save('finalSummary.mat','finalTABLE')


  


end

end % END of function



function [allScoreDat] = processAlldat(allcell)


allScoreDat = cell(size(allcell));
for ii = 1:length(allcell)
    tmpDay = allcell{ii};
    finScore = cell(height(tmpDay),1);
    perAgree = zeros(height(tmpDay),1);
    for ei = 1:height(tmpDay)
        tmpE = table2cell(tmpDay(ei,1:width(tmpDay)));

        tmpE(cellfun(@(x) isempty(x), tmpE)) = {'U'};
%             finScore{ei} = nan;
%             perAgree(ei) = nan;
%             continue
%         end
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
            tmpTABo(newRow,:) = [st2add(ssi), 0, 0];
        end
    end

    if sum(~ismember(tmpTABo.Stage,stageIDS)) ~= 0
        missStage = stageIDS(~ismember(stageIDS,tmpTABo.Stage));
        tmpTABo.Stage(~ismember(tmpTABo.Stage,stageIDS)) = missStage;
    end


    stageProcess{is} = tmpTABo;
end


end


