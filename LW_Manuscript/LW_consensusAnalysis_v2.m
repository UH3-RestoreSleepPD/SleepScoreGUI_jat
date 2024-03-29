function [] = LW_consensusAnalysis_v2(stagE , dirID)
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
saveDIR = [dirID , ':\01_Coding_Datasets\LW_ConsensusStudy\Summary2']; % NeuroWork
saveDIR3 = [dirID , ':\01_Coding_Datasets\LW_ConsensusStudy\Summary3']; % NeuroWork
saveDirFc = [dirID   , ':\01_Coding_Datasets\LW_ConsensusStudy\SummaryF'];

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

        tmpDat2 = cell(1,length(allLIST));
        for ci = 1:length(allLIST)

            load(allLIST{ci},'TT')
            tmpTab2 = TT(:,"FINALSCORE");
            tmpDat2{ci} = tmpTab2;

        end

        ConFDat = tmpDat2;

        cd(saveDIR)
        save('ConsensusF.mat','ConFDat','allLIST');

        disp('Consensus data sort done!')

        return

    case 3 % Analysis

        % Go to summary folder
        cd(saveDIR)

        % Load in each
        load('InitialReview.mat','initialDat','allLIST')
        initLIST = allLIST;
        load('ConsensusF.mat','ConFDat','allLIST')
        consenFLIST = allLIST;

        for di = 1:2

            switch di
                case 1
                    initialScores = processAlldat(initialDat);
                    initialprocess = epochSummarize(initialScores,1);
                case 2
                    consenFprocess = epochSummarize(ConFDat,3);

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
        writetable(finalTABLE, 'ConsensusScoreSummary2.xlsx')

    case 4 % extract initial scores
        stageDir = [mainLOC , filesep , 'ConsensusFF'];
        cd(stageDir)
        allLIST = getFList(1);

        tmpDat = cell(1,length(allLIST));
        for ci = 1:length(allLIST)

            load(allLIST{ci},'TT')
            tmpTab = TT(:,{'LW','ST','MS','CK'});
            tmpDat{ci} = tmpTab;

        end

        conff1 = tmpDat;

%         cd(saveDIR3)
%         save('InitialReview.mat','initialDat','allLIST');
% 
%         disp('Initial data sort done!')

        stageDir2 = [mainLOC , filesep , 'ConFFcreate'];
        cd(stageDir2)
        foldLIST = getFList(2);

        tmpDat2 = cell(1,length(foldLIST)*4);
        nightLIST = cell(1,length(foldLIST)*4);
        nightCount = 1;
        fileCount = 1;
        for fi = 1:length(foldLIST) % 8 subjects

            tmpFld = [stageDir2 , filesep , foldLIST{fi}];
            cd(tmpFld)

            revNs = getFList(2);

            % preload to get cell height
            chkhdir = [tmpFld , filesep , revNs{1}];
            cd(chkhdir)
            chkHei = getFList(1);
            load(chkHei{1})
%             epochSs = cell(height(TT),length(revNs));
            tmpRscores = cell(length(revNs),3);
            for ri = 1:length(revNs) % 2-4 Reviewers
                tmpScF = [tmpFld , filesep , revNs{1}];
                cd(tmpScF)
                tmpNights = getFList(1);
                
                for ni2 = 1:length(tmpNights) % 3 nights
                    tmpSngNiht = tmpNights{ni2};
                    
                    load(tmpSngNiht,'TT')
                    tmpScore = TT.FINALSCORE;
                    tmpRscores{ri,ni2} = tmpScore;
                    if ri == length(revNs)
                        nightLIST{nightCount} = tmpSngNiht;
                        nightCount = nightCount + 1;
                    end
                    disp(['Night ',num2str(ni2), ' done!'])
                end
                disp(['Reviewr ',num2str(ri), ' done!'])
            end
            % Create tables and save files

            for tabi = 1:width(tmpRscores)

                flatCell = transpose(tmpRscores(:,tabi));
                readyTab = unPkc(flatCell);
                outStab = cell2table(readyTab,'VariableNames',revNs);
                tmpDat2{fileCount} = outStab;
                fileCount = fileCount + 1;

            end
            disp(['Subject ',num2str(fi), ' done!'])
        end

        nightLIST = nightLIST(cellfun(@(x) ~isempty(x), nightLIST));
        conff2 = tmpDat2(cellfun(@(x) ~isempty(x), tmpDat2));
        finalConfL = [conff1 , conff2];
        finalNlist = [allLIST , nightLIST];

        % Save location
        cd(saveDirFc)
        % Save name
        saveName = 'finalReview.mat';
        % Save
        save(saveName,"finalConfL","finalNlist");



        return



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



function [unpacked] = unPkc(inpCell)


unpacked = cell(length(inpCell{1}),length(inpCell));
for ci = 1:length(inpCell)

    unpacked(:,ci) = inpCell{ci};

end


end