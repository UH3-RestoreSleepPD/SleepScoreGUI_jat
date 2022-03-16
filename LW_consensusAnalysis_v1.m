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

        % Load in each

        % Fraction of epochs with 50%

        % Fraction of epochs with 75% 

        % Fraction of epochs with 100%

        % ID of each epoch 



end







% Get folder dir
dir3 = getFList(2);

% Load all Nights and create final save file
finTTlist = getFList(1);
for fi = 1:length(finTTlist)

    fElems = split(finTTlist{fi},{'_','(',')'});

    subID = replace(fElems{1},' ','');
    instit = replace(fElems{2},' ','');

    % Load
    load(finTTlist{fi},'TT')

    % Clear sleep scores
    TT = removevars(TT,{'STNF', 'UNMC'});

    % Add var columns
    LW = TT.notes;
    ST = TT.notes;
    CK = TT.notes;
    MS = TT.notes;
    TT = addvars(TT,LW,ST,CK,MS,'After','ChinZ');

    % Rename
    switch fi
        case 1
            TT1 = TT;
        case 2
            TT2 = TT;
        case 3
            TT3 = TT;
    end
end

% Loop through unique score data
for di = 1:length(dir3)
    % Temp Score Dir
    cd([folderLOC, filesep , dir3{di}])
    % Get list of Nights
    nightLIST = getFList(1);

    % Get Scorer ID
    scoreID = split(dir3{di},'_');
%     scoreID = replace(foldINfo{2},' ','');
    if ismember(scoreID,{'LW','ST','CK'})
        cOLs = 'STNF';
    else
        cOLs = 'UNMC';
    end

    % Loop through nights
    for ni = 1:length(nightLIST)

        tmpNight = nightLIST{ni};
        load(tmpNight,'TT')

        % Get file info - Subject, Institution, Night, Scorer
        fElems = split(tmpNight,{'_','(',')','.'});

        nighTT = replace(fElems{3},' ','');

        switch nighTT
            case '1'
                if ismember(scoreID{1},TT.Properties.VariableNames)
                    TT1.(scoreID{1}) = TT.(scoreID{1});
                else
                    TT1.(scoreID{1}) = TT.(cOLs);
                end

            case '2'
                if ismember(scoreID{1},TT.Properties.VariableNames)
                    TT2.(scoreID{1}) = TT.(scoreID{1});
                else
                    TT2.(scoreID{1}) = TT.(cOLs);
                end

            case '3'
                if ismember(scoreID{1},TT.Properties.VariableNames)
                    TT3.(scoreID{1}) = TT.(scoreID{1});
                else
                    TT3.(scoreID{1}) = TT.(cOLs);
                end
        end

    end

end

cd(saveLOC)
for sti = 1:3
    switch sti
        case 1
            TT = TT1;
            saveN = [subID,'_',instit,'_1_CS.mat'];
        case 2
            TT = TT2;
            saveN = [subID,'_',instit,'_2_CS.mat'];
        case 3
            TT = TT3;
            saveN = [subID,'_',instit,'_3_CS.mat'];
    end
    save(saveN,"TT");
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