function [] = createConsensusTT_FINAL(folderLOC , saveLOC)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% loadF = 'J:\01_Coding_Datasets\SLEEP_Score_Consensus\3_UNMC_FinalCon\Final'
% saveF = 'J:\01_Coding_Datasets\SLEEP_Score_Consensus\3_UNMC_FinalCon\FinalCC'
%
%
% Folder structure
% Consensus Folder
% Separate cases Folder
% CK
% ----MAT
% ----MAT
% ----MAT
% LW
% ----MAT
% ----MAT
% ----MAT
% MS etc

cd(folderLOC)

% Get folder dir
dir3 = getFList(2);

% Load all Nights and create final save file
% Loop through nights
for fi = 1:3

    % Loop through folders Scores
    for ffi = 1:length(dir3)

        tmpDir = [folderLOC , filesep , dir3{ffi}];
        cd(tmpDir)
        curList = getFList(1);

        fElems = split(curList,{'_'});
        nightS = fElems(:,:,3);

        nighInd = matches(nightS,num2str(fi));
        curMAt = curList(nighInd);
        curMAtels = strsplit(curMAt{1},'_');

        subID = replace(curMAtels{1},' ','');
        instit = replace(curMAtels{2},' ','');

        % Load
        load(curMAt{1},'TT')

        if ffi == 1
            fTT = TT;
            fTT.FINALSCORE = [];
        end

        % Clean up TT
        emLOC = cellfun(@(x) isempty(x), TT.FINALSCORE, 'UniformOutput',true);

        if any(emLOC)
            fixX = 1;
        else
            fixX = 0;
        end

        switch dir3{ffi}
            case 'CK'
                if fixX
                    TT.FINALSCORE(emLOC) = TT.CK(emLOC);
                end
                fTT.CK = TT.FINALSCORE;
            case 'LW'
                if fixX
                    TT.FINALSCORE(emLOC) = TT.LW(emLOC);
                end
                fTT.LW = TT.FINALSCORE;
            case 'ST'
                if fixX
                    TT.FINALSCORE(emLOC) = TT.ST(emLOC);
                end
                fTT.ST = TT.FINALSCORE;
            case 'MS'
                if fixX
                    TT.FINALSCORE(emLOC) = TT.MS(emLOC);
                end
                fTT.MS = TT.FINALSCORE;
        end


    end

    switch fi
        case 1
            TT1 = fTT;
        case 2
            TT2 = fTT;
        case 3
            TT3 = fTT;
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