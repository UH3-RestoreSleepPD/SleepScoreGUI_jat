function [] = createConsensusTT_FINAL_modified_v2(folderLOC , saveLOC)
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

% Create sub function to figure out the number of NIGHTS
nightDIR = [folderLOC , filesep , dir3{1}];
cd(nightDIR)
[fLIST] = getFList(1);

cd(folderLOC)

% CORRECTLY IDENTIFY NIGHT

% Load all Nights and create final save file
% Loop through nights
for fi = 1:length(fLIST)

    % Loop through folders Scores
    for ffi = 1:length(dir3)

        tmpDir = [folderLOC , filesep , dir3{ffi}];
        cd(tmpDir)
        curList = getFList(1);

        fElems = split(curList,{'_','.'});
        if isscalar(curList)
            nightS = fElems(3);
        else
            nightS = fElems(:,:,3);
        end
        nightSn = cellfun(@(x) str2double(x), nightS,'UniformOutput',true);

        curMAt = curList(fi);
        curMAtels = strsplit(curMAt{1},{'_', '.'});

        subID = replace(curMAtels{1},' ','');
        instit = replace(curMAtels{2},' ','');
        nighInd = str2double(curMAtels{3});

        % Load
        load(curMAt{1},'TT')

        if ffi == 1
            fTT = TT;
            fTT.FINALSCORE = zeros(height(TT),1);
        end

        % Clean up TT
        if matches('FINALSCORE',TT.Properties.VariableNames)
            emLOC = cellfun(@(x) isempty(x), TT.FINALSCORE, 'UniformOutput',true);
        else
            emLOC = 0;
        end

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

    switch nighInd
        case 1
            TT1 = fTT;
        case 2
            TT2 = fTT;
        case 3
            TT3 = fTT;
    end
end



cd(saveLOC)
for sti = 1:length(nightSn)
    tmpNIGHT = nightSn(sti);
    switch tmpNIGHT
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