function [] = createConsensusTTFirst_newFormat(folderLOC , saveLOC)
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

cd(folderLOC)

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

    % Loop through nights
    for ni = 1:length(nightLIST)

        tmpNight = nightLIST{ni};
        load(tmpNight,'TT')

        % Get file info - Subject, Institution, Night, Scorer
        fElems = split(tmpNight,{'_','(',')','.'});

        nighTT = replace(fElems{3},' ','');

        switch nighTT
            case '1'
                if sum(matches(TT.Properties.VariableNames,'UNMC'))
                    TT1.(scoreID{1}) = TT.UNMC;
                else
                    TT1.(scoreID{1}) = TT.(scoreID{1});
                end
            case '2'
                if sum(matches(TT.Properties.VariableNames,'UNMC'))
                    TT2.(scoreID{1}) = TT.UNMC;
                else
                    TT2.(scoreID{1}) = TT.(scoreID{1});
                end
            case '3'
                if sum(matches(TT.Properties.VariableNames,'UNMC'))
                    TT3.(scoreID{1}) = TT.UNMC;
                else
                    TT3.(scoreID{1}) = TT.(scoreID{1});
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