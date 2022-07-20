function [] = finalConsensusSingleScore(folderLOC , saveLOC , subNUM , institut)
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
% Determine which are score folders
for ffi = 1:length(dir3)
    tmpLOC = [folderLOC , filesep , dir3{ffi}];

    switch dir3{ffi}
        case {'CK','ST','LW','MS'}
            % save tmp TT of night / subject / institution with only ID
            % column
            cd(tmpLOC)
            curList = getFList(1);

            fElems = split(curList,{'_','.'});
            if length(curList) == 1
                nightS = fElems(3);
            else
                nightS = fElems(:,:,3);
            end

            for ci = 1:length(curList)

                load(curList{ci},'TT')

                TTn.(['N', nightS{ci}]).(dir3{ffi}) = TT.FINALSCORE;
           
            end

        case 'ALL'
            cd(tmpLOC)
            curList = getFList(1);
            fElems = split(curList,{'_','.'});
            if length(curList) == 1
                nightS = fElems(3);
            else
                nightS = fElems(:,:,3);
            end

            for ci = 1:length(curList)

                load(curList{ci},'TT')

                TTn.(['N', nightS{ci}]).FINAL = TT.FINALSCORE;
           
            end


    end

end



% Determine final from CC
nightStn = fieldnames(TTn);

for ni = 1:length(nightStn)

    tmpTT = TTn.(nightStn{ni});

    subTNf = fieldnames(tmpTT);

    if length(subTNf) > 1

        allComb = struct2table(tmpTT);

        for ccf = 1:height(allComb)

            tROW = table2cell(allComb(ccf,:));
            emCheck = cellfun(@(x) isempty(x), tROW, 'UniformOutput',true);
            ui = tROW(~emCheck);
            uiF = unique(ui);
            
            
            if length(uiF) == 1
                TTn.(nightStn{ni}).FINAL{ccf,1} = uiF{1};
            else
                TTn.(nightStn{ni}).FINAL{ccf,1} = 'U';
            end
        end
    else
        continue
    end
end

% Creat final cell
ttnNights = fieldnames(TTn);

finalCS = struct;
for ti = 1:length(ttnNights)

    finalCS.(ttnNights{ti}) = TTn.(ttnNights{ti}).FINAL;

end


saveID = [subNUM, '_' , institut ,'_Fscores.mat'];

cd(saveLOC)
save(saveID,"finalCS")





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