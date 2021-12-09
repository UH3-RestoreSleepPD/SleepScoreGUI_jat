function [] = transfer_Score_NN(oldlfpLoc,transLoc)

cd(oldlfpLoc)

% Get list of case folders
foldList = dir();
foldList2 = {foldList.name};
foldList3 = foldList2(~ismember(foldList2,{'..','.'}));

for fi = 1:length(foldList3)

    tmpFold = [oldlfpLoc , filesep , foldList3{fi}];
    cd(tmpFold);

    % Get list of mats
    matD1 = dir('*.mat');
    matD2 = {matD1.name};
    
    % Get sub ID
    subID = matD2{1}(1);

    % Get List of non-raw
    NNcheck = contains(matD2,'NN');
    if any(NNcheck)
        continue
    end

    matD3 = matD2(~contains(matD2,'raw'));

    % Loop through
    for di = 1:length(matD3)
        cd(tmpFold);
        % Get night
        tmpNight = extractBetween(matD3(di),'UNMC_','_LFP');

        % Load data and rename
        load(matD3{di},'LFPTT');

        % Load Transfer Loc + SubID + _UNMC_ + Night
        transDat = [transLoc , filesep , subID , '_UNMC_', tmpNight{1},'.mat'];

        load(transDat,'TT');

        % Transfer Scores to FSScore

        % Resize 
        if height(TT) > height(LFPTT) % if scored is taller
            tmpScore = TT.UNMC(1:height(LFPTT));
        elseif height(TT) < height(LFPTT) % if scored is shorter
            tmpScore = TT.UNMC;
            LFPTT = LFPTT(1:height(TT),:);
        else
            tmpScore = TT.UNMC;
        end

        LFPTT.FSScore = tmpScore;

        % Rename table and save
        saveName = [subID,'_UNMC_',tmpNight{1},'_NN.mat'];

        save(saveName,"LFPTT")

    end
end