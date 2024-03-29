function [] = organizeHypoData()

final_1_LOC = 'D:\VRAWLFP\FinalScores';
proc_LOC = 'D:\VRAWLFP\ProcLFP';
hypno_LOC = 'D:\VRAWLFP\FinalHypnogram';

cd(final_1_LOC)
fmat1 = dir('*.mat');
fmata = {fmat1.name};

for fi = 1:length(fmata)

    cd(final_1_LOC)
    ctmp = fmata{fi};
    caseName = ctmp;
    nameInfo = split(caseName,'_');
    tmpCfind = [nameInfo{2},'_',nameInfo{1}];
    load(ctmp,'finalCS');

    % Add proc_LOC and sub folder
    [caseFolder] = getFolderC(proc_LOC,tmpCfind);
    if isnan(caseFolder)
        continue
    end

    proc_LOCsub = [proc_LOC , filesep , caseFolder];
    % Loop throug LFP nights
    [nightNum , nightName] = getNightInfo(proc_LOCsub);

    cd(proc_LOCsub)

    nightAll = struct;
    for nni = 1:length(nightNum)
        load(nightName{nni},'LFPTTRaw')
        switch nni
            case 1

                if height(LFPTTRaw) > height(finalCS.N1)
                    LFPTTRaw = LFPTTRaw(1:height(finalCS.N1),:);
                elseif height(LFPTTRaw) < height(finalCS.N1)
                    finalCS.N1 = finalCS.N1(1:height(LFPTTRaw),:);
                end

                LFPTTRaw.SleepScore = finalCS.N1;
                LFPTTRaw = LFPTTRaw(:,5:6);
                nightAll.(['Night',num2str(nni)]) = LFPTTRaw;

            case 2

                if height(LFPTTRaw) > height(finalCS.N2)
                    LFPTTRaw = LFPTTRaw(1:height(finalCS.N2),:);
                elseif height(LFPTTRaw) < height(finalCS.N2)
                    finalCS.N2 = finalCS.N2(1:height(LFPTTRaw),:);
                end

                LFPTTRaw.SleepScore = finalCS.N2;
                LFPTTRaw = LFPTTRaw(:,5:6);
                nightAll.(['Night',num2str(nni)]) = LFPTTRaw;

            case 3

                if height(LFPTTRaw) > height(finalCS.N3)
                    LFPTTRaw = LFPTTRaw(1:height(finalCS.N3),:);
                elseif height(LFPTTRaw) < height(finalCS.N3)
                    finalCS.N3 = finalCS.N3(1:height(LFPTTRaw),:);
                end

                LFPTTRaw.SleepScore = finalCS.N3;
                LFPTTRaw = LFPTTRaw(:,5:6);
                nightAll.(['Night',num2str(nni)]) = LFPTTRaw;

        end
    end
    %-- Save the night %%

    % Cd to save loc
    cd(hypno_LOC)

    % Save name with case
    saveName = [tmpCfind , '_hypnoTime.mat'];

    % Save data
    save(saveName , "nightAll");


end
end



function [caseFolder] = getFolderC(procLOCf,caseID)

cd(procLOCf)
folddir1 = dir();
folddir2 = {folddir1.name};
folddir3 = folddir2(~ismember(folddir2,{'..','.'}));

if sum(matches(folddir3,caseID)) == 0
    caseFolder = nan;
else
    caseFolder = folddir3{matches(folddir3,caseID)};
end

end


function [nightNum , nightName] = getNightInfo(proc_LOCsub)

cd(proc_LOCsub)
matdir1 = dir('*.mat');
matdir2 = {matdir1.name};
firstSet = matdir2(contains(matdir2,'LFPrawNF'));

nightNum = transpose(1:3);
nightName = cell(3,1);

for fi = 1:length(firstSet)

    tmpN = firstSet{fi};
    namePARTS = split(tmpN,'_');
    nightNNustr = str2double(namePARTS{3});

    switch nightNNustr
        case 1
            nightName{fi} = tmpN;
        case 2
            nightName{fi} = tmpN;
        case 3
            nightName{fi} = tmpN;
    end

end
end