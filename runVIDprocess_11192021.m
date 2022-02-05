%% Steps for automation

% 1. Folder with data
% 2. Folder for epoch data

% i. Check Raw folder
%    a. if more than one .mp4
%        1. Concatenate
%        2. Save in FullNight folder
%        3. Remove .api extension and save
%    b. if one .mp4
%        1. Save in FullNight folder

% ii. Run Epoch create on data in FullNight folder
vidLOC = 'I:\01_Coding_Datasets\SLEEP VIDEO\UNMC_11\allVid';
epFold = 'I:\01_Coding_Datasets\SLEEP VIDEO\UNMC_11\allVid\VideoDATA';

cd(vidLOC)
vidDir1 = dir('*.mp4');
vidDir2 = {vidDir1.name};
if length(vidDir2) > 1
    [saveDIR] = concatVIDs(vidLOC);
else
    oldLOC =  [vidLOC , filesep , vidDir2{1}];
    saveDIR = [vidLOC , filesep , 'FullNight'];
    if ~exist(saveDIR, 'dir')
        mkdir(saveDIR)
    end
    newLOC =  [saveDIR , filesep , vidDir2{1}];
    movefile(oldLOC , newLOC)
end


%% Create Epochs
VideoEpochCreate_01272022(saveDIR , epFold , 7)

%% Create LFP/EEG .mat files

%% Create final struct synced to video
clear
eegLOC = 'I:\01_Coding_Datasets\SLEEP_Score_Consensus\FinalDATA\NIGHT_2';
lfpLOC = [eegLOC , filesep , 'LfpDATA'];
vidLOC = [eegLOC , filesep , 'VideoDATA'];

% Load Night 1 - EEG Mat
cd(eegLOC)
load("2_UNMC_2.mat","TT")
% Load Night 1 - LFP Mat
cd(lfpLOC)
load("2_UNMC_2_LFP.mat","LFPTT")
% Load number of epochs in Video Folder [Night 1]
cd(vidLOC)
mList = dir('*.mat');
VidNumInd = {mList.name};


%%

maxVIDePOCH = length(VidNumInd);
maxEEGePOCH = height(TT);

if maxVIDePOCH < maxEEGePOCH
    TT = TT(1:maxVIDePOCH,:);
    LFPTT = LFPTT(1:maxVIDePOCH,:);
    % resort vid names
    numVals = cellfun(@(x) str2double(extractBetween(x,'_','.')),...
        VidNumInd,'UniformOutput',true);
    [~,newINDv] = sort(numVals);
    VidNumUse = VidNumInd(newINDv);

    TT.VideoE = transpose(VidNumUse);
    LFPTT.VideoE = transpose(VidNumUse);

else
    remIND = transpose(maxEEGePOCH+1:1:maxVIDePOCH);
    remBool = ~ismember(VidNumInd,remIND);
    VidNumInd = VidNumInd(remBool,:);
    VidSelIndLcs = VidSelIndLcs(remBool,:);
end

%% Save data
% Load Night 1 - EEG Mat
cd(eegLOC)
save("2_UNMC_1.mat","TT")
% Load Night 1 - LFP Mat
cd(lfpLOC)
save("2_UNMC_1_LFP.mat","LFPTT")

%% Create downsampled LFP [250Hz]
clear
eegLOC = 'I:\01_Coding_Datasets\SLEEP_Score_Consensus\FinalDATA\NIGHT_3';
lfpLOC = [eegLOC , filesep , 'LfpDATA'];
vidLOC = [eegLOC , filesep , 'VideoDATA'];
cd(lfpLOC)
load("2_UNMC_3_LFP.mat","LFPTT")

%%
LFPTT2 = LFPTT;
conTACTS = ["0","1","2","3"];
for ci = 1:length(conTACTS)

    for ei = 1:height(LFPTT)

        LFPTT2.(conTACTS{ci}){ei} = resample(LFPTT.(conTACTS{ci}){ei},1,2);
        disp(['Epoch ', num2str(ei), ' out of ',...
            num2str(height(LFPTT)), ' - Contact ' , num2str(ci), ' done!'])

    end
end
cd(lfpLOC)
save("2_UNMC_3_LFP.mat","LFPTT")

%% Check OLD v NEW
clear
newLOC = 'I:\01_Coding_Datasets\SLEEP_Score_Consensus\FinalDATA\NIGHT_1';
oldLOC = 'I:\01_Coding_Datasets\SLEEP_Score_Consensus\FinalDATA\NIGHT_1\oldV';

cd(newLOC)
load("2_UNMC_1.mat","TT");
newTT = TT;

cd(oldLOC)
load("2_UNMC_1.mat","TT")
oldTT = TT;

