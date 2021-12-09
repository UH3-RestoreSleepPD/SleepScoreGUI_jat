function [] = convertLFP_OldMINN(folderLOC,patID,saveLoc)

cd(folderLOC)

mdir = dir('*.mat');
mdir2 = {mdir.name};
app.fileLIST = mdir2;

% Score data
sleepScoreR_fn = app.fileLIST{contains(app.fileLIST,'Index')};
matob_SlSc = matfile(sleepScoreR_fn);
varlist.slsc = who(matob_SlSc);
load(sleepScoreR_fn,varlist.slsc{1});
sleepScoreR = eval(varlist.slsc{1});

% Raw lfp data
lfpDataR_fn = app.fileLIST{contains(app.fileLIST,'LFP')};
matob_lfp = matfile(lfpDataR_fn);
varlist.lfp = who(matob_lfp);
load(lfpDataR_fn,varlist.lfp{1});
lfpDataR = eval(varlist.lfp{1});

disp('data loaded!')

% Get DBS data
labInds = ismember(lfpDataR.montage,{'LFP0','LFP1','LFP2','LFP3'});
lfpCols = lfpDataR.data(:,labInds);
lfpFs = lfpDataR.Fs;

% DETERMINE ID for LFP ----------------------------------------
chanLab_LFP = {'DBS','DBS','DBS','DBS'};
chanID_LFP = {'0','1','2','3'};
aoID_LFP = {'LFP0','LFP1','LFP2','LFP3'};

labelTab_LFP = table(transpose(aoID_LFP),transpose(chanLab_LFP),...
    transpose(chanID_LFP),'VariableNames',{'AO_ID',...
    'Label','ID'});

% Convert to Time tables with appropriate sample frequencies
chansleep_LFP = {'DBS'};
chanNUMb_LFP = zeros(1,1);
for cNi = 1:length(chansleep_LFP)
    chanNUMb_LFP(cNi) = sum(ismember(labelTab_LFP.Label,chansleep_LFP{cNi}));
end

lfpMODsf = lfpFs;

tmpLFParray = zeros(size(lfpCols)); % 14 Hours
chanTabLabs_LFP = cell(1,1);
aoTabLabs_LFP = cell(1,1);

chanChooselfp = chansleep_LFP{1};
chanTabLabs_LFP{1} = labelTab_LFP.ID(ismember(labelTab_LFP.Label,chanChooselfp));
aoTabLabs_LFP{1} = labelTab_LFP.AO_ID(ismember(labelTab_LFP.Label,chanChooselfp));

for chi = 1:length(aoTabLabs_LFP{1})

    tmpChanlfp = aoTabLabs_LFP{1}(ismember(aoTabLabs_LFP{1},aoTabLabs_LFP{1}(chi)));
    tmpRowslfp = transpose(matches(aoTabLabs_LFP{1},tmpChanlfp));
    tmpLFParray(:,chi) = lfpCols(:,tmpRowslfp);

end

tmpLFPtab = array2table(tmpLFParray);
tmpLFPtab.Properties.VariableNames = chanTabLabs_LFP{1};
lfpTTraw = table2timetable(tmpLFPtab,'SampleRate',lfpMODsf);
disp('Clean up done!')

% Normalize
disp('Normalize Start');
lfpTTrawN = normalize(lfpTTraw);
disp('Normalize Done');

% Notch Filter
disp('Notch Filter Start');
lfpNOTCH = lfpTTrawN;
for li = 1:size(lfpTTrawN,2)
    tmpLFPraw = table2array(lfpTTrawN(:,li));
    tmpLFPnotch = spectrumInterpolation(tmpLFPraw, lfpMODsf, 60, 3, 1);
    tmpLFPntab = array2table(tmpLFPnotch);
    lfpNOTCH(:,li) = tmpLFPntab;
end
disp('Notch Filter Done');

% High pass Filter
disp('High Pass Filter Start');
lfpTTppHp = lfpNOTCH;
for li = 1:size(lfpTTppHp,2)
    lfpTTpTmp = table2array(lfpTTppHp(:,li));
    lfpTTppHp{:,li} = highpass(lfpTTpTmp,0.2,1024,'ImpulseResponse','iir','Steepness',0.8);
end
disp('High Pass Filter Done');

% Low pass Filter
disp('Low Pass Filter Start');
lfpTTppLp = lowpass(lfpTTppHp,250); % For algorithm
lfpTTppSm = lfpTTppLp;
disp('Low Pass Filter Done');

% Downsample
% CREATE RC+S version
% CREATE RAW version - no downsampling
disp('Downsample 250Hz Start');
lfpTTppDs = retime(lfpTTppSm,'regular','mean','SampleRate',250);
disp('Downsample 250Hz Done');

% Preprocess Signal
app.lfpTTpp = lfpTTppDs;
app.lfpTTNdnS = lfpTTppHp;

% Create Bins
lfpTABf = createBinsNonGUI(app,'lfp',chanTabLabs_LFP{1},250);
lfpTABfND = createBinsNonGUI(app,'lfpr',chanTabLabs_LFP{1},1024);
disp('Create Bins Done!')

% Create sleep score bins
sleepScoreCA = cell(height(lfpTABf) , 1);
if height(sleepScoreCA) > height(sleepScoreR)
    lfpTABf = lfpTABf(1:height(sleepScoreR),:);
    lfpTABfND = lfpTABfND(1:height(sleepScoreR),:);
    sleepScoreCA = cell(height(lfpTABf) , 1);
elseif height(sleepScoreCA) < height(sleepScoreR)
    sleepScoreR = sleepScoreR(1:height(sleepScoreCA),:);
end

for si = 1:height(lfpTABf)

    tmpSS = sleepScoreR(si,3);
    switch tmpSS
        case -1
            sleepScoreCA{si} = 'W';
        case 0
            sleepScoreCA{si} = 'W';
        case 1
            sleepScoreCA{si} = 'N1';
        case 2
            sleepScoreCA{si} = 'N2';
        case 3
            sleepScoreCA{si} = 'N3';
        case 5
            sleepScoreCA{si} = 'R';
    end

end
% Create FINAL LFP Timetables -------------------------
LFPTT = lfpTABf;
LFPTT.FSScore = sleepScoreCA;  
LFPTTRaw = lfpTABfND;
LFPTTRaw.FSScore = sleepScoreCA;  

% SAVE OUT
cd(saveLoc)

app.fname = [patID,'_UMin_1'];
app.fnameMAT = [app.fname,'.mat'];
app.fnameMATlfp = [app.fname,'_LFP.mat'];
app.fnameMATlfpR = [app.fname,'_LFPraw.mat'];

save(app.fnameMATlfp,'LFPTT')
save(app.fnameMATlfpR,'LFPTTRaw')

end





