function [] = process_ACT_Phillips(actLOC , subID)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

cd(actLOC)
load([subID ,'_ACT_DATA.mat'],'dataTable');

% First first row of non nan
onWristData = cellfun(@(x) str2double(x), dataTable.("Off-Wrist Status"));
actTABLE = dataTable(~onWristData,:);

% Extract Unique dates
actDATES = unique(actTABLE.Date);
actDATES2 = sortrows(datetime(actDATES,'InputFormat','M/dd/yyyy'));

monthS = zeros(144,length(actDATES2));
dayS = zeros(144,length(actDATES2));
hourS = nan(144,length(actDATES2));
minuteS = nan(144,length(actDATES2));
actDAYtm = NaT(144,length(actDATES2));
LFPall = zeros(144,length(actDATES2));
stimAll = zeros(144,length(actDATES2));
timXax = [];

for li = 1:length(actDATES2)

    [monthOI,dayOI,hourOI,minuteOI,actDayT] = getDT(timeD , tzOffsN);

    % convert minute column to floor round
    minuteOIc = floor(minuteOI/10)*10;
    % combine hour , minute , second
    durFind = duration(hourOI,minuteOIc,zeros(length(minuteOIc),1));

    % search for where to align times;
    [alignIND , allBlok] = alignTime(durFind);

    monthS(alignIND,li) = monthOI;
    dayS(alignIND,li) = dayOI;
    hourS(alignIND,li) = hourOI;
    minuteS(alignIND,li) = minuteOI;
    actDAYtm(alignIND,li) = actDayT;
    LFPall(alignIND,li) = tLFP;
    stimAll(alignIND,li) = tstim_mA;
    timXax = allBlok;

end


% Process Actigraphy
% Convert into 10 minute bins
% Loop through date / hour / 10 min bins

activityMean = zeros(144,length(uniMDate));
activitySTD = zeros(144,length(uniMDate));
WLMean = zeros(144,length(uniMDate));
WLSTD = zeros(144,length(uniMDate));
RLMean = zeros(144,length(uniMDate));
RLSTD = zeros(144,length(uniMDate));
GLMean = zeros(144,length(uniMDate));
GLSTD = zeros(144,length(uniMDate));
BLMean = zeros(144,length(uniMDate));
BLSTD = zeros(144,length(uniMDate));
WakeFrac = zeros(144,length(uniMDate));
MaxInter = cell(144,length(uniMDate));

for di = 1:length(uniMDate)
    % Unique date


    % temp act table
    dconvert = cellfun(@(x) replace(x,'/','-'),actTABLE.Date,"UniformOutput",false);
    actDATEn = datetime(dconvert,'InputFormat','M-dd-uuuu');

    % Logical date index - LFP



    % 1: date, 2: Time, 3: ActivityMean, 4: ActivitySTD, 5: WhitelightMean
    % 6: WhitelightSTD, 7: RedLightMean, 8: RedLightSTD, 9:
    % GreenlightMean, 10: GreenlightSTD, 11: BluelightMean, 12:
    % BlulightSTD, 13: SleepWakefrac, 14: IntervalstatusMax
    act10minBin = cell(length(hourMINlfp2),14);
    act10minBinAll = cell(length(hourMINlfp2),1);

    % Loop through LFP times
    for lfpT = 1:length(hourMINlfp2)
        tmpLFPbin = hourMINlfp2{lfpT};
        % Remove blanks
        tmpLFPbin = replace(tmpLFPbin,' ','');
        % Start time bin
        binStart = find(matches(actINDtime,tmpLFPbin),1,'first');

        % If bin is not present because actigraphy was not recorded
        % I.E., off-wrist
        if isempty(binStart)
            act10minBin{lfpT,3} = NaN;
            act10minBin{lfpT,4} = NaN;
            act10minBin{lfpT,5} = NaN;
            act10minBin{lfpT,6} = NaN;
            act10minBin{lfpT,7} = NaN;
            act10minBin{lfpT,8} = NaN;
            act10minBin{lfpT,9} = NaN;
            act10minBin{lfpT,10} = NaN;
            act10minBin{lfpT,11} = NaN;
            act10minBin{lfpT,12} = NaN;
            act10minBin{lfpT,13} = NaN;
            continue
        else

            binSIZE = (2*10) - 1; % Bins per min * num of minutes
            binEnd = binStart + binSIZE;

            if binEnd > length(actINDtime)
                binEnd = length(actINDtime);
            end

            binTable = actDATEtab(binStart:binEnd,:);
            act10minBinAll{lfpT} = binTable;

            % date
            act10minBin{lfpT,1} = tmpDI;
            % time for 10 minute bin
            act10minBin{lfpT,2} = actHOUR24h{binStart};
            % average activity
            tmpC1 = cellfun(@(x) str2double(x),...
                binTable.Activity , 'UniformOutput',true);
            act10minBin{lfpT,3} = mean(tmpC1(~isnan(tmpC1)));
            act10minBin{lfpT,4} = std(tmpC1(~isnan(tmpC1)));
            % whitelightM
            tmpC2 = cellfun(@(x) str2double(x),...
                binTable.("White Light") , 'UniformOutput',true);
            act10minBin{lfpT,5} = mean(tmpC2(~isnan(tmpC2)));
            act10minBin{lfpT,6} = std(tmpC2(~isnan(tmpC2)));
            % redlightM
            tmpC3 = cellfun(@(x) str2double(x),...
                binTable.("Red Light") , 'UniformOutput',true);
            act10minBin{lfpT,7} = mean(tmpC3(~isnan(tmpC3)));
            act10minBin{lfpT,8} = std(tmpC3(~isnan(tmpC3)));
            % greenlightM
            tmpC4 = cellfun(@(x) str2double(x),...
                binTable.("Green Light") , 'UniformOutput',true);
            act10minBin{lfpT,9} = mean(tmpC4(~isnan(tmpC4)));
            act10minBin{lfpT,10} = std(tmpC4(~isnan(tmpC4)));
            % bluelightM
            tmpC5 = cellfun(@(x) str2double(x),...
                binTable.("Blue Light") , 'UniformOutput',true);
            act10minBin{lfpT,11} = mean(tmpC5(~isnan(tmpC5)));
            act10minBin{lfpT,12} = std(tmpC5(~isnan(tmpC5)));
            % WakeFrac
            sleepwakeT = cellfun(@(x) str2double(x),...
                binTable.("Sleep/Wake") , 'UniformOutput',true);
            act10minBin{lfpT,13} = sum(sleepwakeT)/length(sleepwakeT);
            % IntervalstatusMax
            uniStates = unique(binTable.("Interval Status"));
            numStates = zeros(length(uniStates));
            for ni = 1:length(uniStates)
                numStates(ni) = sum(matches(binTable.("Interval Status"),uniStates{ni}));
            end
            [~,maxI] = max(numStates);
            maxSTATE = uniStates{maxI};

            act10minBin{lfpT,14} = maxSTATE;

        end % End of if statement checking bin presence
    end % End of bin loop
    % Finalize table
    activityMean(:,di) = cell2mat(act10minBin(:,3));
    activitySTD(:,di) = cell2mat(act10minBin(:,4));
    WLMean(:,di) = cell2mat(act10minBin(:,5));
    WLSTD(:,di) = cell2mat(act10minBin(:,6));
    RLMean(:,di) = cell2mat(act10minBin(:,7));
    RLSTD(:,di) = cell2mat(act10minBin(:,8));
    GLMean(:,di) = cell2mat(act10minBin(:,9));
    GLSTD(:,di) = cell2mat(act10minBin(:,10));
    BLMean(:,di) = cell2mat(act10minBin(:,11));
    BLSTD(:,di) = cell2mat(act10minBin(:,12));
    WakeFrac(:,di) = cell2mat(act10minBin(:,13));
    MaxInter(:,di) = act10minBin(:,14);

end % End of date loop

cd(saveLOC)
if matches(inPS.hemiS,'L')
    fileNAME = ['SPPD',num2str(inPS.subID),'_L_TimeLine.csv'];
else
    fileNAME = ['SPPD',num2str(inPS.subID),'_R_TimeLine.csv'];
end
writetable(outTable,fileNAME);

LFPaMAT = reshape(LFPaCOL,144,size(hourS,2));

outMAT.actTime = actDAYtm;
outMAT.month = monthS;
outMAT.day = dayS;
outMAT.hour = hourS;
outMAT.minu = minuteS;
outMAT.LFP = LFPaMAT;
outMAT.Stim = stimAll;
outMAT.TimeX = cellstr(datestr(timXax));
outMAT.senseChan = activeSenseChan;
outMAT.senseFreq = activeSenseFreq;
outMAT.ActMean = activityMean;
outMAT.ActSTD = activitySTD;
outMAT.WLMean = WLMean;
outMAT.WLSTD = WLSTD;
outMAT.RLMean = RLMean;
outMAT.RLSTD = RLSTD;
outMAT.GLMean = GLMean;
outMAT.GLSTD = GLSTD;
outMAT.BLMean = BLMean;
outMAT.BLSTD = BLSTD;
outMAT.WakeFrac = WakeFrac;
outMAT.MaxState = MaxInter;

if matches(inPS.hemiS,'L')
    fileNAMEm = ['SPPD',num2str(inPS.subID),'_L_TimeLine.mat'];
else
    fileNAMEm = ['SPPD',num2str(inPS.subID),'_R_TimeLine.mat'];
end
save(fileNAMEm,'outMAT');

























end