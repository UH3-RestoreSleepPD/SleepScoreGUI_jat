function [] = blockTrialExtract_as(recDataLoc , movIDp , plotCheckfl , runAlign)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


cd(recDataLoc)

matDirall = dir('*.mat');
matDirnames = {matDirall.name};
dataFileName = matDirnames{contains(matDirnames,'Data')};
timeFileName = matDirnames{contains(matDirnames,'Times')};


% Load time file
load(timeFileName,'Acctable');

% Trial Names [start and stop]
% trialNames = Acctable.Properties.VariableNames;

% trialNums = extractAfter(trialNames,'0');
% startStops = extractBefore(trialNames,'0');
% totalNumTrials = numel(unique(trialNums));

% trialTypesN = sum(cellfun(@(x) ~isempty(x), table2cell(Acctable(:,1)),...
% 'UniformOutput',true));

% Load data file
load(dataFileName,'GeneratedData');

% Get all time vector
allTime = GeneratedData.localTime;
% Get all timezome
Timezonechange = allTime.TimeZone;

%%%% GET TIME VECTOR IN SAME DIMENMSION
tempNAN = ~isnan(GeneratedData.Accel_XSamples);
allTime2 = allTime(tempNAN);

% Calculate indicies for time offsets for trials/blocks
useTimeTab = createTABtimes(Acctable , Timezonechange , allTime2);

% Loop through and create cell array of raw data 3D array
[outRaw] = extractRAWdata(useTimeTab , GeneratedData);

if plotCheckfl
    plotCHECK(outRaw , movIDp , useTimeTab)
end

% Align trials within block
if runAlign
    [alignRaw] = reAlignTrials(outRaw);
end







end





function [outTable] = createTABtimes(inTABLE , inTimeZone , dataTIMES)

tab2celVec = reshape(table2cell(inTABLE),numel(inTABLE),1);
fillCellC = sum(cellfun(@(x) ~isempty(x), tab2celVec,...
    'UniformOutput',true));

% Loop through columns
moveIDS = 1;
allCount = 0;
allTIMES = NaT(fillCellC/2,2);
allrowNames = cell(fillCellC/2,1);

% REORDER by TRIAL and

for ti = 1:width(inTABLE)
    tmpTabRow = inTABLE(ti,:);
    tmpTabData = table2array(inTABLE(ti,:));
    trialTabIndicies = cellfun(@(x) ~isempty(x), tmpTabData,...
        'UniformOutput',true);
    rowName = tmpTabRow.Properties.VariableNames(trialTabIndicies);
    timeEvents = extractBefore(rowName,'0');
    trialNums = extractAfter(rowName,'0');
    for bi = 1:length(rowName)

        if mod(bi,2) == 1
            allCount = allCount + 1;
        end
        switch timeEvents{bi}
            case 'Start'

                allTIMES(allCount,1) = tmpTabData{bi};

            case 'Stop'
                allTIMES(allCount,2) = tmpTabData{bi};

        end
        moveRowName = ['M',num2str(moveIDS),'_T',trialNums{bi}];
        allrowNames{allCount} = moveRowName;
    end
    moveIDS = moveIDS + 1;
end


allTimeINDs = zeros(size(allTIMES));
for rri = 1:height(allTIMES)

    for cci = 1:2
        tmpTime = allTIMES(rri,cci);

        tmpTime.TimeZone = inTimeZone;
        [~ , tmpTimeInd] = min(abs(dataTIMES - tmpTime));
        allTimeINDs(rri,cci) = tmpTimeInd;

    end
end

tabTimes = array2table(allTIMES,'VariableNames',{'StartTime','StopTime'});
tabInds = array2table(allTimeINDs,'VariableNames',{'StartInd','StopTInd'});
tabTrID = table(allrowNames,'VariableNames',{'TrialID'});


outTable = [tabTimes , tabInds , tabTrID];


end % function end





function [outRaw] = extractRAWdata(inTABLE , RAWdata)


outRaw = cell(height(inTABLE),1);
for ii = 1:height(inTABLE)

    startIND = inTABLE.StartInd(ii);
    stopIND = inTABLE.StopTInd(ii);

    xSample = RAWdata.Accel_XSamples(startIND:stopIND);
    ySample = RAWdata.Accel_YSamples(startIND:stopIND);
    zSample = RAWdata.Accel_ZSamples(startIND:stopIND);

    finAccelMat = [xSample , ySample , zSample];

    cleanNans = ~isnan(xSample);

    finAccelMatO = finAccelMat(cleanNans,:);

    outRaw{ii} = finAccelMatO;

end




end



function [outRAW] = reAlignTrials(inRAW)



%%% FIRST ATTEMPT - Realign based on first peak
% test block 1
blockLENGTHS = [160, 110, 160, 110];
blockStrT = [49, 39, 49, 39];
blockEndT = [110, 70, 110, 70];

allblocks = cell(1,4);
for bxb = 1:4
    blocktemp = epochDATA{bxb};
    % initial peak
    block1x = zeros(3,blockLENGTHS(bxb));
    for x3x = 1:3
        xtemp = blocktemp{x3x}(1,:);
        xtm = abs(abs(xtemp) - mean(abs(xtemp)));
        % plot(xtm)

        [~,peakLoc] = findpeaks(xtm,"NPeaks",1,"MinPeakHeight",30);
        trimStart = peakLoc - blockStrT(bxb);
        trimEnd = peakLoc + blockEndT(bxb);
        block1x(x3x,:) = xtm(trimStart:trimEnd);
    end

    % close all
    figure;
    plot(transpose(block1x))
    allblocks{bxb} = block1x;


end




end















function [] = plotCHECK(inDATA , moveNUM , inTable)


indiciesForMove = contains(extractBefore(inTable.TrialID,'_'),num2str(moveNUM));

tmpDATA = inDATA(indiciesForMove);

for ti = 1:height(tmpDATA)

    tmpX = tmpDATA{ti}(:,1);
    hold on
    plot(tmpX)

end

end
