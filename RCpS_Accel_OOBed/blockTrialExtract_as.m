function [] = blockTrialExtract_as(recDataLoc , movIDp , plotCheckfl , runAlign , UinTable)


arguments

    recDataLoc  (1,:) char     = 'Default'
    movIDp      (1,:) double   = NaN;
    plotCheckfl (1,:) logical  = 0
    runAlign    (1,:) logical  = 0
    UinTable    (:,:) table    = Nan;

end


if matches(recDataLoc,'Default')
    cd('C:\Users\Admin\Documents\Github\SleepScoreGUI_jat\RCpS_Accel_OOBed\Thompson_0622\Trial_1')
else
    cd(recDataLoc)
end

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
if isnan(UinTable{1,1})
    useTimeTab = createTABtimes(Acctable , Timezonechange , allTime2);
else
    useTimeTab = UinTable;
end

% Loop through and create cell array of raw data 3D array
[outRaw] = extractRAWdata(useTimeTab , GeneratedData);

if plotCheckfl
    plotCHECK(outRaw , movIDp , useTimeTab)
end

% Align trials within block
if runAlign
    [alignRaw] = reAlignTrials(outRaw , useTimeTab);
end


[move1_pca] = getPCAloadings(alignRaw , 1); % Standing basline
[move2_pca] = getPCAloadings(alignRaw , 2); % Standing to Sitting
[move3_pca] = getPCAloadings(alignRaw , 3); % Sitting to Laying
[move4_pca] = getPCAloadings(alignRaw , 4); % Standing to Walking 
% Slice down X Y Z down to event and stack into column vectors
close all

[pcaCO_12] = comparePCAS(move1_pca , move2_pca); % 1 should be different to 2
[pcaCO_13] = comparePCAS(move1_pca , move3_pca); % 1 should be different to 3
[pcaCO_14] = comparePCAS(move1_pca , move4_pca); % 1 should be similar to 4
[pcaCO_23] = comparePCAS(move2_pca , move3_pca); % 2 should be similar to 3
[pcaCO_24] = comparePCAS(move2_pca , move4_pca); % 2 should be different from 4
[pcaCO_34] = comparePCAS(move3_pca , move4_pca); % 3 should be different from 4


figure;
bar([pcaCO_12,...
     pcaCO_13,...
     pcaCO_14,...
     pcaCO_23,...
     pcaCO_24,...
     pcaCO_34])
ylabel('Disparity')
xticklabels({'Stand vs Sit2Stand',...
             'Stand vs LieDSit',...
             'Stand vs Walking',...
             'Sit2Stand vs LieDSit',...
             'Sit2Stand vs Walking',...
             'LieDSit vs Walking'});
axis square





figure;
plot3(move1_pca.accelDataDBc(:,1),move1_pca.accelDataDBc(:,2),move1_pca.accelDataDBc(:,3),'ko')
hold on
text(mean(move1_pca.accelDataDBc(:,1))+10,mean(move1_pca.accelDataDBc(:,2)),...
    mean(move1_pca.accelDataDBc(:,3)),'STANDING','FontSize',20,'FontWeight','bold','Color','k')


plot3(move2_pca.accelDataDBc(:,1),move2_pca.accelDataDBc(:,2),move2_pca.accelDataDBc(:,3),'ro')
text(mean(move2_pca.accelDataDBc(:,1))+10,mean(move2_pca.accelDataDBc(:,2)),...
    mean(move2_pca.accelDataDBc(:,3)),'Sit2Stand','FontSize',20,'FontWeight','bold','Color','r')


plot3(move3_pca.accelDataDBc(:,1),move3_pca.accelDataDBc(:,2),move3_pca.accelDataDBc(:,3),'go')
text(mean(move3_pca.accelDataDBc(:,1))+10,mean(move3_pca.accelDataDBc(:,2)),...
    mean(move3_pca.accelDataDBc(:,3)),'LieDown2Sit','FontSize',20,'FontWeight','bold','Color','g')


plot3(move4_pca.accelDataDBc(:,1),move4_pca.accelDataDBc(:,2),move4_pca.accelDataDBc(:,3),'bo')
text(mean(move4_pca.accelDataDBc(:,1))+10,mean(move4_pca.accelDataDBc(:,2)),...
    mean(move4_pca.accelDataDBc(:,3)),'Walking','FontSize',20,'FontWeight','bold','Color','b')




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
tabInds = array2table(allTimeINDs,'VariableNames',{'StartIndraw','StopTIndraw'});
movementIDSall = extractBefore(allrowNames,'_');
trialIDSall = cellfun(@(x) str2double(extractAfter(extractAfter(x,'_'),'T')),...
    allrowNames,'UniformOutput',true);
tabTrID = table(allrowNames, movementIDSall, trialIDSall,'VariableNames',{'FullTrialID','MoveID','TrialID'});


outTable = [tabTimes , tabInds , tabTrID];


end % function end





function [outRaw] = extractRAWdata(inTABLE , RAWdata)


outRaw = cell(height(inTABLE),1);
for ii = 1:height(inTABLE)

    startIND = inTABLE.StartIndraw(ii);
    stopIND = inTABLE.StopTIndraw(ii);

    noneNANrawX = RAWdata.Accel_XSamples(~isnan(RAWdata.Accel_XSamples));
    noneNANrawY = RAWdata.Accel_YSamples(~isnan(RAWdata.Accel_YSamples));
    noneNANrawZ = RAWdata.Accel_ZSamples(~isnan(RAWdata.Accel_ZSamples));

    xSample = noneNANrawX(startIND:stopIND);
    ySample = noneNANrawY(startIND:stopIND);
    zSample = noneNANrawZ(startIND:stopIND);

    finAccelMat = [xSample , ySample , zSample];

    outRaw{ii} = finAccelMat;

end




end



function [allblocksfinal] = reAlignTrials(inRAW , inInfoTab)

% TRIM by MIN

% Combine and Fix STANDING for 1 minute

% BLOCKS equals types of MOVEMENTs
mUNInum = tabulate(inInfoTab.MoveID);

allblocks = cell(3,height(mUNInum));
for bxb = 1:height(mUNInum)
    % initial peak
    tmpMoveID = mUNInum{bxb,1};
    numTrials = mUNInum{bxb,2};
    rawDataCell = inRAW(matches(inInfoTab.MoveID,tmpMoveID));

    if numel(rawDataCell) > 1
        blockLENGTHS = max(cellfun(@(x) height(x), rawDataCell, 'UniformOutput',true));
    else
        blockLENGTHS = height(rawDataCell{1});
    end

    block1x = nan(numTrials,blockLENGTHS);
    block1y = nan(numTrials,blockLENGTHS);
    block1z = nan(numTrials,blockLENGTHS);

    for x3x = 1:numTrials
        xtemp = rawDataCell{x3x}(:,1);
        ytemp = rawDataCell{x3x}(:,2);
        ztemp = rawDataCell{x3x}(:,3);

        xtm = abs(xtemp - mean(xtemp)) - min(abs(xtemp - mean(xtemp)));
        % plot(xtm)
        % plot(xtemp)

        [~,peakLoc] = findpeaks(xtm,"NPeaks",1,"MinPeakHeight",mean(xtm)*1.2);
        if isempty(peakLoc)
            trimStart = 1;
        else
            trimStart = peakLoc;
        end
        trimEnd = length(xtemp);
        trimLength = numel(trimStart:trimEnd);
        block1x(x3x,1:trimLength) = transpose(xtemp(trimStart:trimEnd));
        block1y(x3x,1:trimLength) = transpose(ytemp(trimStart:trimEnd));
        block1z(x3x,1:trimLength) = transpose(ztemp(trimStart:trimEnd));
    end

    % close all
    % figure;
    % plot(transpose(block1x))
    % allblocks{bxb} = block1x;

    allblocks{1,bxb} = block1x;
    allblocks{2,bxb} = block1y;
    allblocks{3,bxb} = block1z;

end

% FIX M1 and M5
% hold on
% plot(allblocks{1,5})
% plot(allblocks{2,5})
% plot(allblocks{3,5})
% plot(mean(cell2mat(allblocks(:,5))),'k')
move1and5 = [1,5];
tableMat = zeros(2,3);
for m15 = 1:2

    % meanLine = mean(cell2mat(allblocks(:,move1and5(m15))));
    % noMeanLine = meanLine(~isnan(meanLine));
    % rmsThresh = (rms(noMeanLine)*0.1) + rms(noMeanLine);
    % 
    % startI = 1;
    % stopI = 5;
    % stePSize = 5;
    % stePNum = floor(length(noMeanLine)/stePSize);
    % stepIds = zeros(1,stePNum);
    % allRMS = zeros(1,stePNum);
    % for si = 1:stePNum
    % 
    %     tmpRMS = rms(noMeanLine(startI:stopI));
    %     allRMS(si) = tmpRMS;
    %     if tmpRMS > rmsThresh
    %         stepIds(si) = 1;
    %     end
    % 
    %     startI = startI + stePSize;
    %     stopI = stopI + stePSize;
    % 
    % end
    % 
    % offsetThresh = diff(stepIds);
    % numOffset = sum(numel(find(offsetThresh)));
    % if numOffset > 5
    %     firstStable = find(stepIds,1,'first') * stePSize;
    %     lastStable = find(stepIds,1,'last') * stePSize;
    %     % plot(noMeanLine(firstStable:lastStable))
    % 
    %     tableMat(m15,1) = firstStable;
    %     tableMat(m15,2) = lastStable;
    %     tableMat(m15,3) = length(firstStable:lastStable);
    % 
    % else
    %     offsetLOCS = find(offsetThresh);
    %     offSETfracs = offsetLOCS/stePNum;
    %     if max(offSETfracs) < 0.2
    %         firstStable = max(offsetLOCS) * stePSize;
    %         lastStable = stePNum*stePSize;
    %         % plot(noMeanLine(firstStable:end));
    %     else
    %         firstStable = 1;
    %         lastStable = min(offsetLOCS) * stePSize;
    %         % plot(noMeanLine(1:lastStable))
    %     end

        tableMat(m15,1) = 10;%firstStable;
        tableMat(m15,2) = inInfoTab.StopTIndraw(m15) - inInfoTab.StartIndraw(m15) - 10;
        tableMat(m15,3) = length(tableMat(m15,1):tableMat(m15,2));

    % end

end

m15Infotab = array2table(tableMat,'VariableNames',{'startIND','stopIND','Numels'});

[minLENGTH , ~] = min(m15Infotab.Numels);

bothBlocks = cell(3,2);
for bii = 1:2
    tmpBlock = allblocks(:,move1and5(bii));
    startINDi = m15Infotab.startIND(bii);
    stopINDi = startINDi + minLENGTH - 1;
    tmpBlockt = cellfun(@(x) x(startINDi:stopINDi) , tmpBlock, 'UniformOutput' , false);
    bothBlocks(:,bii) = tmpBlockt;
end

combineBlocks = cellfun(@(x,y) [x ; y], bothBlocks(:,1) , bothBlocks(:,2),...
    'UniformOutput' , false);

allblocksfinal = cell(3,4);
allblocksfinal(:,1) = combineBlocks;
allblocksfinal(:,2:4) = allblocks(:,2:4);


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









function [pcaOUTdata] = getPCAloadings(allBLOCKs , movementNUM)
% Assuming accelData is your matrix with three columns (X, Y, Z)

% Convert to X , Y , Z - with 3 trials along rows

tempBlOCK = allBLOCKs(:,movementNUM);

Xaccel = reshape(transpose(tempBlOCK{1}),numel(tempBlOCK{1}),1);
Yaccel = reshape(transpose(tempBlOCK{2}),numel(tempBlOCK{2}),1);
Zaccel = reshape(transpose(tempBlOCK{3}),numel(tempBlOCK{3}),1);

accelData = [Xaccel , Yaccel , Zaccel];

accelDataNan = accelData(~isnan(accelData(:,1)),:);

pdistAll = zeros(height(accelDataNan),1);
for pi = 1:height(accelDataNan)

    D = pdist2(accelDataNan(pi,:),accelDataNan);

    minD = min(D(D ~= 0));
    pdistAll(pi) = minD;
end

pdistALLs = sort(pdistAll,'descend');
plot(pdistALLs)

%%%% THIS WORKS for 1 and 2
epsilon = mean(pdistALLs,'omitnan') + (std(pdistALLs,'omitnan'));

if numel(pdistALLs) < 100
    numGroupMems = 5;
else
    numGroupMems = 15;
end

% close all
idx_DBSCAN = dbscan(accelDataNan,epsilon,numGroupMems); 
% gscatter(accelData(:,1),accelData(:,2),idx);
% plot3(accelData(:,1),accelData(:,2),accelData(:,3),'k.')
% hold on
% plot3(accelData(idx == 1,1),accelData(idx == 1,2),accelData(idx == 1,3),'r.')
% tabulate(idx)
% title('DBSCAN Using Euclidean Distance Metric')
% Standardize the data to have zero mean and unit variance
tabtab = tabulate(idx_DBSCAN);
% extractNon-1 
clusterTab = tabtab(tabtab(:,1) ~= -1,:);
% find max ID
[~, maxRow] = max(clusterTab(:,2));
maxID = clusterTab(maxRow,1);

accelDATA_dbscan = accelDataNan(idx_DBSCAN == maxID,:);
accelData_standardized = zscore(accelDATA_dbscan);

% Perform PCA
[coeff,score,~,~,explained] = pca(accelData_standardized);

pcaOUTdata.movNum = movementNUM;
pcaOUTdata.accelData = accelDataNan;
pcaOUTdata.accelDataDBc = accelDATA_dbscan;
pcaOUTdata.loading = coeff;
pcaOUTdata.score = score;
pcaOUTdata.explained = explained;

% The output 'coeff' contains the principal component vectors (loadings).

end








function [pcaCOMPaRE] = comparePCAS(move1 , move2)


% Perform Procrustes Analysis
[pcaCOMPaRE, ~, ~] = procrustes(move1.loading, move2.loading);


end

