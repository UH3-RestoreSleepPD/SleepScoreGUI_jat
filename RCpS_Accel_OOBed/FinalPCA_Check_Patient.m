%% Accelerometry prediction

% 1. Find cluster IDs from clinic testing


% 2. Move through collected accel overnight 


% 3. In 30 second chunks determine if clusterd data falls near one of the
% predefined centroids

%% Clinic data
cd('E:\Dropbox\Platt_accel_testing_09172023\UN12-PrelimData')
% load("UN12_07-Aug-2023_setting.mat")

%% 

close all
rLC2 = 'E:\Dropbox\Platt_accel_testing_09172023\UN12-PrelimData';

movePCAs = blockTrialExtract_as(rLC2,2,0,1);
close all
% accelDataDBc

% Get centroids
moveCentroids = zeros(4,3);
for mi = 1:4

    moveID = ['move',num2str(mi)];

    moveCentroids(mi,:) = mean(movePCAs.(moveID).accelDataDBc);

end

%%
accelDATAall = [X_axisData , Y_axisData , Z_axisData];
% 30 hz
% 30 x 30 = 30 seconds

% clean all nans
moveMoves = {'baseline','Sit2Stand','LieDown2Sit','Walking'};
accelDAll = accelDATAall(~isnan(accelDATAall(:,1)),:);

% get steps
numSamples = 30*30;
startsI = transpose(floor(linspace(1,height(accelDAll),floor(height(accelDAll)/numSamples))));
starts = startsI(1:length(startsI)-1);
stops = [starts(2:length(starts))-1 ; startsI(length(startsI))];
indexStSts = [starts , stops];

% Process each step
moveAllsegs = cell(height(indexStSts),1);
moveAllids = zeros(height(indexStSts),1);
for ii = 1:height(indexStSts)

    stepIi = accelDAll(indexStSts(ii,1):indexStSts(ii,2),:);
    [pcaOUTdata] = getPCAloadings(stepIi);

    pcaTEMPcen = mean(pcaOUTdata.accelDataDBc);

    moveDists = pdist2(pcaTEMPcen,moveCentroids);
    
    
    [~, move2id] = min(moveDists);

    moveAllsegs{ii} = moveMoves{move2id};
    moveAllids(ii) = move2id;
    disp([num2str(ii) , ' ' , moveMoves{move2id}])
end



plot(moveAllids)
test = 1;





















function [pcaOUTdata] = getPCAloadings(accelData)
% Assuming accelData is your matrix with three columns (X, Y, Z)

% Convert to X , Y , Z - with 3 trials along rows

accelDataNan = accelData(~isnan(accelData(:,1)),:);

pdistAll = zeros(height(accelDataNan),1);
for pi = 1:height(accelDataNan)

    D = pdist2(accelDataNan(pi,:),accelDataNan);

    minD = min(D(D ~= 0));
    pdistAll(pi) = minD;
end

pdistALLs = sort(pdistAll,'descend');
% plot(pdistALLs)

%%%% THIS WORKS for 1 and 2
epsilon = mean(pdistALLs,'omitnan') + (std(pdistALLs,'omitnan'));

if numel(pdistALLs) < 100
    numGroupMems = 5;
else
    numGroupMems = 15;
end

% close all
idx_DBSCAN = dbscan(accelDataNan,epsilon,numGroupMems); 
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

pcaOUTdata.accelData = accelDataNan;
pcaOUTdata.accelDataDBc = accelDATA_dbscan;
pcaOUTdata.loading = coeff;
pcaOUTdata.score = score;
pcaOUTdata.explained = explained;


end







