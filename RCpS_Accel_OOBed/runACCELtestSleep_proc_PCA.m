% Load the accelerometry data into Matlab
function [] = runACCELtestSleep_proc_PCA(plottest)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


load("experBlockTrialRaw.mat",'blockSc');
epochDATA = blockSc;

load("jattimes3.mat",'Acctable');
allDATA = Acctable;

[xSampleE , ySampleE , zSampleE] = cleanSamples(allDATA);


if plottest
    runPLOT(epochDATA)
end

%%% FIRST ATTEMPT - Realign based on first peak
% test block 1
blockLENGTHS = [160, 110, 160, 110];
blockStrT = [49, 39, 49, 39];
blockEndT = [110, 70, 110, 70];

allblocks = cell(3,4); % BLOCKS equal to different movements - one PCA for each movement
for bxb = 1:4
    blocktemp = epochDATA{bxb};
    % initial peak
    block1x = zeros(3,blockLENGTHS(bxb));
    block1y = zeros(3,blockLENGTHS(bxb));
    block1z = zeros(3,blockLENGTHS(bxb));
    for x3x = 1:3
        xtemp = blocktemp{x3x}(1,:);
        ytemp = blocktemp{x3x}(2,:);
        ztemp = blocktemp{x3x}(3,:);
        xtm = abs(abs(xtemp) - mean(abs(xtemp)));
        % plot(xtm)

        [~,peakLoc] = findpeaks(xtm,"NPeaks",1,"MinPeakHeight",30);
        trimStart = peakLoc - blockStrT(bxb);
        trimEnd = peakLoc + blockEndT(bxb);
        block1x(x3x,:) = xtemp(trimStart:trimEnd);
        block1y(x3x,:) = ytemp(trimStart:trimEnd);
        block1z(x3x,:) = ztemp(trimStart:trimEnd);
    end

    % close all
    allblocks{1,bxb} = block1x;
    allblocks{2,bxb} = block1y;
    allblocks{3,bxb} = block1z;


end






[move1_pca] = getPCAloadings(allblocks , 1); % sitting to standing
[move2_pca] = getPCAloadings(allblocks , 2); % laying down
[move3_pca] = getPCAloadings(allblocks , 3); % laying down to sitting
[move4_pca] = getPCAloadings(allblocks , 4); % standing to sitting 
% Slice down X Y Z down to event and stack into column vectors


% MDS_plot(move1_pca)





[pcaCO_12] = comparePCAS(move1_pca , move2_pca); % 1 should be different to 2
[pcaCO_13] = comparePCAS(move1_pca , move3_pca); % 1 should be different to 3
[pcaCO_14] = comparePCAS(move1_pca , move4_pca); % 1 should be similar to 4
[pcaCO_23] = comparePCAS(move2_pca , move3_pca); % 2 should be similar to 3
[pcaCO_24] = comparePCAS(move2_pca , move4_pca); % 2 should be different from 4
[pcaCO_34] = comparePCAS(move3_pca , move4_pca); % 3 should be different from 4

figure;
bar([pcaCO_12,pcaCO_13,pcaCO_14,pcaCO_23,pcaCO_24,pcaCO_34])
ylabel('Disparity')
xticklabels({'Sit1 vs LieD','Sit1 vs LieDSit','Sit1 vs Sit2','LieD vs LieDSit','LieD vs Sit2', 'LieDSit vs Sit2'});
axis square


figure;
plot3(move1_pca.accelDataDBc(:,1),move1_pca.accelDataDBc(:,2),move1_pca.accelDataDBc(:,3),'ko')
hold on
text(mean(move1_pca.accelDataDBc(:,1))+10,mean(move1_pca.accelDataDBc(:,2)),...
    mean(move1_pca.accelDataDBc(:,3)),'Sit1','FontSize',20,'FontWeight','bold','Color','k')
plot3(move2_pca.accelDataDBc(:,1),move2_pca.accelDataDBc(:,2),move2_pca.accelDataDBc(:,3),'ro')
text(mean(move2_pca.accelDataDBc(:,1))+10,mean(move2_pca.accelDataDBc(:,2)),...
    mean(move2_pca.accelDataDBc(:,3)),'LieDown','FontSize',20,'FontWeight','bold','Color','r')
plot3(move3_pca.accelDataDBc(:,1),move3_pca.accelDataDBc(:,2),move3_pca.accelDataDBc(:,3),'go')
text(mean(move3_pca.accelDataDBc(:,1))+10,mean(move3_pca.accelDataDBc(:,2)),...
    mean(move3_pca.accelDataDBc(:,3)),'LieDown2Sit','FontSize',20,'FontWeight','bold','Color','g')
plot3(move4_pca.accelDataDBc(:,1),move4_pca.accelDataDBc(:,2),move4_pca.accelDataDBc(:,3),'bo')
text(mean(move4_pca.accelDataDBc(:,1))+10,mean(move4_pca.accelDataDBc(:,2)),...
    mean(move4_pca.accelDataDBc(:,3)),'Sit2','FontSize',20,'FontWeight','bold','Color','b')

end






function [] = runPLOT(blockSc)


for bi = 1:length(blockSc)
    blocki = blockSc{bi};
    figure;

    for ti = 1:3
        subplot(3,1,1)
        hold on
        xtemp = blocki{ti}(1,:);
        plot(xtemp)

        subplot(3,1,2)
        hold on
        ytemp = blocki{ti}(2,:);
        plot(ytemp)

        subplot(3,1,3)
        hold on
        ztemp = blocki{ti}(3,:);
        plot(ztemp)

    end
end



end





function [xOut, yOut , zOut] = cleanSamples(Acctable)

xOut = Acctable.epsilon{2,1}{:,8};
xOut = xOut(~isnan(xOut));
yOut = Acctable.epsilon{2,1}{:,9};
yOut = yOut(~isnan(yOut));
zOut = Acctable.epsilon{2,1}{:,10};
zOut = zOut(~isnan(zOut));

end

function [pcaOUTdata] = getPCAloadings(allBLOCKs , movementNUM)
% Assuming accelData is your matrix with three columns (X, Y, Z)

% Convert to X , Y , Z - with 3 trials along rows

tempBlOCK = allBLOCKs(:,movementNUM);

Xaccel = reshape(transpose(tempBlOCK{1}),numel(tempBlOCK{1}),1);
Yaccel = reshape(transpose(tempBlOCK{2}),numel(tempBlOCK{2}),1);
Zaccel = reshape(transpose(tempBlOCK{3}),numel(tempBlOCK{3}),1);

accelData = [Xaccel , Yaccel , Zaccel];

pdistAll = zeros(height(accelData),1);
for pi = 1:height(accelData)

    D = pdist2(accelData(pi,:),accelData);

    minD = min(D(D ~= 0));
    pdistAll(pi) = minD;
end

pdistALLs = sort(pdistAll,'descend');
plot(pdistALLs)

%%%% THIS WORKS for 1 and 2
epsilon = mean(pdistALLs) + (std(pdistALLs));

% close all
idx_DBSCAN = dbscan(accelData,epsilon,15); 
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

accelDATA_dbscan = accelData(idx_DBSCAN == maxID,:);
accelData_standardized = zscore(accelDATA_dbscan);

% Perform PCA
[coeff,score,~,~,explained] = pca(accelData_standardized);

pcaOUTdata.movNum = movementNUM;
pcaOUTdata.accelData = accelData;
pcaOUTdata.accelDataDBc = accelDATA_dbscan;
pcaOUTdata.loading = coeff;
pcaOUTdata.score = score;
pcaOUTdata.explained = explained;

% The output 'coeff' contains the principal component vectors (loadings).

end






function [pcaCOMPaRE] = comparePCAS(move1 , move2)


% Assuming loadings1 and loadings2 are the loadings of your two PCA datasets



% Perform Procrustes Analysis
[pcaCOMPaRE, ~, ~] = procrustes(move1.loading, move2.loading);

% d is the disparity between loadings1 and the transformed loadings2.
% Z is loadings2 transformed to best match loadings1.
% transform is a structure describing the optimal transformation.



% Assuming accelData is your matrix with three columns (X, Y, Z)




end



function [] = MDS_plot(moveIN)

% Compute the Euclidean distance matrix
D = pdist(moveIN.accelData, 'euclidean');

% Apply MDS
% [Y, stress] = mdscale(D, 2);
[Y, ~] = mdscale(D, 2);
% The output 'Y' is a configuration of points in a two-dimensional space.
% 'stress' is Kruskal's stress criterion that measures goodness of fit.

% Plot the result
figure;
scatter(Y(:,1), Y(:,2));
xlabel('First Dimension');
ylabel('Second Dimension');
title('Multidimensional Scaling (MDS) of Accelerometry Data');






end