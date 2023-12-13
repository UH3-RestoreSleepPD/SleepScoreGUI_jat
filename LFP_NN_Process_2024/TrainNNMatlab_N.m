%Run and evaluate Tensorflow model

%%%%%%%%%%%%%%%% First get the new data
%addpath('/Users/joelzylberberg/Dropbox/DBS_2022/New2022Data')

%load LFP data and sleep stages
load LFPPow_Bestbipolars_StanPennUNMC.mat

%for each patient and night, normalize (z-score) the LFP and pass it through the
%model. 
x = []
t = []

%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Then get the old data
clear LFPPow SleepStages
addpath('/Users/joelzylberberg/Dropbox/DBS_2022/OldUMinData')
load LFPPow_new
load OLDSleepStages
for pt=1:10
    %get LFP and z-score it
    lff = LFPPow{pt};
    numepochs = length(lff);
    lff = lff - mean(lff); 
    lff = lff./repmat(std(lff),numepochs,1);
    x = [x; lff];

    %clean up the labels 
    ss = SleepStages{pt};
    sst = ss;
    %sst(ss==2) = 1; 
    % sst(ss==3) = 2;
    sst(ss==5) = 3;
    sst(isnan(ss))=nan;
    sst = sst + 1;

    %sst(ss==-1,1) = 1; %AWAKE
    %sst(ss==0,1) = 1;  %AWAKE
    %sst(ss==1,2) = 1;  %nREM 1/2/3
    %sst(ss==2,2) = 1;  %nREM 1/2/3
    %sst(ss==3,2) = 1;  %nREM 1/2/3
    %sst(ss==5,3) = 1;  %REM 

    t = [t ; sst];

    if length(lff) ~= length(sst)
        pt  
    end


end

%%%%%%%%%%%%%%%% 

x = x(~isnan(t),:);
t = categorical(t(~isnan(t)));
picker = randperm(numel(t));
numtrain = ceil(0.95*numel(t)); %95% training data

trainx = x(picker(1:numtrain),:);
trainy = t(picker(1:numtrain));

%define network
layers = [
    featureInputLayer(7)
    fullyConnectedLayer(1028)
    batchNormalizationLayer
    reluLayer
    dropoutLayer(0.5)
    %fullyConnectedLayer(128)
    %batchNormalizationLayer
    %reluLayer   
    fullyConnectedLayer(1028)
    batchNormalizationLayer
    reluLayer   
    dropoutLayer(0.5)
    fullyConnectedLayer(1028)
    batchNormalizationLayer
    reluLayer
    dropoutLayer(0.5)
    fullyConnectedLayer(4)
    softmaxLayer
    classificationLayer];


options = trainingOptions('sgdm', ...
    'InitialLearnRate',0.05, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropFactor',0.5, ...
    'LearnRateDropPeriod',50, ...
    'L2Regularization', 2.0000e-03, ...
    'MaxEpochs',500, ...
    'Shuffle','every-epoch', ...
    'Verbose',false, ...
    'Plots','training-progress');

    %'LearnRateSchedule','piecewise', ...
    %'LearnRateDropFactor',0.5, ...
    %'LearnRateDropPeriod',50, ...



net = trainNetwork(trainx,trainy,layers,options);
ypred = classify(net,x(picker(numtrain+1:end),:));
accuracy = sum(ypred == t(picker(numtrain+1:end)))/numel((picker(numtrain+1:end)))

save 4catnet_new net