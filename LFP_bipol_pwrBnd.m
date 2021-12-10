function [] = LFP_bipol_pwrBnd(folderLoc , saveLoc)

% cd('J:\01_Coding_Datasets\NeuralNetwork_Sleep\NewData\UNMC_2');
% J:\01_Coding_Datasets\NeuralNetwork_Sleep\PreFinal
cd(folderLoc)

% Get list
matDIR = dir('*.mat');
fileDIR = {matDIR.name};

sleepData = struct;
% Loop through
for fi = 1:length(fileDIR)
    cd(folderLoc)
    tmpFile = fileDIR{fi};
    % Extract relevant info [Subject, Institution, Night]
    tmpItems = split(tmpFile,{'_','.'});
    tmpSub = tmpItems{1};
    tmpNight = tmpItems{3};
    tmpInstit = tmpItems{2};

    % Save out struct with biopolar table and Norm Power table and key
    load(tmpFile,'LFPTT')
    % Extract contacts as table
    leadTable = table2array(LFPTT(:,1:4));
    % Mean subtract
    lead0 = leadTable(:,1);
    lead1 = leadTable(:,2);
    lead2 = leadTable(:,3);
    lead3 = leadTable(:,4);

    lead01 = cellfun(@(x , y) x - mean([x,y],2), lead0, lead1, 'UniformOutput',false );
    lead12 = cellfun(@(x , y) x - mean([x,y],2), lead1, lead2, 'UniformOutput',false );
    lead23 = cellfun(@(x , y) x - mean([x,y],2), lead2, lead3, 'UniformOutput',false );
    allBipol = [lead01 , lead12 , lead23];
    % % Compute normalized power

    freqIdx = [0,3; ...     %Delta
        3,8; ...            %Theta
        8,13; ...           %Alpha
        13,22; ...          %low-beta
        22,30; ...          %high-beta
        30,50; ...          %low-gamma
        50,100; ...         %mid-gamma
        100,125];           %high-gamma
    freqList = {'Delta';'Theta'; 'Alpha'; 'LowBeta'; ...
        'HighBeta'; 'LowGamma'; 'MidGamma'; 'HighGamma'};

    freqPwrAv = cell(size(allBipol));
    for bpi = 1:size(allBipol,2) % Number of bipolar leads
        for epi = 1:size(allBipol,1) % Number of epochs
            % Temporary bin
            tmpBin = allBipol{epi,bpi};
            % Get power and bins
            [pwrR , freqI] = pspectrum(tmpBin,250,'FrequencyResolution',1);
            pwRnorm = normalize(pwrR);
            freqAveTmp = zeros(8,1);
            for pbi = 1:size(freqIdx,1) % Number of frequency bins
                freqAveTmp(pbi) = mean(pwRnorm(freqI > freqIdx(pbi,1) & freqI < freqIdx(pbi, 2)));
            end
            freqPwrAv{epi,bpi} = freqAveTmp;
        end
    end

    % Convert bipolar back into Timetable with Labels
    % Add bipolar
    BPrLFP = cell2table(allBipol, 'VariableNames',{'L01','L12','L23'});
    BPrLFP.FSScore = LFPTT.FSScore;

    % Save stuff
    sleepData.bipolarLFP = BPrLFP;
    sleepData.pwrFreq = freqPwrAv;
    sleepData.info.freqIndx = freqIdx;
    sleepData.info.freqList = freqList;
    sleepData.info.sub = tmpSub;
    sleepData.info.night = tmpNight;
    sleepData.info.institute = tmpInstit;

    cd(saveLoc)
    saveName = [tmpSub , '_' , tmpInstit , '_' , tmpNight , '_fNN.mat'];
    save(saveName, 'sleepData');

end

