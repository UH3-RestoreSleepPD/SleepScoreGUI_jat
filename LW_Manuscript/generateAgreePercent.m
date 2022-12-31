function [] = generateAgreePercent()
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here

% Compute fraction agreement by date

curname = getenv('COMPUTERNAME');

switch curname
    case 'DESKTOP-FAGRV5G' % home pc

        maindir.data = ['E:\Dropbox\Publications_Meta\InProgress\' ...
            'LWest_ScoreConsensus2022\Extra material\StatsAnalysis'];
        maindir.save = ['E:\Dropbox\Publications_Meta\InProgress\' ...
            'LWest_ScoreConsensus2022\Extra material\StatsAnalysis'];

    case 'OTHER PC' % Erin add your computer name and directory locations [data is on box]

end


cd(maindir.data)

% Initial
load("InitialReview.mat","allLIST","initialDat");

initialList = allLIST;
initialData = initialDat;

% initialSummary = %
inSum = cell(size(initialData));
for dayII = 1:width(initialData)

    tmpIday = initialData{dayII};

    tmpIdata = cell(height(tmpIday),1);
    for epochII = 1:height(tmpIday)

        tmpEpoch = tabulate(table2cell(tmpIday(epochII,:)));
        numEpochs = sum(cellfun(@(x) ~isempty(x), table2cell(tmpIday(epochII,:))));

        if numEpochs == 0
            tmpIdata{epochII,1} = 'U';
            continue
        end

        [maxPer,maxInd] = max(cell2mat(tmpEpoch(:,3)));

        if maxPer >= 65
            tmpIdata{epochII,1} = tmpEpoch{maxInd,1};
        else
            tmpIdata{epochII,1} = 'U';
        end

    end
    inSum{1,dayII} = tmpIdata;
end

% Final
load("finalReview.mat","finalConfL","finalNlist");

finalList = finalNlist;
finalData = finalConfL;

% finalSummary
finSum = cell(size(finalData));
for dayII = 1:width(finalData)

    tmpIday = finalData{dayII};

    tmpIdata = cell(height(tmpIday),1);
    for epochII = 1:height(tmpIday)

        numEpochs = sum(cellfun(@(x) ~isempty(x), table2cell(tmpIday(epochII,:))));
        if numEpochs == 0 || numEpochs == 1
            tmpIdata{epochII,1} = 'U';
            continue
        end

        tmpEpoch = tabulate(table2cell(tmpIday(epochII,:)));
        [maxPer,maxInd] = max(cell2mat(tmpEpoch(:,3)));

        if maxPer >= 65
            tmpIdata{epochII,1} = tmpEpoch{maxInd,1};
        else
            tmpIdata{epochII,1} = 'U';
        end

    end
    finSum{1,dayII} = tmpIdata;
end

% Save out ScoreSummary
cd(maindir.save)
save("Final_InitialAgreementRaw.mat","finSum","inSum","finalList","initialList")


% Generate Percent Agreement

% Initial
sleepStage = {'W','N1','N2','N3','R','U'};
istage = {};
ifile = {};
icount = [];
ipercent = [];
cCount = 1;
for dayC = 1:width(inSum)

    tmpIN = inSum{dayC};
    tmpTAB = tabulate(tmpIN);

    for sii = 1:length(sleepStage)
        tmpSS = sleepStage{sii};
        if sum(matches(tmpTAB(:,1),tmpSS)) == 0
            istage{cCount,1} = tmpSS;
            ifile{cCount,1} = initialList{dayC};
            icount(cCount,1) = 0;
            ipercent(cCount,1) = 0;
        else
            rowIND = matches(tmpTAB(:,1),tmpSS);
            istage{cCount,1} = tmpSS;
            ifile{cCount,1} = initialList{dayC};
            icount(cCount,1) = tmpTAB{rowIND,2};
            ipercent(cCount,1) = round(tmpTAB{rowIND,3},2);
        end
        cCount = cCount + 1;
    end
end

timepI = repmat({'I'},size(ipercent));

% Initial
fstage = {};
ffile = {};
fcount = [];
fpercent = [];
cCount = 1;
for dayF = 1:width(finSum)

    tmpIN = finSum{dayF};
    tmpTAB = tabulate(tmpIN);

    for sii = 1:length(sleepStage)
        tmpSS = sleepStage{sii};
        if sum(matches(tmpTAB(:,1),tmpSS)) == 0
            fstage{cCount,1} = tmpSS;
            ffile{cCount,1} = finalList{dayF};
            fcount(cCount,1) = 0;
            fpercent(cCount,1) = 0;
        else
            rowIND = matches(tmpTAB(:,1),tmpSS);
            fstage{cCount,1} = tmpSS;
            ffile{cCount,1} = finalList{dayF};
            fcount(cCount,1) = tmpTAB{rowIND,2};
            fpercent(cCount,1) = round(tmpTAB{rowIND,3},2);
        end
        cCount = cCount + 1;
    end
end

timepF = repmat({'F'},size(fpercent));

allstage = [istage ; fstage];
allfile = [ifile ; ffile];
allcount = [icount ; fcount];
allpercent = [ipercent ; fpercent];
alltime = [timepI ; timepF];

sumAgreeTab = table(allstage,allcount,allpercent,allfile,alltime,'VariableNames',...
    {'Stage','Count','Percent','CaseID','StageID'});

% Save out Sum Agree ScoreSummary
cd(maindir.save)
save("AllAgreeSummary.mat","sumAgreeTab")
writetable(sumAgreeTab,"AllAgreeSummary.csv")














end