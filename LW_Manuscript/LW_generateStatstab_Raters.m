function [dataTable] = LW_generateStatstab_Raters()

% 3way anova by rater and timepoint for stage count

cd('E:\Dropbox\Publications_Meta\InProgress\LWest_ScoreConsensus2022\Extra material\StatsAnalysis')
% cd('F:\01_Coding_Datasets\LW_ConsensusStudy\SummaryAnalysis_12_2022\AllData')

% Load Initial data
load("InitialReview.mat","allLIST","initialDat");

initTables = initialDat;
initLIST = allLIST;

% Load Final data
load("finalReview.mat","finalNlist","finalConfL");

finalTables = finalConfL;
finalLIST = finalNlist;

% Create Initial missing List
[allMissName,allMissIND] = generateIMiss(initLIST,finalLIST,initTables,finalTables);

fNamei = {};
raterNi = {};
slStagei = {};
timePi = {};
counTi = [];

% INITIAL


for ii = 1:length(initLIST)

    tmpTab = initTables{ii};
    % Fix none character values

    % Loop through rater
    rNms = tmpTab.Properties.VariableNames;
    for ri = 1:length(rNms)

        tmpCounts = tmpTab.(rNms{ri});

        % Save index of missing
        searchID = [rNms{ri},'=',initLIST{ii}];
        matLOC = matches(allMissName,searchID);
        SmissVEC = allMissIND{matLOC};
 
        missLOG = ones(size(tmpCounts),'logical');
        missLOG(SmissVEC) = false; 

        newTab = tmpCounts(missLOG);

        tabCounts = tabulate(newTab);
        tmpTP = repmat({'I'},height(tabCounts),1);
        tmpRat = repmat(rNms(ri),height(tabCounts),1);
        tmpFN = repmat(initLIST(ii),height(tabCounts),1);

        fNamei = [fNamei ; tmpFN];
        raterNi = [raterNi ; tmpRat];
        timePi = [timePi ; tmpTP];
        slStagei = [slStagei ; tabCounts(:,1)];
        counTi = [counTi ; cell2mat(tabCounts(:,2))];

    end

end


% FINAL

fNamef = {};
raterNf = {};
slStagef = {};
timePf = {};
counTf = [];

for fi = 1:length(finalLIST)

    tmpTab = finalTables{fi};

    % Fix none character values

    % Loop through rater
    rNms = tmpTab.Properties.VariableNames;
    for ri = 1:length(rNms)

        tmpCounts = tmpTab.(rNms{ri});

        % Save index of missing
        searchID = [rNms{ri},'=',finalLIST{fi}];
        matLOC = matches(allMissName,searchID);
        SmissVEC = allMissIND{matLOC};
 
        missLOG = ones(size(tmpCounts),'logical');
        missLOG(SmissVEC) = false; 

        % Remove index of missing
        newTab = tmpCounts(missLOG);

        tabCounts = tabulate(newTab);
        tmpTP = repmat({'F'},height(tabCounts),1);
        tmpRat = repmat(rNms(ri),height(tabCounts),1);
        tmpFN = repmat(finalLIST(fi),height(tabCounts),1);

        fNamef = [fNamef ; tmpFN];
        raterNf = [raterNf ; tmpRat];
        timePf = [timePf ; tmpTP];
        slStagef = [slStagef ; tabCounts(:,1)];
        counTf = [counTf ; cell2mat(tabCounts(:,2))];

    end

end

allfName = [fNamei ; fNamef];
allrater = [raterNi ; raterNf];
alltimeP = [timePi ; timePf];
allSStage = [slStagei ; slStagef];
allcountS = [counTi ; counTf];

dataTable = table(allfName,allrater,alltimeP,allSStage,allcountS,...
    'VariableNames',{'FileN','RaterID','TimeP','SStage','CountS'});

writetable(dataTable,'SummaryTable.csv')
save('summaryTable.mat',"dataTable")

end





function [missName,missIND] = generateIMiss(inLIST,fnLIST,intab,fntab)

missIND = cell(length(inLIST),4);
missName = cell(length(inLIST),4);

for ii = 1:length(inLIST)

    tmpTab = intab{ii};

    % Loop through rater
    rNms = tmpTab.Properties.VariableNames;
    for ri = 1:length(rNms)

        tmpCounts = tmpTab.(rNms{ri});
        % Save index of missing
        missInd = find(matches(tmpCounts,''));
        missNme = [rNms{ri},'=',inLIST{ii}];

        switch rNms{ri}
            case 'LW'

                missIND{ii,1} = missInd;
                missName{ii,1} = missNme;

            case 'ST'

                missIND{ii,2} = missInd;
                missName{ii,2} = missNme;

            case 'MS'

                missIND{ii,3} = missInd;
                missName{ii,3} = missNme;

            case 'CK'

                missIND{ii,4} = missInd;
                missName{ii,4} = missNme;

        end

    end

end


for fi = 1:length(fnLIST)

    tmpTab = fntab{fi};

    % Loop through rater
    rNms = tmpTab.Properties.VariableNames;
    for ri = 1:length(rNms)

        tmpCounts = tmpTab.(rNms{ri});
        % Save index of missing
        if sum(cellfun(@(x) isfloat(x), tmpCounts)) ~= 0
            missIndf2 = find(cellfun(@(x) isfloat(x), tmpCounts));
            missIndf1 = [];
        else
            missIndf2 = [];
            missIndf1 = find(matches(tmpCounts,''));
        end
        missInd = [missIndf1 ; missIndf2];
        searchID = [rNms{ri},'=',fnLIST{fi}];
        matLOC = matches(missName,searchID);
        missVEC = missIND{matLOC};

        allMiss1 = [missVEC ; missInd];
        allMiss2 = unique(allMiss1);

        missIND{matLOC} = allMiss2;

    end

end



end


