% Rater stats

% 3way anova by rater and timepoint for stage count

cd('J:\01_Coding_Datasets\LW_ConsensusStudy\SummaryF')

% Load Initial data
load("InitialReview.mat","allLIST","initialDat");

initTables = initialDat;
initLIST = allLIST;

% Load Final data
load("finalReview.mat","finalNlist","finalConfL");

finalTables = finalConfL;
finalLIST = finalNlist;

% Extract stage counts
% Date 

fNamei = {};
raterNi = {};
slStagei = {};
timePi = {};
counTi = [];

% INITIAL

indexMiss_I = cell(length(initLIST),4);
nameMiss_I = cell(length(initLIST),4);

for ii = 1:length(initLIST)

    tmpFn = initLIST{ii};
    tmpTab = initTables{ii};

    % Fix none character values

    % Loop through rater
    rNms = tmpTab.Properties.VariableNames;
    for ri = 1:length(rNms)

        tmpCounts = tmpTab.(rNms{ri});

        % Save index of missing
        missInd = find(matches(tmpCounts,''));
        missNme = [rNms{ri},'=',initLIST{ii}];
        indexMiss_I{ii,ri} = missInd;
        nameMiss_I{ii,ri} = missNme;
        % Remove index of missing
        tmpConMr = tmpCounts(~matches(tmpCounts,''));

        tabCounts = tabulate(tmpConMr);
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

    tmpFn = finalLIST{fi};
    tmpTab = finalTables{fi};

    % Fix none character values

    % Loop through rater
    rNms = tmpTab.Properties.VariableNames;
    for ri = 1:length(rNms)

        tmpCounts = tmpTab.(rNms{ri});

        % Save index of missing
        missLOG = ones(size(tmpCounts),'logical');
        % Find missIndex
        searchID = [rNms{ri},'=',finalLIST{fi}];
        matLOC = matches(nameMiss_I,searchID);
        missVEC = indexMiss_I{matLOC};
        missLOG(missVEC) = false; 

        % Remove index of missing
        newTab = tmpCounts(missLOG);

        if sum(cellfun(@(x) isfloat(x), newTab)) ~= 0
            newTab = newTab(cellfun(@(x) ~isfloat(x), newTab));
        end

        tabCounts = tabulate(newTab);
        tmpTP = repmat({'F'},height(tabCounts),1);
        tmpRat = repmat(rNms(ri),height(tabCounts),1);
        tmpFN = repmat(initLIST(fi),height(tabCounts),1);

        fNamef = [fNamef ; tmpFN];
        raterNf = [raterNf ; tmpRat];
        timePf = [timePf ; tmpTP];
        slStagef = [slStagef ; tabCounts(:,1)];
        counTf = [counTf ; cell2mat(tabCounts(:,2))];

    end

end

test = 1;

