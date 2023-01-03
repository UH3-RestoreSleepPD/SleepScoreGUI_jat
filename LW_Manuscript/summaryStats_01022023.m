% Glm negative possion assessment

% cd('E:\Dropbox\Publications_Meta\InProgress\LWest_ScoreConsensus2022\Extra material\StatsAnalysis')
cd('D:\Dropbox\Publications_Meta\InProgress\LWest_ScoreConsensus2022\Extra material\StatsAnalysis')
load('Final_InitialAgreementRaw.mat')

[finLS,fiLi] = sort(finalList);
finSumS = finSum(fiLi);
[inLS,iLi] = sort(initialList);
inSumS = inSum(iLi);

% Find rows in finSumS and determine the fraction of endpoints
inSumAllsrt = [];
for ini = 1:length(inSumS)
    inSumAllsrt = [inSumAllsrt ; inSumS{ini}];
end

finSumAllsrt = [];
for fni = 1:length(finSumS)
    finSumAllsrt = [finSumAllsrt ; finSumS{fni}];
end


intSSu = {'W';'N1';'N2';'N3';'R'};
% Get overall Final Fractions of Sleep states

nightPredictor = nan(length(finSumS),1);
stagePredictor = nan(length(finSumS),1);
undecideCnt = nan(length(finSumS),1);
nightStageC = 1;
for iniS = 1:length(finSumS)

    % get a night
    tmpNNight = finSumS{iniS};
    % get initial U index
    tmpUUind = inSumS{iniS};
    getUlocs = matches(tmpUUind,'U');

    % get night num
    tmpNname = inLS{iniS};
    tmpNNparts = split(tmpNname,{'_','.'});
    tmpNightU = tmpNNparts{3};


    for epi = 1:length(intSSu)
        nightPredictor(nightStageC) = str2double(tmpNightU);
        stagePredictor(nightStageC) = epi;
        undecideCnt(nightStageC) = sum(matches(tmpNNight(getUlocs),intSSu{epi}));
        nightStageC = nightStageC + 1;
    end

end

%%

% Load the data
data = undecideCnt;

% Extract the count outcome and predictor variables
y = data;
x = [nightPredictor , stagePredictor];

% Fit the negative binomial regression model
% [b, dev, stats] = glmfit(x, y, 'poisson', 'log', 'on');
% fitO = glmfit(x, y, 'poisson', 'log', 'on');

mdl =  fitglm(x,data,'linear','Distribution','poisson')
confint = coefCI(mdl)
% Make predictions using the model
plotDiagnostics(mdl)
legend('show') % Show the legend
plotDiagnostics(mdl,'cookd')

plotPartialDependence(mdl,2,"Conditional","centered")
% plotPartialDependence(mdl,2,mdl.ClassNames);

imp = predictorImportance(Mdl);
figure
bar(imp)
title("Predictor Importance Estimates")
ylabel("Estimates")
xlabel("Predictors")
ax = gca;
ax.XTickLabel = Mdl.PredictorNames;


%% 

% Extract the initial counts and final counts
iempt = matches(inSumAllsrt,' ');
fempt = matches(finSumAllsrt,' ');
keepIND = ~(iempt == 1 | fempt == 1);
allstages = [inSumAllsrt(keepIND) ; finSumAllsrt(keepIND)];
conStage = [repmat({'I'},size(inSumAllsrt(keepIND))) ; repmat({'F'},size(inSumAllsrt(keepIND)))];
tablEE = table(allstages,conStage,'VariableNames',{'SleepSt','ConSt'});


% Create a contingency table of the counts
[conttbl,chi2,p,labels] = crosstab(tablEE.SleepSt,tablEE.ConSt);
heatmap(tablEE,'SleepSt','ConSt');
% % Perform the chi-square test
% [h, stats] = chi2test(ct);
% 
% % Display the p-value
% disp(p)

%% icc

% loop through subjects starting with Initial
% inSumS % finSumS
sleepStarInd = zeros(length(inSumS),1);
for iii = 1:length(inSumS)

    % Get index for the start of contiguous sleep block (5 min = 10 epochs)
    tmpINsub = inSumS{iii};

    sleepStC = 0;
    for eii = 1:length(tmpINsub)

        tmpEpochei = tmpINsub{eii};
        if matches(tmpEpochei,{'N1','N2','N3','R'})
            sleepStC = sleepStC + 1;
        else
            sleepStC = 0;
        end

        if sleepStC >= 10
            disp(['INDEX found! ', num2str(eii-9)])
            sleepStarInd(iii) = eii - 9;
            break
        end

    end
end

% cd('E:\Dropbox\Publications_Meta\InProgress\LWest_ScoreConsensus2022\Extra material\StatsAnalysis')
cd('D:\Dropbox\Publications_Meta\InProgress\LWest_ScoreConsensus2022\Extra material\StatsAnalysis')
load('InitialReview.mat')
intSSu = {'N1';'N2';'N3';'R'};
% Get overall Final Fractions of Sleep states
raterTST = zeros(length(intSSu),4,length(initialDat));
for iniS = 1:length(initialDat)
    tmpNNight = initialDat{iniS};
    tmpNNightS = tmpNNight(sleepStarInd(iniS):end,:);

    raterS = {'LW','ST','MS','CK'};
    for ri = 1:4
        tmpRATE = raterS{ri};

        tmpNtab = tabulate(tmpNNight.(tmpRATE));

        for ssi = 1:length(intSSu)

            tabROW = matches(tmpNtab(:,1),intSSu{ssi});
            totStimeI = matches(tmpNtab(:,1),{'R','N1','N2','N3'});
            totStimeE = sum(cell2mat(tmpNtab(totStimeI,2))); % in epochs
            totStimeS = totStimeE*30; % seconds

            if sum(tabROW) == 0
                continue
            else
                sleepStageE = cell2mat(tmpNtab(tabROW,2));
                sleepStageS = sleepStageE*30; % seconds
                fracTST = sleepStageS / totStimeS;
                raterTST(ssi,ri,iniS) = fracTST;
            end
        end
    end
end

meanT = mean(raterTST,3);


ICC = f_ICC(meanT,0.05);

%%

intSSu = {'W';'N1';'N2';'N3';'R'};
% Get overall Final Fractions of Sleep states


stagePredictorI = cell(length(finSumS),1);
undecideCntI = nan(length(finSumS),1);
stagePredictorF = cell(length(finSumS),1);
undecideCntF = nan(length(finSumS),1);
nightStageC = 1;
for iniS = 1:length(finSumS)

    % get a night
    tmpNNightI = inSumS{iniS};
    tmpNNightF = finSumS{iniS};

    for epi = 1:length(intSSu)
        stagePredictorI{nightStageC} = intSSu{epi};
        stagePredictorF{nightStageC} = intSSu{epi};
        undecideCntI(nightStageC) = sum(matches(tmpNNightI,intSSu{epi}));
        undecideCntF(nightStageC) = sum(matches(tmpNNightF,intSSu{epi}));
        nightStageC = nightStageC + 1;
    end

end

stageALL =  [stagePredictorI ; stagePredictorF];
assessALL = [repmat({'I'},size(undecideCntI)) ; repmat({'F'},size(undecideCntI))];
countALL = [undecideCntI ; undecideCntF];


%%

intSSu = {'W';'N1';'N2';'N3';'R'};
% Get overall Final Fractions of Sleep states


stagePredictor = cell(length(finSumS),1);
undecideCnt = nan(length(finSumS),1);
nightStageC = 1;
for iniS = 1:length(finSumS)

    % get a night
    tmpNNight = finSumS{iniS};
    % get initial U index
    tmpUUind = inSumS{iniS};
    getUlocs = matches(tmpUUind,'U');


    for epi = 1:length(intSSu)
        stagePredictor{nightStageC} = intSSu{epi};
        undecideCnt(nightStageC) = sum(matches(tmpNNight(getUlocs),intSSu{epi}));
        nightStageC = nightStageC + 1;
    end

end