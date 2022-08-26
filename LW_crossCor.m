function [] = LW_crossCor()
%UNTITLED3 Summary of this function goes here


cd('E:\Dropbox\Publications_Meta\InProgress\LWest_ScoreConsensus2022\Extra material\StatsAnalysis')

% Load Initial data
load("InitialReview.mat","initialDat");


allCell = {};
for ii = 1:length(initialDat)

    tmpI = initialDat{ii};

    tmpCell = table2cell(tmpI(:,1:4));

    missInd = matches(tmpCell,'') | cellfun(@(x) isfloat(x), tmpCell);

    tmpCell(missInd) = {'U'};

    allCell = [allCell ; tmpCell];

end

% Reclassify
corMat = zeros(size(allCell));

letMap = {'U','W','N1','N2','N3','R'};
numMap = [0, 1, 2, 3, 4, 5];

for ui = 1:length(letMap)

    tmpUni = letMap{ui};
    reMap = matches(allCell,tmpUni);
    corMat(reMap) = numMap(ui);

end

outtabe = array2table(corMat,"VariableNames",{'LW','CK','MS','ST'});

writetable(outtabe,'corMAT.csv')


end