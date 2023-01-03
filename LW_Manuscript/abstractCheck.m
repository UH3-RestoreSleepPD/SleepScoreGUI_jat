% Abstract analysis

cd('E:\Dropbox\Publications_Meta\InProgress\LWest_ScoreConsensus2022\Extra material\StatsAnalysis')
load('Final_InitialAgreementRaw.mat')

inSumAll = [];
for ini = 1:length(inSum)
    inSumAll = [inSumAll ; inSum{ini}];
end



inSumAlluFi = zeros(length(inSum),1);
for iniS = 1:length(inSum)
    inSumAlluFi(iniS) = 1 - (sum(matches(inSum{iniS},'U')) / length(inSum{iniS}));
end

mean(inSumAlluFi)

inSumAlluF = 1 - (sum(matches(inSumAll,'U')) / length(inSumAll));

% Get overall Final Fractions of Sleep states
finSumAll = [];
for fni = 1:length(finSum)
    finSumAll = [finSumAll ; finSum{fni}];
end


finSumAlluFi = zeros(length(inSum),1);
for iniS = 1:length(inSum)
    finSumAlluFi(iniS) = 1 - (sum(matches(finSum{iniS},'U')) / length(finSum{iniS}));
end

mean(finSumAlluFi)

finSumAlluF = 1 - (sum(matches(finSumAll,'U')) / length(finSumAll));






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
iniSSf = zeros(length(intSSu),1);
iniSSc = zeros(length(intSSu),1);
getUlocs = matches(inSumAllsrt,'U');
for iniS = 1:length(intSSu)

    iniSSf(iniS) = sum(matches(finSumAllsrt(getUlocs),intSSu{iniS})) / sum(getUlocs);
    iniSSc(iniS) = sum(matches(finSumAllsrt(getUlocs),intSSu{iniS}));

end


intSSu = {'W';'N1';'N2';'N3';'R'};
% Get overall Final Fractions of Sleep states
ifiniSSc = zeros(length(intSSu),2);
for iniS = 1:length(intSSu)

    ifiniSSc(iniS,1) = sum(matches(finSumAllsrt,intSSu{iniS}));
    ifiniSSc(iniS,2) = sum(matches(inSumAllsrt,intSSu{iniS}));

end










