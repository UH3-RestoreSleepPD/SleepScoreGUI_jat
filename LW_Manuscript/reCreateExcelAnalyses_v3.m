% cd('D:\Dropbox\Publications_Meta\InProgress\LWest_ScoreConsensus2022\Extra material')
cd('E:\Dropbox\Publications_Meta\InProgress\LWest_ScoreConsensus2022\Extra material\StatsAnalysis\')
%%

% dataTab = readtable('DataForAnalysis.xlsx');
dataTab = readtable('AllAgreeSummary.csv');

%% 

uniqueNIGHTs = unique(dataTab.CaseID);

nightPer = zeros(length(uniqueNIGHTs),2);
for ui = 1:length(uniqueNIGHTs)

    for ic = 1:2

        switch ic
            case 1
                tmpTAB = dataTab(matches(dataTab.StageID,'I'),:);
            case 2
                tmpTAB = dataTab(matches(dataTab.StageID,'F'),:);
        end
        tmpUI = uniqueNIGHTs{ui};
        tmpSt = tmpTAB(matches(tmpTAB.CaseID,tmpUI),:);
        nightPer(ui,ic) = sum(tmpSt.Percent(~matches(tmpSt.Stage,{'U'})));

    end

end

nPerTab = table(uniqueNIGHTs,nightPer(:,1),nightPer(:,2),'VariableNames',...
    {'Nights','Initial','Final'});

% Determine total percent success in 2nd batch (shorter process subjects)
% i.e., excluding UNMC2 and UNMC3

grpIDs = split(nPerTab.Nights,'_');
subN = grpIDs(:,1);
instN = grpIDs(:,2);
grpNew = cellfun(@(x,y) [x,'_',y], subN,instN,'UniformOutput', false);

abbreGrp = nPerTab(~matches(grpNew,{'2_UNMC','3_UNMC'}),:);

meanABgp = mean(abbreGrp.Final);
stdABgp = std(abbreGrp.Final);

allSubs = zeros(height(nPerTab),1);
allNs = zeros(height(nPerTab),1);
allInts = cell(height(nPerTab),1);
for nni = 1:height(nPerTab)

    tmpROW = nPerTab.Nights{nni};

    tmpItems = split(tmpROW,{'_','.'});

    allSubs(nni) = str2double(tmpItems{1});
    allNs(nni) = str2double(tmpItems{3});
    allInts{nni} = tmpItems{2};
end

nPerTab.Sub = allSubs;
nPerTab.nightN = allNs;
nPerTab.Instu = allInts;


% uniSub = unique(nPerTab.Sub);
% for uui = 1:length(uniSub)
% 
%     tmpS = matches(nPerTab.Sub,uniSub(uui));
%     nPerTab.SubI(tmpS) = uui;
% 
% end

%%

for yi = 1:height(nPerTab)
    nPerTab.Ylab{yi} = ['S',num2str(nPerTab.Sub(yi)),'_',...
        nPerTab.Instu{yi},'_N',num2str(nPerTab.nightN(yi))];
end

% % 3 Gradient for UNMC
% unmcGrad = [250, 172, 183; 
%             228, 125, 130;   
%             207, 78, 78];
unmcAll = [255, 190, 11]/255;
% % 3 Gradient for Stan
% stanGrad = [172, 250, 183; 
%             117, 215, 134;   
%             63, 181, 85];
stanAll = [255, 0, 110]/255;
% % 3 Gradient for PENN
% upennGrad = [212, 172, 251; 
%             198, 125, 228;   
%             184, 79, 206];
upennAll = [58, 134, 255]/255;

% dumbbell plot
sNperTab = sortrows(nPerTab,2,'ascend');
scYv = transpose(height(sNperTab):-1:1);

coloRMap = zeros(height(sNperTab),3);
for coi = 1:height(sNperTab)

    tmpInst = sNperTab.Instu{coi};

    switch tmpInst
        case 'UNMC'
            %             tmpNight = sNperTab.nightN(coi);
            %             coloRMap(coi,:) = unmcGrad(tmpNight,:);
            coloRMap(coi,:) = unmcAll;
        case 'UPENN'
            %             tmpNight = sNperTab.nightN(coi);
            %             coloRMap(coi,:) = upennGrad(tmpNight,:);
            coloRMap(coi,:) = stanAll;
        case 'STAN'
            %             tmpNight = sNperTab.nightN(coi);
            %             coloRMap(coi,:) = stanGrad(tmpNight,:);
            coloRMap(coi,:) = upennAll;
    end
end

% coloRMapRGB = coloRMap/255;
xVals = [transpose(sNperTab.Initial) ; transpose(sNperTab.Final)];
yVals = [transpose(scYv) ; transpose(scYv)];
line(xVals,yVals,'Color','k')
hold on
scatter(sNperTab.Initial,scYv,40,coloRMap,'filled')
% scatter(sNperTab.Initial,scYv,40,"red",'filled')
% scatter(sNperTab.Initial,scYv,40,"red",'filled')
scatter(sNperTab.Final,scYv,40,"black",'filled')

yticks(1:height(sNperTab))

subID = extractBefore(sNperTab.Ylab,'_');
nightID = extractAfter(extractAfter(sNperTab.Ylab,'_'),'_');

sNperTab.Ylab2 = cellfun(@(x,y) [x , ' ' , y], subID, nightID, 'UniformOutput', false);

yticklabels(flipud(sNperTab.Ylab2));

ax = gca;
ax.TickLabelInterpreter = 'none';

ylim([1 46])
xlim([60 100])
xticks([60 70 80 90 100])
xlabel('Percent of consensus epochs')
ylabel('Subjects and Nights')
axis square

% % Break up subject and night
% for sni = 1:height(nPerTab)
%     tmpRow = nPerTab.Nights{sni};
%     tmpParts = split(tmpRow,{'_','.'});
%     patNUM = [tmpParts{1},tmpParts{2}];
%     nightNUM = str2double(tmpParts{3});
%     nPerTab.Sub{sni} = patNUM;
%     nPerTab.night(sni) = nightNUM;
% end

% Clean up original datatab


%% Transition stack plots

cd('E:\Dropbox\Publications_Meta\InProgress\LWest_ScoreConsensus2022\Extra material\StatsAnalysis')
load('Final_InitialAgreementRaw.mat','finSum','finalList','inSum','initialList')

% Get overall Initial Fractions of Sleep states
inSumAll = [];
for ini = 1:length(inSum)
    inSumAll = [inSumAll ; inSum{ini}];
end
intSSu = {'W';'N1';'N2';'N3';'R';'U'};
% Get overall Final Fractions of Sleep states
iniSSf = zeros(length(intSSu),1);
for iniS = 1:length(intSSu)
    iniSSf(iniS) = sum(matches(inSumAll,intSSu{iniS})) / length(inSumAll);
end

% Get overall Final Fractions of Sleep states
finSumAll = [];
for fni = 1:length(finSum)
    finSumAll = [finSumAll ; finSum{fni}];
end
fintSSu = {'W';'N1';'N2';'N3';'R';'U'};
% Get overall Final Fractions of Sleep states
finiSSf = zeros(length(fintSSu),1);
for finiS = 1:length(fintSSu)
    finiSSf(finiS) = sum(matches(finSumAll,fintSSu{finiS})) / length(finSumAll);
end

colorMap = [255, 200, 44;
    255, 127, 5;
    207, 64, 47;
    170, 32, 76;
    96, 0, 129;
    0, 0, 0];
colorMPrgb = colorMap/255;

% ALL DATA PLOT
figure;
y = [transpose(iniSSf) ; transpose(finiSSf)];
b = bar(y,"stacked");
xticklabels(["Initial","Final"]);
% b.CData(2,:) = [.5 0 .5];
for bi = 1:6
    b(bi).FaceAlpha = 0.3;
    b(bi).FaceColor = colorMPrgb(bi,:);
end
b(1).BarWidth = 0.45;

iniallBs = zeros(length(iniSSf),3);
for si = 1:length(iniSSf)
    if si == 1
        iniallBs(si,1) = 0;
        iniallBs(si,3) = iniSSf(si);
        iniallBs(si,2) = ((iniallBs(si,3) - iniallBs(si,1))/2) + iniallBs(si,1);
    elseif si == length(iniSSf)
        iniallBs(si,1) = iniallBs(si-1,1) + iniSSf(si - 1);
        iniallBs(si,3) = 1;
        iniallBs(si,2) = ((iniallBs(si,3) - iniallBs(si,1))/2) + iniallBs(si,1);
    else
        iniallBs(si,1) = iniallBs(si-1,1) + iniSSf(si - 1);
        iniallBs(si,3) = iniallBs(si-1,3) + iniSSf(si);
        iniallBs(si,2) = ((iniallBs(si,3) - iniallBs(si,1))/2) + iniallBs(si,1);
    end
end

hold on

for li = 1:6

    line([0.5 1],[iniallBs(li,2) iniallBs(li,2)],'Color',colorMPrgb(li,:),'LineWidth',1)
    text(0.48,iniallBs(li,2),[intSSu{li},': ' , num2str(round(iniSSf(li)*100,1)),'%'],...
        'HorizontalAlignment','right','VerticalAlignment','middle',...
        'Color',colorMPrgb(li,:),'FontWeight','bold')
end


finiallBs = zeros(length(finiSSf),3);
for si = 1:length(finiSSf)
    if si == 1
        finiallBs(si,1) = 0;
        finiallBs(si,3) = finiSSf(si);
        finiallBs(si,2) = ((finiallBs(si,3) - finiallBs(si,1))/2) + finiallBs(si,1);
    elseif si == length(finiSSf)
        finiallBs(si,1) = finiallBs(si-1,1) + finiSSf(si - 1);
        finiallBs(si,3) = 1;
        finiallBs(si,2) = ((finiallBs(si,3) - finiallBs(si,1))/2) + finiallBs(si,1);
    else
        finiallBs(si,1) = finiallBs(si-1,1) + finiSSf(si - 1);
        finiallBs(si,3) = finiallBs(si-1,3) + finiSSf(si);
        finiallBs(si,2) = ((finiallBs(si,3) - finiallBs(si,1))/2) + finiallBs(si,1);
    end
end

for li = 1:6

    line([2 2.5],[finiallBs(li,2) finiallBs(li,2)],'Color',colorMPrgb(li,:),'LineWidth',1)
    text(2.52,finiallBs(li,2),[fintSSu{li},': ' , num2str(round(finiSSf(li)*100,1)),'%'],...
        'HorizontalAlignment','left','VerticalAlignment','middle',...
        'Color',colorMPrgb(li,:),'FontWeight','bold')
end

% CREATE SORTED FILES
% Loop through each Stage in inSumS
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

% Individual subjects plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%figure
figure;
% Subject 3 and 6
% Get Bar data 
[iniBar_s3 , finBar_s3, iniC_s3, finC_s3] = getIndivBarData(3 , inLS, inSumS , finSumS , fintSSu);

t1 = tiledlayout(1,2);

nexttile

bi3 = bar(iniBar_s3,"stacked");
xticklabels(["Night 1","Night 2", "Night 3"]);
% b.CData(2,:) = [.5 0 .5];
for bi = 1:6
    bi3(bi).FaceAlpha = 0.5;
    bi3(bi).FaceColor = colorMPrgb(bi,:);
    bi3(bi).EdgeColor = 'none';
end
bi3(1).BarWidth = 1;
xline(1.5)
xline(2.5)
subtitle('Subject 3 - Initial Review')
axi3 = gca;
axi3.TitleHorizontalAlignment = 'left';
ylabel('Fraction of sleep stage')
yticks([0 0.25 0.5 0.75 1])
axis square

nexttile

bf3 = bar(finBar_s3,"stacked");
xticklabels(["Night 1","Night 2", "Night 3"]);
% b.CData(2,:) = [.5 0 .5];
for bi = 1:6
    bf3(bi).FaceAlpha = 0.5;
    bf3(bi).FaceColor = colorMPrgb(bi,:);
    bf3(bi).EdgeColor = 'none';
end
bf3(1).BarWidth = 1;
xline(1.5)
xline(2.5)
subtitle('Subject 3 - Final Review')
axf3 = gca;
axf3.TitleHorizontalAlignment = 'left';
ylabel('Fraction of sleep stage')
yticks([0 0.25 0.5 0.75 1])
axis square

t.TileSpacing = 'compact';
t.Padding = 'compact';

finINchange_s3 = finC_s3 - iniC_s3;
generateIndivEpochChangePlot(finINchange_s3,colorMPrgb, 3)

% Group subject plots

totalCount = zeros(length(intSSu),1);
totalper = zeros(length(intSSu),1);
allBINS = cell(length(intSSu),1);
for ssi2 = 1:length(intSSu)

    % Get Index for empty epochs both F and I
    iempt = matches(inSumAllsrt,' ');
    fempt = matches(finSumAllsrt,' ');
    keepIND = ~(iempt == 1 | fempt == 1);

    inSumS2 = inSumAllsrt(keepIND);
    finSumS2 = finSumAllsrt(keepIND);

    ssIind = matches(inSumS2,intSSu{ssi2});

    % All epochs
    allssEpochs = sum(ssIind);
    totalCount(ssi2) = allssEpochs;
    totalper(ssi2) = allssEpochs/length(inSumS2);

    % All fin SS ids
    allfinIds = finSumS2(ssIind);

    uniFinds = unique(allfinIds);
    nonSSind = ~matches(uniFinds,intSSu{ssi2});
    uniFindsNON = uniFinds(nonSSind);

    finIDvec = cell(length(uniFindsNON),1);
    finIDcount = zeros(length(uniFindsNON),1);
    for ui = 1:length(uniFindsNON)

        tmpUI = sum(matches(allfinIds,uniFindsNON{ui}));
        finIDvec{ui} = uniFindsNON{ui};
        finIDcount(ui) = tmpUI;

    end

    nonSStot = sum(finIDcount);
    finIDper = finIDcount/nonSStot;

    tmpTable = table(finIDvec,finIDcount,finIDper,'VariableNames',...
        {'ChangeSS','Count','Percent'});

    allBINS{ssi2} = tmpTable;

end

% clean up bins
allBINSU = allBINS;
for ssi2 = 1:length(intSSu)

    switch intSSu{ssi2}
        case {'W','N1','N2','N3','R'}
            tmpTab2 = allBINS{ssi2};
            tmpUtab = tmpTab2(matches(tmpTab2.ChangeSS,'U'),:);
            tmpUtab.Percent = tmpUtab.Count/(totalCount(ssi2));
            allBINSU{ssi2} = tmpUtab;

        case {'U'}
            continue;
    end
end
fnonSSind = ~matches(intSSu,'U');
% Creat plot for which stages became Us
% Make the plot

% ONLY PLOT WHERE U epochs go
getUtable = allBINSU{matches(intSSu,'U')};

countSSuu = getUtable.Count;
barDATAUU = transpose(round(getUtable.Percent,2));

figure;
b2 = bar(barDATAUU);
b2.FaceColor = 'flat';

colorMPrgbN = colorMPrgb(fnonSSind,:);

for nli = 1:5

    b2.FaceAlpha = 0.3;
    b2.CData(nli,:) = colorMPrgbN(nli,:);

    tmpCOUs = countSSuu(nli);
    tmpCOUsU = [num2str(round(tmpCOUs))];

    text(nli,barDATAUU(nli),tmpCOUsU,...
        'HorizontalAlignment','center','VerticalAlignment','bottom',...
        'Color',colorMPrgbN(nli,:),'FontWeight','bold')
end

ylim([0 0.4])
yticks([0 0.1 0.2 0.3 0.4])
ylabel("Fraction of undecided converted to consensus")


xticklabels(categorical(intSSu(fnonSSind)));


% ONLY PLOT WHERE U epochs go
% Make the plot
% Single bar for each stage by count

getBardata = zeros(1,5);
allBINSUn = allBINSU(fnonSSind);
countSS = cellfun(@(x) x.Count, allBINSUn,'UniformOutput', false);
countSS{cellfun(@(x) isempty(x), countSS,'UniformOutput', true)} = 0;
countSSu = cell2mat(countSS);
perSS = cellfun(@(x) x.Percent, allBINSUn,'UniformOutput', false);
perSS{cellfun(@(x) isempty(x), perSS,'UniformOutput', true)} = 0;
perSSu = cell2mat(perSS);

figure;
barDATAnon = transpose(countSSu);
b3 = bar(barDATAnon);
b3.FaceColor = 'flat';

colorMPrgbN = colorMPrgb(fnonSSind,:);

for nli = 1:5

    b3.FaceAlpha = 0.3;
    b3.CData(nli,:) = colorMPrgbN(nli,:);

    tmpPERs = perSSu(nli);
    if tmpPERs == 0
        tmpPERsU = '0%';
    elseif tmpPERs < 0.01
        tmpPERsU = '< 0.01%';
    else
        tmpPERsU = [num2str(round(tmpPERs,2)) ,'%'];
    end

    text(nli,barDATAnon(nli),tmpPERsU,...
        'HorizontalAlignment','center','VerticalAlignment','bottom',...
        'Color',colorMPrgbN(nli,:),'FontWeight','bold')
end

xticklabels(categorical(intSSu(fnonSSind)));

ylabel("Number of epochs converted to undecided")


%% Timeline plot

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

        if sleepStC >= 5
            disp(['INDEX found! ', num2str(eii-4)])
            sleepStarInd(iii) = eii - 4;
            break
        end

    end
end


% Loop through all cases and for each - select from 5 mins minus start and
% determine how many steps have changed

for ifi = 1:length(inSumS)

    beginIND = sleepStarInd(ifi) - 10;
    iniTEMP = inSumS{ifi}(beginIND:end);
    fniTEMP = finSumS{ifi}(beginIND:end);







end




























