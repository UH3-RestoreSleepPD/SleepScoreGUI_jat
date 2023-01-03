% cd('D:\Dropbox\Publications_Meta\InProgress\LWest_ScoreConsensus2022\Extra material')
% cd('E:\Dropbox\Publications_Meta\InProgress\LWest_ScoreConsensus2022\Extra material\StatsAnalysis\')
cd('D:\Dropbox\Publications_Meta\InProgress\LWest_ScoreConsensus2022\Extra material\StatsAnalysis\')
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
            coloRMap(coi,:) = upennAll;
        case 'STAN'
            %             tmpNight = sNperTab.nightN(coi);
            %             coloRMap(coi,:) = stanGrad(tmpNight,:);
            coloRMap(coi,:) = stanAll;
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

%%
tiledlayout(3,1)
nexttile

night1tab = sNperTab(ismember(sNperTab.nightN,1),:);

coloRMapN1 = zeros(height(night1tab),3);
for coi = 1:height(night1tab)
    tmpInst = night1tab.Instu{coi};
    switch tmpInst
        case 'UNMC'
            coloRMapN1(coi,:) = unmcAll;
            night1tab.CM{coi} = 'y';
        case 'UPENN'
            coloRMapN1(coi,:) = upennAll;
            night1tab.CM{coi} = 'b';
        case 'STAN'
            coloRMapN1(coi,:) = stanAll;
            night1tab.CM{coi} = 'r';
    end
end

xValsN1 = [transpose(night1tab.Initial) ; transpose(night1tab.Final)];
yValN1a = width(xValsN1):-1:1;
yValsN1 = [yValN1a ; yValN1a];
line(xValsN1,yValsN1,'Color','k')
hold on
sctiN1 = night1tab.Initial;
sctfN1 = night1tab.Final;
sctyN1 = transpose(yValN1a);
colmpN1 = coloRMapN1;
scatter(sctiN1,sctyN1,40,colmpN1,'filled')
scatter(sctfN1,sctyN1,40,"black",'filled')
xline(mean(sctiN1),'--')
xline(mean(sctfN1),'-')

yticks(1:height(sctiN1))
yticklabels(flipud(night1tab.Nights));
subtitle('Night 1')
ax1 = gca;
ax1.TickLabelInterpreter = 'none';
ax1.TitleHorizontalAlignment = 'left';

ylim([0 max(sctyN1)+1])
xlim([60 100])
xticks([60 70 80 90 100])
xlabel('Percent of consensus epochs')
ylabel('Subjects')

axis square

nexttile % Night 2

night2tab = sNperTab(ismember(sNperTab.nightN,2),:);

coloRMapN2 = zeros(height(night2tab),3);
for coi = 1:height(night2tab)
    tmpInst = night2tab.Instu{coi};
    switch tmpInst
        case 'UNMC'
            coloRMapN2(coi,:) = unmcAll;
            night2tab.CM{coi} = 'y';
        case 'UPENN'
            coloRMapN2(coi,:) = upennAll;
            night2tab.CM{coi} = 'b';
        case 'STAN'
            coloRMapN2(coi,:) = stanAll;
            night2tab.CM{coi} = 'r';
    end
end

xValsN2 = [transpose(night2tab.Initial) ; transpose(night2tab.Final)];
yValN2a = width(xValsN2):-1:1;
yValsN2 = [yValN2a ; yValN2a];
line(xValsN2,yValsN2,'Color','k')
hold on
sctiN2 = night2tab.Initial;
sctfN2 = night2tab.Final;
sctyN2 = transpose(yValN2a);
scatter(sctiN2,sctyN2,40,coloRMapN2,'filled')
scatter(sctfN2,sctyN2,40,"black",'filled')
xline(mean(sctiN2),'--')
xline(mean(sctfN2),'-')

yticks(1:height(sctiN2))
yticklabels(flipud(night2tab.Nights));
subtitle('Night 2')
ax1 = gca;
ax1.TickLabelInterpreter = 'none';
ax1.TitleHorizontalAlignment = 'left';

ylim([0 max(sctyN2)+1])
xlim([60 100])
xticks([60 70 80 90 100])
xlabel('Percent of consensus epochs')
ylabel('Subjects')

axis square

nexttile % Night 3

night3tab = sNperTab(ismember(sNperTab.nightN,3),:);

coloRMapN3 = zeros(height(night3tab),3);
for coi = 1:height(night3tab)
    tmpInst = night3tab.Instu{coi};
    switch tmpInst
        case 'UNMC'
            coloRMapN3(coi,:) = unmcAll;
            night3tab.CM{coi} = 'y';
        case 'UPENN'
            coloRMapN3(coi,:) = upennAll;
            night3tab.CM{coi} = 'b';
        case 'STAN'
            coloRMapN3(coi,:) = stanAll;
            night3tab.CM{coi} = 'r';
    end
end

xValsN3 = [transpose(night3tab.Initial) ; transpose(night3tab.Final)];
yValN3a = width(xValsN3):-1:1;
yValsN3 = [yValN3a ; yValN3a];
line(xValsN3,yValsN3,'Color','k')
hold on
sctiN3 = night3tab.Initial;
sctfN3 = night3tab.Final;
sctyN3 = transpose(yValN3a);
scatter(sctiN3,sctyN3,40,coloRMapN3,'filled')
scatter(sctfN3,sctyN3,40,"black",'filled')
xline(mean(sctiN3),'--')
xline(mean(sctfN3),'-')

yticks(1:height(sctiN3))
yticklabels(flipud(night3tab.Nights));
subtitle('Night 3')
ax3 = gca;
ax3.TickLabelInterpreter = 'none';
ax3.TitleHorizontalAlignment = 'left';

ylim([0 max(sctyN3)+1])
xlim([60 100])
xticks([60 70 80 90 100])
xlabel('Percent of consensus epochs')
ylabel('Subjects')

axis square


set(gcf,'Position',[439 184 773 1109])


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


axis square




%% CREATE SORTED FILES
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
subNUM = 6;
[iniBar_s3 , finBar_s3, iniC_s3, finC_s3] = getIndivBarData(subNUM , inLS, inSumS , finSumS , fintSSu);
subTITLE = ['Subject ',num2str(subNUM)];
t1 = tiledlayout(1,3);

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
subtitle([subTITLE,' - Initial Review'])
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
subtitle([subTITLE,' - Final Review'])
axf3 = gca;
axf3.TitleHorizontalAlignment = 'left';
ylabel('Fraction of sleep stage')
yticks([0 0.25 0.5 0.75 1])
axis square

t.TileSpacing = 'compact';
t.Padding = 'compact';

finINchange_s3 = finC_s3 - iniC_s3;
nexttile
generateIndivEpochChangePlot(finINchange_s3,colorMPrgb, subNUM)

%% Group subject plots

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
%

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

axis square

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
axis square

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

        if sleepStC >= 10
            disp(['INDEX found! ', num2str(eii-9)])
            sleepStarInd(iii) = eii - 9;
            break
        end

    end
end


% Loop through all cases and for each - select from 5 mins minus start and
% determine how many steps have changed
allconCkcell = cell(length(inSumS),1);
allconCkcellnoOFF = cell(length(inSumS),1);
for ifi = 1:length(inSumS)

    beginIND = sleepStarInd(ifi) - 10;
    iniTEMP = inSumS{ifi}(beginIND:end);
    fniTEMP = finSumS{ifi}(beginIND:end);

    % Logical
    conSenCk = zeros(1,length(fniTEMP),'logical');
    % What was final stage
    for checkC = 1:length(fniTEMP)

        conSenCk(checkC) = ~matches(iniTEMP{checkC},fniTEMP{checkC});

    end
    
    conSenCkNO = zeros(1,length(finSumS{ifi}),'logical');
    for checkC = 1:length(finSumS{ifi})

        conSenCkNO(checkC) = ~matches(inSumS{ifi}{checkC},finSumS{ifi}{checkC});

    end


    allconCkcell{ifi} = conSenCk;
    allconCkcellnoOFF{ifi} = conSenCkNO;
end


allconCkcell2 = allconCkcell(sleepStarInd < 900,:);
allconCkcellnoOFF2 = allconCkcellnoOFF(sleepStarInd < 900,:);
sleepStarInd2 = sleepStarInd(sleepStarInd< 900);



% Stack up in an array
longestcell = max(cellfun(@(x) numel(x), allconCkcell2));

probCCmat = nan(length(allconCkcell2),longestcell);
for ppi = 1:length(allconCkcell2)

    probCCmat(ppi,1:numel(allconCkcell2{ppi})) = allconCkcell2{ppi};

end

% Stack no offset

longestcellNO = max(cellfun(@(x) numel(x), allconCkcellnoOFF2));

probCCmatNO = nan(length(allconCkcellnoOFF2),longestcellNO);
for ppi = 1:length(allconCkcellnoOFF2)

    probCCmatNO(ppi,1:numel(allconCkcellnoOFF2{ppi})) = allconCkcellnoOFF2{ppi};

end

numChange = sum(probCCmat,'omitnan');
numPresent = sum(~isnan(probCCmat));

thirdNightAvail = find(numPresent < 15, 1, 'first');
binChange = probCCmatNO(:,1:700);
numChange2 = numChange(1:thirdNightAvail);
numPresent2 = numPresent(1:thirdNightAvail);
%%
figure;
binChange2 = binChange;

% for bib = 1:length(sleepStarInd2)
% 
%     binChange2(bib,sleepStarInd2(bib)) = 10;
% 
% end

% imagesc(binChange2);colormap('parula')

for liNe = 1:44

    line([find(binChange2(liNe,:) == 1) ; find(binChange2(liNe,:) == 1)],...
        [repmat(liNe - 0.5,1, sum(binChange2(liNe,:) == 1));...
        repmat(liNe + 0.5,1, sum(binChange2(liNe,:) == 1))],'Color',[0.5 0.5 0.5],...
        'LineWidth',0.5)

    hold on

    line([sleepStarInd2(liNe) ; sleepStarInd2(liNe)],...
        [repmat(liNe - 0.5,1,1);...
        repmat(liNe + 0.5,1,1)],'Color','r','LineWidth',1.5)
end

axis square
yticks(1:2:44)
ylabel('Subject/Night #')
xticks([1 175 350 525 700])
xlabel('Epoch #')
text(3,44,'| Start of sleep onset','Color','r')

%%

figure;
probChange = numChange2./numPresent2;

plot(probChange,'Color',[0.7176 0.2745 1.0000],'LineWidth',2)
xline(10,'-', 'Start of sleep onset')
yline(mean(probChange),'-','Mean fraction of epochs with disagreement')
xlim([1 300])
ylabel('Fraction of nights epoch changed with consensus')
xlabel('Epoch #: aligned for all nights by first contiguous 5 minutes of sleep')

axis square

%% Sanity check don't align and look for dead zones at sleep onset

for ifi = 6

    beginIND = sleepStarInd2(ifi);
    iniTEMP = inSumS{ifi};
    fniTEMP = finSumS{ifi};

    % Logical
    conSenCk = zeros(1,length(fniTEMP),'logical');
    % What was final stage
    for checkC = 1:length(fniTEMP)

        conSenCk(checkC) = ~matches(iniTEMP{checkC},fniTEMP{checkC});

    end

    plot(conSenCk,'Color',[0.7176 0.2745 1.0000],'LineWidth',1)
    xlim([1 800])
    xline(beginIND,'-', 'Start of sleep onset','LineWidth',3,'Color','k')
    title(num2str(ifi))
    pause

end



























