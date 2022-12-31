% cd('D:\Dropbox\Publications_Meta\InProgress\LWest_ScoreConsensus2022\Extra material')
cd('E:\Dropbox\Publications_Meta\InProgress\LWest_ScoreConsensus2022\Extra material')
%%

% dataTab = readtable('DataForAnalysis.xlsx');
dataTab = readtable('DataForAnalysis2.xlsx');

%% 

uniqueNIGHTs = unique(dataTab.CaseID);

nightPer = zeros(length(uniqueNIGHTs),2);
for ui = 1:length(uniqueNIGHTs)

    for ic = 1:2

        switch ic
            case 1
                tmpTAB = dataTab(matches(dataTab.StageID,'Initial'),:);
            case 2
                tmpTAB = dataTab(matches(dataTab.StageID,'ConsensusF'),:);
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

% dumbbell plot
sNperTab = sortrows(nPerTab,2,'ascend');
scYv = transpose(height(sNperTab):-1:1);
scatter(sNperTab.Initial,scYv,40,"red",'filled')
hold on
scatter(sNperTab.Final,scYv,40,"blue",'filled')
xVals = [transpose(sNperTab.Initial) ; transpose(sNperTab.Final)];
yVals = [transpose(scYv) ; transpose(scYv)];
line(xVals,yVals,'Color','k')
yticks(1:height(sNperTab))
yticklabels(sNperTab.Ylab);

ax = gca;
ax.TickLabelInterpreter = 'none';

ylim([1 22])
xlim([70 100])
xticks([70 80 90 100])
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
%%

dataTabGGa = dataTab;
% dataTabGGa = dataTabGGa(:,[1,2,4,5]);

% Break up subject and night
for sni2 = 1:height(dataTabGGa)
    tmpRow = dataTabGGa.CaseID{sni2};
    tmpParts = split(tmpRow,{'_','.'});
    patNUM = [tmpParts{1},tmpParts{2}];
    nightNUM = tmpParts{3};
    dataTabGGa.SubA{sni2} = patNUM;
    dataTabGGa.Sub{sni2} = tmpParts{1};
    dataTabGGa.night{sni2} = nightNUM;
end

% uniSUBS = unique(dataTabGGa.Sub);
% for ui = 1:length(unique(dataTabGGa.Sub))
% 
%     tmpUNI = uniSUBS{ui};
%     dataTabGGa.Sub(matches(dataTabGGa.Sub,tmpUNI)) = {num2str(ui)};
% 
% end

intiLOCs = matches(dataTabGGa.StageID,'Initial');
dataTabGGa.SubIC = num2cell(num2str(zeros(height(dataTabGGa),1)));
dataTabGGa.SubIC(intiLOCs) = cellfun(@(x) num2str(x), transpose(num2cell(1:sum(intiLOCs))),...
    "UniformOutput",false);
conFLocs = matches(dataTabGGa.StageID,'ConsensusF');
dataTabGGa.SubIC(conFLocs) = cellfun(@(x) num2str(x), transpose(num2cell(1:sum(conFLocs))),...
    "UniformOutput",false);
% dataTabGGa = dataTabGGa(:,[1,2,4,5,6,7]);

dataTabGGa.StageID(matches(dataTabGGa.StageID,"ConsensusF")) = {'Final'};


%%
allTable = table;
for ci = 1:3

    switch ci
        case 1 % Consensus

            tmpTTable = dataTabGGa(:,[2,5]); % ADD CASEID
%             tmpTTable.StageID = replace(tmpTTable.StageID,'ConsensusF','Final');
            tmpTTable.x = repmat({'Consensus'},height(tmpTTable),1);
            tmpTTable = renamevars(tmpTTable,["Count","StageID"],["Freq","Stratum"]);
            tmpTTable.Cohort = transpose(1:height(tmpTTable));

        case 2 % Night

            tmpTTable = dataTabGGa(:,[2,8]);
            tmpTTable.x = repmat({'Night'},height(tmpTTable),1);
            tmpTTable = renamevars(tmpTTable,["Count","night"],["Freq","Stratum"]);
            tmpTTable.Cohort = transpose(1:height(tmpTTable));

        case 3 % Stage

            tmpTTable = dataTabGGa(:,[2,1]);
            tmpTTable.x = repmat({'Stage'},height(tmpTTable),1);
            tmpTTable = renamevars(tmpTTable,["Count","Stage"],["Freq","Stratum"]);
            tmpTTable.Cohort = transpose(1:height(tmpTTable));


    end
allTable = [allTable ; tmpTTable];



end

%%


% Save out for sankey plot in R
writetable(nPerTab,'subjectSankey.csv')
writetable(allTable,'subjectAlluvial.csv')
writetable(dataTabGGa,'subjectAlluvial2.csv')

% Subtract Final from Initial
% Create table - remove U
dataTabNU = dataTab(~matches(dataTab.Stage,'U'),:);
dataTabNU_F = dataTabNU(matches(dataTabNU.StageID,'ConsensusF'),:);
dataTabNU_I = dataTabNU(matches(dataTabNU.StageID,'Initial'),:);
dataTabnuFdI = dataTabNU_I(:,{'Stage','CaseID'});

stageAll = cell(height(dataTabnuFdI),1);
caseAll = cell(height(dataTabnuFdI),1);
fper = zeros(height(dataTabnuFdI),1);
iper = zeros(height(dataTabnuFdI),1);
dper = zeros(height(dataTabnuFdI),1);
for sI = 1:height(dataTabnuFdI)

    tmpCaseF = dataTabNU_F.CaseID{sI};
    caseTABF = dataTabNU_F(matches(dataTabNU_F.CaseID,tmpCaseF),:);
    caseTABI = dataTabNU_I(matches(dataTabNU_I.CaseID,tmpCaseF),:);

    tmpStage = dataTabNU_F.Stage{sI};
    stageTF = caseTABF.Percent(matches(caseTABF.Stage,tmpStage),:);
    stageTI = caseTABI.Percent(matches(caseTABI.Stage,tmpStage),:);

    stageAll{sI} = tmpStage;
    caseAll{sI} = tmpCaseF;
    fper(sI) = stageTF;
    iper(sI) = stageTI;
    dper(sI) = stageTF - stageTI;

end

% dataTabnuFdI.F = dataTabNU_F.Percent;
% dataTabnuFdI.I = dataTabNU_I.Percent;
% dataTabnuFdI.Diff = dataTabNU_F.Percent - dataTabNU_I.Percent;

perTableStage = table(caseAll,stageAll,fper,iper,dper,'VariableNames',...
    {'CaseID','StageID','Final','Initial','Differ'});

gpSum = groupsummary(perTableStage,'StageID','mean','Differ');


save('finalSummaryCon.mat','perTableStage','gpSum')



% 
