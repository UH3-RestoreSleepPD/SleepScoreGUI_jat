function [] = LW_processSumTab()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

cd('E:\Dropbox\Publications_Meta\InProgress\LWest_ScoreConsensus2022\Extra material\StatsAnalysis')

% Load Initial data
sumtab = readtable("SummaryTable.csv");

subTrim = split(sumtab.FileN,'_');
subNNu = subTrim(:,1);
subIIu = subTrim(:,2);
subStack = cellfun(@(x,y) [x , '_' , y], subNNu, subIIu,'UniformOutput',false);
uniSub = unique(subStack);
sumtab.SStack = subStack;


%%
newTable = table;

for ui = 1:length(uniSub)

    tmpSub = uniSub{ui};
    subTab = sumtab(matches(sumtab.SStack,tmpSub),:);

    tps = {'I','F'};
    for ti = 1:2

        tmpTP = tps{ti};
        tpTab = subTab(matches(subTab.TimeP,tmpTP),:);

        uniRats = unique(tpTab.RaterID);
        for ri = 1:length(uniRats)

            tmpRat = uniRats{ri};
            ratTab = tpTab(matches(tpTab.RaterID,tmpRat),:);

            uniSSt = unique(ratTab.SStage);
            for si = 1:length(uniSSt)

                tmpSStg = uniSSt{si};
                ssTab = ratTab(matches(ratTab.SStage,tmpSStg),:);

                tmptable = ssTab(1,:);
                tmptable.CountS = sum(ssTab.CountS);

                newTable = [newTable ; tmptable];



            end
        end
    end
end


newTable2 = removevars(newTable,'FileN');

writetable(newTable2,'SumSummryTab.csv'); 


end