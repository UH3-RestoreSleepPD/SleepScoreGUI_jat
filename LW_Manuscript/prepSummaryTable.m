function [] = prepSummaryTable()
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

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
load('summaryTable.mat','dataTable')

% Create Instit / Sub / Date columns

tmpRem = replace(dataTable.FileN,'.mat','');
tmpSplit = split(tmpRem,'_');

dataTable.SubNum = tmpSplit(:,1);
dataTable.Institut = tmpSplit(:,2);
dataTable.NightNum = tmpSplit(:,3);

% Compute percent by Scorer and Date

% Unique File ids
fnIds = unique(dataTable.FileN);

for fi = 1:height(fnIds)

    % Use FileN
    % Get Unique file N
    filnTable = dataTable(matches(dataTable.FileN,fnIds{fi}),:);
    filnInd = matches(dataTable.FileN,fnIds{fi});

    % Get Unique TimePoint
    tpS = {'I','F'};
    for tpi = 1:2
        tpTable = filnTable(matches(filnTable.TimeP,tpS{tpi}),:);
        tpInd = matches(dataTable.TimeP,tpS{tpi});

        scOrErs = unique(tpTable.RaterID);
        for scI = 1:length(scOrErs)
            scTable = tpTable(matches(tpTable.RaterID,scOrErs{scI}),:);
            scInd = matches(dataTable.RaterID,scOrErs{scI});

            totCount = sum(scTable.CountS);
            perCent = round((scTable.CountS / totCount)*100,2);

            finalInd = filnInd & tpInd & scInd;

            dataTable.PercentS(finalInd) = perCent;

        end
    end
end

% dataTable2 = renamevars(dataTable,["TimeP","FileN","CountS","PercentS","SStage"],...
%     ["StageID","CaseID","Count","Percent","Stage"]);




cd(maindir.save)
writetable(dataTable,'DataForAnalysis2.csv')
disp('DONE!')


end