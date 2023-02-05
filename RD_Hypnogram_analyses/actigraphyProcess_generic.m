function [] = actigraphyProcess_generic(subID,dirEL)


switch dirEL
    case 1
        mainLOC = 'D:\Dropbox\Publications_Meta\InProgress\ABaumgartner_Percept2020\Data\';
        % subID = 2;
        subNUM = ['SPPD' , num2str(subID)];
        actLOC = [mainLOC,subNUM,'\ACT_data'];

        saveDIR = [actLOC , '\Summary'];
        if ~exist(saveDIR,'dir')
            mkdir(saveDIR)
        end
        cd(actLOC)
    case 2
        mainLOC = 'C:\Users\Admin\Dropbox\Publications_Meta\InProgress\ABaumgartner_Percept2020\Data\';
        % subID = 2;
        subNUM = ['SPPD' , num2str(subID)];
        actLOC = [mainLOC,subNUM,'\ACT_data'];

        saveDIR = [actLOC , '\Summary'];
        if ~exist(saveDIR,'dir')
            mkdir(saveDIR)
        end
        cd(actLOC)
    case 3
        mainLOC = uigetdir;
        cd(mainLOC)
        subNUM = subID;
        saveDIR = 'D:\VRAWLFP\ActigraphyProcessALL\PRE_ALL_Mat';
end

allCELL = csvFIND();

%--------STATS
statsT = find(contains(allCELL , 'Interval Type'));
statsB1 = contains(allCELL , 'Column Title');
statsB2 = find(statsB1, 1, 'first') - 2;

statsCell1 = allCELL(statsT:statsB2);

statsCell2 = statsCell1(~contains(statsCell1,'EXCLUDED'));
statsCellt = cell(length(statsCell2),38);
for si = 1:length(statsCell2)
    tmpSC = statsCell2{si};
    tmpSCsp = strsplit(tmpSC , ',');
    statsCellt(si,1:length(tmpSCsp)) = tmpSCsp;
end

% Fix Titles
titleRow1 = statsCellt(1,:);
titleRow2 = statsCellt(2,:);

if any(cellfun(@(x) isfloat(x), titleRow1, 'UniformOutput',true))
    numbCells1 = repmat({''}, 1, sum(cellfun(@(x) isfloat(x), titleRow1, 'UniformOutput',true)));
    titleRow1(cellfun(@(x) isfloat(x), titleRow1, 'UniformOutput',true))  = numbCells1;
end

if any(cellfun(@(x) isfloat(x), titleRow2, 'UniformOutput',true))
    numbCells2 = repmat({''}, 1, sum(cellfun(@(x) isfloat(x), titleRow2, 'UniformOutput',true)));
    titleRow2(cellfun(@(x) isfloat(x), titleRow2, 'UniformOutput',true))  = numbCells2;
end

titleRow3 = cellfun(@(x) [strrep(x,'"','') , ' minutes'] ,...
    titleRow1(contains(titleRow2,'minutes')), 'UniformOutput' , false);

statsCellt(1,contains(titleRow2,'minutes')) = titleRow3;
% Remove the 2nd row and last column
statsCellt2 = statsCellt([1,3:size(statsCellt,1)],1:size(statsCellt,2)-1);
% Remove quotations
statsCellt3 = cellfun(@(x) strrep(x,'"',''),statsCellt2,...
    'UniformOutput' , false);
% Create a Table
dataFt_S = statsCellt3(2:end,:);

% clean empty columns
emptyCols = cellfun(@(x) isempty(x), dataFt_S(1,:), 'UniformOutput',true);
dataFt_S2 = dataFt_S(:,~emptyCols);

col_S = statsCellt3(1,~emptyCols);
statsTable = cell2table(dataFt_S2 , 'VariableNames',col_S);
% Save Stats Table
saveFileMAT_S = [saveDIR , filesep , subNUM , '_ACT_STATS.mat'];
saveFileCSV_S = [saveDIR , filesep , subNUM , '_ACT_STATS.csv'];

save(saveFileMAT_S , 'statsTable');
writetable(statsTable , saveFileCSV_S);


%--------EVENTS
eventsT1 = contains(allCELL , 'Column Title');
eventsT2 = find(eventsT1, 1, 'first') + 7;
eventsB1 = contains(allCELL , 'Column Title');
eventsB2 = find(eventsB1, 1, 'last') - 2;

eventsCell1 = allCELL(eventsT2:eventsB2);

eventsCellt = cell(length(eventsCell1),6);
for ei = 1:length(eventsCell1)
    tmpEC = eventsCell1{ei};
    tmpSCev = strsplit(tmpEC , ',');
    eventsCellt(ei,:) = tmpSCev;
end

% Remove last column
eventsCellt1 = eventsCellt(:,1:size(eventsCellt,2)-1);

% Remove quotations
eventsCellt2 = cellfun(@(x) strrep(x,'"',''),eventsCellt1,...
    'UniformOutput' , false);

% Create a Table
dataFt_E = eventsCellt2(2:end,:);
col_E = eventsCellt2(1,:);
eventsTable = cell2table(dataFt_E , 'VariableNames' , col_E);
% Save Stats Table
saveFileMAT_E = [saveDIR , filesep , subNUM , '_ACT_EVENTS.mat'];
saveFileCSV_E = [saveDIR , filesep , subNUM , '_ACT_EVENTS.csv'];

save(saveFileMAT_E , 'eventsTable');
writetable(eventsTable,saveFileCSV_E);

%--------DATA
dataT1 = contains(allCELL , 'Column Title');
dataT2 = find(dataT1, 1, 'last') + 14;

dataCell1 = allCELL(dataT2:length(allCELL));

dataCellt = cell(length(dataCell1),13);
for di = 1:length(dataCell1)
    tmpDC = dataCell1{di};
    tmpDCev = strsplit(tmpDC , ',');
    dataCellt(di,:) = tmpDCev;
end

% Remove last column
dataCellt1 = dataCellt(:,1:size(dataCellt,2)-1);

% Remove quotations
dataCellt2 = cellfun(@(x) strrep(x,'"',''),dataCellt1,...
    'UniformOutput' , false);

% Create a Table
dataFt_D = dataCellt2(2:end,:);
col_D = dataCellt2(1,:);
dataTable = cell2table(dataFt_D , 'VariableNames' , col_D);
% Save Stats Table
saveFileMAT_D = [saveDIR , filesep , subNUM , '_ACT_DATA.mat'];
saveFileCSV_D = [saveDIR , filesep , subNUM , '_ACT_DATA.csv'];

save(saveFileMAT_D , 'dataTable');
writetable(dataTable,saveFileCSV_D);

end





function [outCELL] = csvFIND()

csvDir = dir('*.csv');
csvNames = {csvDir.name};

fid = fopen(csvNames{1});
tline = fgetl(fid);
lineCount = 1;
outCELL = cell(1,1000000);
while ischar(tline)
    outCELL{lineCount} = tline;
    lineCount = lineCount + 1;
    tline = fgetl(fid);
end
fclose(fid);

outCELL = outCELL(cellfun(@(x) ~isempty(x) , outCELL , 'UniformOutput',true));

end

