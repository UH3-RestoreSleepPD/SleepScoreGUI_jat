function [] = actigraphyProcess_RD_pre(dirEL)

mainLOC = dirEL;
cd(mainLOC)

saveDIR = 'D:\VRAWLFP\ActigraphyProcessALL\PRE_ALL_Mat';

csvALL = dir('*.csv');
csvALL2 = {csvALL.name};

for cci = 13:length(csvALL2)

    csvName = csvALL2{cci};

    nameParts = split(csvName,'_');
    saveNAMEi = [nameParts{1}, '_' nameParts{2}];

    allCELL = csvFIND(csvName);

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
    saveFileMAT_S = [saveDIR , filesep , saveNAMEi , '_ACT_STATS.mat'];
    saveFileCSV_S = [saveDIR , filesep , saveNAMEi , '_ACT_STATS.csv'];

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
    saveFileMAT_E = [saveDIR , filesep , saveNAMEi , '_ACT_EVENTS.mat'];
    saveFileCSV_E = [saveDIR , filesep , saveNAMEi , '_ACT_EVENTS.csv'];

    save(saveFileMAT_E , 'eventsTable');
    writetable(eventsTable,saveFileCSV_E);

    %--------DATA
    dataT1 = contains(allCELL , 'Column Title');

    if matches(saveNAMEi,{'UNMC_11','UNMC_12','UNMC_13','UNMC_17','UNMC_1',...
            'UNMC_20','UNMC_21','UNMC_22','UNMC_2','UNMC_3','UNMC_6',...
            'UNMC_8','UNMC_9','UPEN_1','UPEN_2','UPEN_3'})
        dataT2 = find(dataT1, 1, 'last') + 15;
    else
        dataT2 = find(dataT1, 1, 'last') + 14;
    end

    dataCell1 = allCELL(dataT2:length(allCELL));

    dataCellt = cell(length(dataCell1),13);
    for di = 1:length(dataCell1)
        tmpDC = dataCell1{di};
        tmpDCev = strsplit(tmpDC , ',');

%         if any(cellfun(@(x) isempty(x), tmpDCev, 'UniformOutput',true))
%             tmpDCev = tmpDCev(~cellfun(@(x) isempty(x), tmpDCev, 'UniformOutput',true));
%         end
        if length(tmpDCev) ~= 13
            tmpDCev = tmpDCev(1:13);
        end

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
    saveFileMAT_D = [saveDIR , filesep , saveNAMEi , '_ACT_DATA.mat'];
    saveFileCSV_D = [saveDIR , filesep , saveNAMEi , '_ACT_DATA.csv'];

    save(saveFileMAT_D , 'dataTable');
    writetable(dataTable,saveFileCSV_D);

end

end





function [outCELL] = csvFIND(csvName)

fid = fopen(csvName);
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

