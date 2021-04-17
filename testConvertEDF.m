%% Unpack data

fs = 100;
dataELs = height(TT)*length(TT.A1{1});

rawDat = zeros(dataELs , width(TT) - 1);
start = 1;
stop = 3000;
for t = 1:height(TT)
    
    for t2 = 1:width(TT) - 1
        rawDat(start:stop,t2) = TT{t,t2}{1};
    end
    
    start = stop + 1; 
    stop = start + 2999;
    
end



%% Annotations
start_index = transpose(round(1:3000:dataELs));
Onset = seconds(start_index./fs);
Duration = seconds(ones(length(Onset),1)*30);
Annotations = string(transpose(1:length(Duration)));

%%
annotationslist = timetable(Onset,Annotations,Duration);

%%

hdr = edfheader("EDF+");

hdr.NumDataRecords = 1;
hdr.DataRecordDuration = seconds(dataELs/fs);
hdr.NumSignals = width(TT) - 1;
hdr.SignalLabels = string(TT.Properties.VariableNames(1:end-1));

hdr.PhysicalDimensions = repelem("mV",hdr.NumSignals);

max2u = zeros(1,16);
min2u = zeros(1,16);
for i = 1:16
   
    min2u(i) = min(cellfun(@(x) min(x), TT.(hdr.SignalLabels{i}), 'UniformOutput',true));
    max2u(i) = max(cellfun(@(x) max(x), TT.(hdr.SignalLabels{i}), 'UniformOutput',true));
    
end

hdr.PhysicalMin = round(min2u);
hdr.PhysicalMax = round(max2u);
hdr.DigitalMin = repmat(-32768,1,16);
hdr.DigitalMax = repmat(32767,1,16);

%%
edfw = edfwrite("2_UNMC_1.edf",hdr,rawDat,annotationslist);



