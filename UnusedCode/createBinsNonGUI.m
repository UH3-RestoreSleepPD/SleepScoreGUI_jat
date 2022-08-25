function outTTab = createBinsNonGUI(app, chan,chanLABS,SF)

switch chan
    case 'emg'
        rawDAT = app.emgTTpp;
    case 'eog'
        rawDAT = app.eogTTpp;
    case 'eeg'
        rawDAT = app.eegTTpp;
    case 'lfp'
        rawDAT = app.lfpTTpp;
    case 'lfpr'
        rawDAT = app.lfpTTNdnS;
end

num30bins = floor(height(rawDAT)/SF/30);
tstart = seconds(0);
tend = seconds(30);
outBINS = cell(num30bins,size(rawDAT,2));
for si = 1:num30bins

    timeLimits = timerange(tstart,tend);

    tmpT = rawDAT(timeLimits,:);
    tt2t = timetable2table(tmpT);
    t2c = table2cell(tt2t);
    cleanCel = cell2mat(t2c(:,2:end));

    for cchi = 1:size(rawDAT,2)
        outBINS{si,cchi} = cleanCel(:,cchi);
    end

    tstart = tend;
    tend = tstart + seconds(30);

    disp([num2str(si) ,' out of ' , num2str(num30bins)])

end

% Create new TimeTable
% Time
timeBINS = transpose(seconds(0:30:30*num30bins-1));
% Data
finalTT = cell2table(outBINS);
% Variable names
varNAMES = chanLABS;
% Create Table
finalTT.Properties.VariableNames = varNAMES;
finalTT.Time = timeBINS;

outTTab = table2timetable(finalTT);


end