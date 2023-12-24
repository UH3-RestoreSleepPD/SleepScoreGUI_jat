function [] = aim2_rawJSONcheck(nameValueInputs)


arguments
    nameValueInputs.JsonDir = 0;
    nameValueInputs.Aim = 1;
    nameValueInputs.Save = 0;
end


argsCell = namedargs2cell(nameValueInputs);
argsTable = cell2table(argsCell(2:2:length(argsCell)),'VariableNames',...
    argsCell(1:2:length(argsCell)));

if ~iscell(argsTable.JsonDir)
    mainDIR = uigetdir();
else
    mainDIR = argsTable.JsonDir{1};
end

cd(mainDIR)

% Dig down to find JSON files
tmpFname = 'RawDataTD.json';
dataRead = jsondecode(fileread(tmpFname));

test = 1;

dataReadTD = dataRead.TimeDomainData;
sample1 = [];
tmpSecs = zeros(height(dataReadTD),1);
fracUnqiue = zeros(height(dataReadTD),1);
sampleNUM = zeros(height(dataReadTD),1)
for di = 1:height(dataReadTD)

    tmpSecs(di) = dataReadTD(di).Header.timestamp.seconds;

    sample1 = [sample1 ; dataReadTD(di).ChannelSamples(1).Value];

    % plot(dataReadTD(di).ChannelSamples(1).Value)
    % hold on

    sampleNUM(di) = numel(dataReadTD(di).ChannelSamples(1).Value);

    allVALS = dataReadTD(di).ChannelSamples(1).Value;

    fracUnqiue(di) = numel(unique(allVALS))/length(allVALS);

end
plot(sampleNUM)
ylabel('Number of elements in each packet')
xlabel('Packet ID')


numel(unique(sample1))/length(sample1)

% plot(tmpSecs)

x = dataReadTD(1).ChannelSamples(1).Value;
y = dataReadTD(3).ChannelSamples(1).Value;

[c,lags] = xcorr(x,y,10,'normalized');
stem(lags,c)








end