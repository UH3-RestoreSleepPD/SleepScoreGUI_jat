% main location
% C:\Users\Admin\Downloads\AccelerometerTesting\AccelerometerTesting


trialNames = Acctable.Properties.VariableNames;

trialNums = extractAfter(trialNames,'0');
startStops = extractBefore(trialNames,'0');

trialTypesN = sum(cellfun(@(x) ~isempty(x), table2cell(Acctable(:,1)),...
    'UniformOutput',true));

allTime = GeneratedData.localTime;
Timezonechange = allTime.TimeZone;

startTime1 = Acctable.Start01{1,1};
startTime1.TimeZone = Timezonechange;

[~ , x] = min(abs(allTime - startTime1));




%%



xSample = Acctable.epsilon{2,1}{:,8};
xSample = xSample(~isnan(xSample));
ySample = Acctable.epsilon{2,1}{:,9};
ySample = ySample(~isnan(ySample));
zSample = Acctable.epsilon{2,1}{:,10};
zSample = zSample(~isnan(zSample));

timeSample = size( Acctable.epsilon{2,1}{:,8},1);
timeSample2 = (timeSample/250)/60;


rcTime = Acctable.epsilon{2,1}{:,1};
rcTime2 =rcTime(~isnan(Acctable.epsilon{2,1}{:,8}));
rcTimezone = Acctable.epsilon{2,1}{1,1};
Timezonechange = rcTimezone.TimeZone;
Tone = Acctable.alpha{1,1};
Ttwo = Acctable.beta{1,1};
Tthree = Acctable.gamma{1,1};
Tfour = Acctable.delta{1,1};
Tfive = Acctable.epsilon{1,1};
Tone.TimeZone = Timezonechange;
Ttwo.TimeZone = Timezonechange;
Tthree.TimeZone = Timezonechange;
Tfour.TimeZone = Timezonechange;
Tfive.TimeZone = Timezonechange;

time1find = find(rcTime2 > Tone,1,'first');
time2find = find(rcTime2 > Ttwo,1,'first');
time3find = find(rcTime2 > Tthree,1,'first');
time4find = find(rcTime2 > Tfour,1,'first');
time5find = find(rcTime2 > Tfive,1,'first');

figure()
subplot(3,1,1)
plot(xSample)

xline(time1find)
xline(time2find)
xline(time3find)
xline(time4find)
xline(time5find)

subplot(3,1,2)
plot(ySample)

xline(time1find)
xline(time2find)
xline(time3find)
xline(time4find)
xline(time5find)

subplot(3,1,3)
plot(zSample)
xline(time1find)
xline(time2find)
xline(time3find)
xline(time4find)
xline(time5find)








