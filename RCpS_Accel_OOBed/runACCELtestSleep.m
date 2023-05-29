function [] = runACCELtestSleep()

% main location
% D:\Dropbox\SleepAccel_Platt
% C:\Users\Admin\Downloads\AccelerometerTesting\AccelerometerTesting

% baseline data
load("baseLINE.mat",'Acctable');
baseDATA = Acctable;

load("jattimes3.mat",'Acctable');
epochDATA = Acctable;

[xSampleB , ySampleB , zSampleB] = cleanSamples(baseDATA);
[xSampleE , ySampleE , zSampleE] = cleanSamples(epochDATA);

% timeSample = size( Acctable.epsilon{2,1}{:,8},1);
% timeSample2 = (timeSample/250)/60;

[timesOUT] = getMarkTimes(Acctable);


figure()
subplot(3,1,1)
plot(xSampleE)
xline(timesOUT(:,1))

subplot(3,1,2)
plot(ySampleE)
xline(timesOUT(:,1))

subplot(3,1,3)
plot(zSampleE)
xline(timesOUT(:,1))


end








function [timesOUT] = getMarkTimes(Acctable)

rcTime = Acctable.epsilon{2,1}{:,1};
rcTime2 =rcTime(~isnan(Acctable.epsilon{2,1}{:,8}));
rcTimezone = Acctable.epsilon{2,1}{1,1};
Timezonechange = rcTimezone.TimeZone;

trialONE = Acctable.alpha;
trialTWO = Acctable.beta;
trialTHREE = Acctable.gamma;
trialFOUR = Acctable.delta;
endOFtrials = Acctable.epsilon{1,1};

alltimeInds = zeros(1,13);

startTIMES = [trialONE{1} , trialONE{2} , trialONE{3} , trialTWO{1}, ...
    trialTWO{2}, trialTWO{3}, trialTHREE{1}, trialTHREE{2}, trialTHREE{3},...
    trialFOUR{1}, trialFOUR{2}, trialFOUR{3}, endOFtrials];

for si = 1:length(startTIMES)

    ttEMPT = startTIMES(si);
    ttEMPT.TimeZone = Timezonechange;

    tsI = find(rcTime2 > ttEMPT,1,'first');
    alltimeInds(si) = tsI;

end

alltimeSta = alltimeInds(1:12);
allTimeEnd = alltimeInds(2:end);

timesOUT = transpose([alltimeSta ; allTimeEnd]);

end



function [xOut, yOut , zOut] = cleanSamples(Acctable)

xOut = Acctable.epsilon{2,1}{:,8};
xOut = xOut(~isnan(xOut));
yOut = Acctable.epsilon{2,1}{:,9};
yOut = yOut(~isnan(yOut));
zOut = Acctable.epsilon{2,1}{:,10};
zOut = zOut(~isnan(zOut));


end


