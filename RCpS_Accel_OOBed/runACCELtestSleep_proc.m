function [] = runACCELtestSleep_proc(plottest)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


load("experBlockTrialRaw.mat",'blockSc');
epochDATA = blockSc;

load("jattimes3.mat",'Acctable');
allDATA = Acctable;

[xSampleE , ySampleE , zSampleE] = cleanSamples(allDATA);


if plottest
    runPLOT(epochDATA)
end

%%% FIRST ATTEMPT - Realign based on first peak
% test block 1
blockLENGTHS = [160, 110, 160, 110];
blockStrT = [49, 39, 49, 39];
blockEndT = [110, 70, 110, 70];

allblocks = cell(1,4);
for bxb = 1:4
    blocktemp = epochDATA{bxb};
    % initial peak
    block1x = zeros(3,blockLENGTHS(bxb));
    for x3x = 1:3
        xtemp = blocktemp{x3x}(1,:);
        xtm = abs(abs(xtemp) - mean(abs(xtemp)));
        % plot(xtm)

        [~,peakLoc] = findpeaks(xtm,"NPeaks",1,"MinPeakHeight",30);
        trimStart = peakLoc - blockStrT(bxb);
        trimEnd = peakLoc + blockEndT(bxb);
        block1x(x3x,:) = xtm(trimStart:trimEnd);
    end

    % close all
    figure;
    plot(transpose(block1x))
    allblocks{bxb} = block1x;


end

templateWaveform = mean(allblocks{1});
crossCorrelation = xcorr(xSampleE, templateWaveform);

[peaks, templateOccurrences] = findpeaks(crossCorrelation);


figure;
plot(xSampleE);
hold on;
plot(templateOccurrences, repmat(150,size(templateOccurrences)), 'r', 'LineWidth', 2);

end






function [] = runPLOT(blockSc)


for bi = 1:length(blockSc)
    blocki = blockSc{bi};
    figure;

    for ti = 1:3
        subplot(3,1,1)
        hold on
        xtemp = blocki{ti}(1,:);
        plot(xtemp)

        subplot(3,1,2)
        hold on
        ytemp = blocki{ti}(2,:);
        plot(ytemp)

        subplot(3,1,3)
        hold on
        ztemp = blocki{ti}(3,:);
        plot(ztemp)

    end
end



end





function [xOut, yOut , zOut] = cleanSamples(Acctable)

xOut = Acctable.epsilon{2,1}{:,8};
xOut = xOut(~isnan(xOut));
yOut = Acctable.epsilon{2,1}{:,9};
yOut = yOut(~isnan(yOut));
zOut = Acctable.epsilon{2,1}{:,10};
zOut = zOut(~isnan(zOut));

end