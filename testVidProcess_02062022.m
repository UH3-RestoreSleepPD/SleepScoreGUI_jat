currAxes = axes;
for i = 1:300
    vidFrame = read(v,i);
    image(vidFrame, 'Parent', currAxes);
    currAxes.Visible = 'off';
    pause(1/(v.FrameRate*5));
end


%% CREATE BINS
nFrames = v.NumFrames;
nBins = round(nFrames/300);
totFrames = 300*nBins;
bStart = transpose(1:300:totFrames);
bStop = bStart(2:end)-1;
bStart = bStart(1:end-1);
bBoth = [bStart,bStop];