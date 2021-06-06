%% CD to video folder
cd('D:\01 Coding Datasets\SLEEP VIDEO\UNMC2_N1')


%% Load video file in

v3 = VideoReader('3_UNMC_1.mp4');

%%

currAxes = axes;
while hasFrame(v3)
    vidFrame = readFrame(v3);
    image(vidFrame, 'Parent', currAxes);
    currAxes.Visible = 'off';
    pause(1/7);
end

%%

%% CD to video folder
cd('D:\01 Coding Datasets\SLEEP VIDEO\UNMC2_N1')


%% Load video file in
v3 = VideoReader('3_UNMC_1.mp4');

%% Numbers

% Number of frames
numF = v3.NumFrames;
% Frame rate
fr = 7;
% Number of frames per 30 second bin
fepoch = fr*30;
% Total Number of 30 second bins
numSecs = v3.NumFrames/fr;
numBins = ceil(numF/fepoch);
% Cell array at number of bins with 3D matrices with image sets gray

%% Load temp frame
vidFrame = read(v3,f);

%% Edit frame
close 
currAxes = axes;
bwVF = im2gray(vidFrame);
image(bwVF);
colormap('gray')
currAxes.Visible = 'off';
%% Create save structure

% app.videoFile;
% videoFile = VideoReader('3_UNMC_1.mp4');
% app.videoFile = videoFile;
% app.epochIndex [start and stop]

% Number of frames
numF = v3.NumFrames;
% Frame rate
fr = 7;

fepoch = fr*30;

numBins = ceil(numF/fepoch);
binSTOP = numBins - 1;
maxWoBin = fepoch*binSTOP;

startS = transpose(1:fepoch:maxWoBin);
stopS = [startS(2:end) - 1 ; maxWoBin];
allVinds = [startS , stopS];
allVinds(end+1,1) = allVinds(end,2) + 1;
allVinds(end,2) = numF;
app.VidINDs = allVinds;

% function for processing
% function [vidEpoch] = grabVidEpoch(app,epoch)
% 
% vidEpoch = zeros(app.Fwidth,app.Fheight,app.fepoch,'uint8');
% 
% startF = app.VidINDs(epoch,1);
% stopF = app.VidINDs(epoch,2);
% 
% tmpFrames = read(v3,[startF stopF]);
% 
% for ti = 1:app.fepoch
%     
%     ttFrm = im2gray(tmpFrames(:,:,:,ti));
%     vidEpoch(:,:,ti) = ttFrm;
%     
% end
% 
% end


%%



mainBINS = cell(numBins,1);

start = 1;
stop = fepoch;
stepSize = fepoch;
for en = 1:numBins % loop through bins
   
    disp(['frame ', num2str(en), ' out of ', num2str(numBins)])
    tmpBin = zeros(size(vidFrame,1),size(vidFrame,2),fepoch,'uint8');
    
    tmpFrames = read(v3,[start stop]);
    
    for ti = 1:size(tmpFrames,4)
       
        ttFrm = im2gray(tmpFrames(:,:,:,ti));
        tmpBin(:,:,ti) = ttFrm;
        
    end
    
    mainBINS{en} = tmpBin;
    
    start = stop + 1;
    stop = start + stepSize - 1;
    
     
end