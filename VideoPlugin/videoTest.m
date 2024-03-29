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
cd('D:\01_Coding_Datasets\SLEEP VIDEO\UNMC2_N1')


%% Load video file in
v3 = VideoReader('3_UNMC_1.mp4');


%%

vidWidth = v3.Width;
vidHeight = v3.Height;

mov = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),...
    'colormap','gray');

%%

k = 1;
while hasFrame(v3)
    mov(k).cdata = readFrame(v3);
    k = k+1;
    disp([num2str(k) ' out of ' num2str(v3.NumFrames)])
end

%%
currAxes = axes;
for i = 1
    vidFrame = readFrame(v3 , 'Grayscale');
    image(vidFrame, 'Parent', currAxes);
    currAxes.Visible = 'off';
    pause(1/7);
    %     text(1,100,['This is frame ', num2str(i)],'FontSize',40)
end





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
Jtest = imresize(bwVF, 0.5);
image(Jtest);
colormap('gray')
currAxes.Visible = 'off';

%% Test Epoch structure - 30 seconds = 210 frames


mov = struct('cdata',zeros(round(vidHeight/2),round(vidWidth/2),1,'uint8'),...
    'colormap','gray');

kki = 1210:1420;

for k = 1:210
    tmpF = read(v3,kki(k));
    bwtmp = im2gray(tmpF);
    reStmp = imresize(bwtmp, 0.5);
    mov(k).cdata = reStmp;
    
    disp([num2str(k) ' out of ' num2str(v3.NumFrames)])
end


%% play through

currAxes = axes;
for i = 1:length(mov)
%     image(mov(i).cdata);
    image(mov(i).cdata, 'Parent', currAxes);
    currAxes.Visible = 'off';
    colormap('gray')
    pause(1/7)
end

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