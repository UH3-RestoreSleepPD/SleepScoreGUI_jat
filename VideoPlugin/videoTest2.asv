%% NEED TO ADD TTL START
load('F210219-0001.mat');
numSecsTS = CDIG_IN_1_Down/44000;
numSampsTS = CDIG_IN_1_Down;

if numSecsTS < 30
    % remove first epoch
    
end

%% CD to video folder
cd('D:\01_Coding_Datasets\SLEEP VIDEO\UNMC2_N1')

%% Load video file in
tic
v3 = VideoReader('3_UNMC_1.mp4');
toc
%% Numbers
% Number of frames
tic
numF = v3.NumFrames;
toc
% Frame rate
fr = 7;
% Number of frames per 30 second bin
fepoch = fr*30;
% Total Number of 30 second bins
numSecs = v3.NumFrames/fr;
numBins = ceil(numF/fepoch);
% Cell array at number of bins with 3D matrices with image sets gray
binSTOP = numBins - 1;
maxWoBin = fepoch*binSTOP;

startS = transpose(1:fepoch:maxWoBin);
stopS = [startS(2:end) - 1 ; maxWoBin];
allVinds = [startS , stopS];
allVinds(end+1,1) = allVinds(end,2) + 1;
allVinds(end,2) = numF;
app.VidINDs = allVinds;


%% Test Epoch structure - 30 seconds = 210 frames
vidWidth = v3.Width;
vidHeight = v3.Height;

mov = struct('cdata',zeros(round(vidHeight/2),round(vidWidth/2),1,'uint8'),...
    'colormap','gray');

start = 1;
stop = fepoch;
stepSize = fepoch;

for en = 1:20 % loop through bins
    
    disp(['frame ', num2str(en), ' out of ', num2str(numBins)])
    tic
    tmpFrames = read(v3,[start stop]);
    
    for k = 1:size(tmpFrames,4)
        tmpF = tmpFrames(:,:,:,k);
        bwtmp = im2gray(tmpF);
        reStmp = imresize(bwtmp, 0.5);
        mov(k).cdata = reStmp;
        
        disp([num2str(k) ' out of ' num2str(v3.NumFrames)])
    end
    
    saveLoc = ['D:\01_Coding_Datasets\SLEEP VIDEO\UNMC2_N1\VideoEpochMat\',...
        'epoch_',num2str(en),'.mat'];
    save(saveLoc,'mov');
    toc
    start = stop + 1;
    stop = start + stepSize - 1;
    
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