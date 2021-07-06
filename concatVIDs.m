function [outputVideo] = concatVIDs(vidLOC)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

cd(vidLOC)
vidDir1 = dir('*.mp4');
vidDir2 = {vidDir1.name};

% Load first video to extract FrameRate

% vid1 = VideoReader(vidDir2{1});
% new video
outputVideo = VideoWriter('newvideo.mp4');
outputVideo.FrameRate = 7;
open(outputVideo);

for vi = 1:length(vidDir2)
    
    tmpVid = VideoReader(vidDir2{vi});
    numFrames = tmpVid.NumFrames;
    
    tframe = 1;
    while hasFrame(tmpVid)
        img1 = readFrame(tmpVid);    % play video
        imgt = img1;
        
        % record new video
        writeVideo(outputVideo, imgt);
        disp(['Frame ', num2str(tframe), ' of out ' , num2str(numFrames)])
        tframe = tframe + 1;
    end
    disp(['Video file ',num2str(vi),' of ' ,num2str(length(vidDir2)), ' done'])
    
    
    
end


% while hasFrame(vid1)
%     img1 = readFrame(vid1);    % play video
%     imgt = img1;
%     
%     % record new video
%     writeVideo(outputVideo, imgt);
% end
% [rows1, columns1, numColors1] = size(img1);
% vid2 = VideoReader('second.avi');
% while hasFrame(vid2)
%     img2 = readFrame(vid2);    % play video
%     imgt = imresize(img2, rows1, columns1);
%     
%     % record new video
%     writeVideo(outputVideo, imgt);
% end
close(outputVideo);




end

