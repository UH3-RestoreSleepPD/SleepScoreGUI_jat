function [] = epoch2avi(epochFold,aviFold)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

cd(epochFold)

epochLIST1 = dir('*.mat');
epochLIST2 = {epochLIST1.name};

for ei = 1:length(epochLIST2)

    cd(epochFold)
    epName = epochLIST2{ei};
    epNparts = strsplit(epName,'.');

    load(epName,'mov')

    cd(aviFold)
    aviName = [epNparts{1},'.avi'];
    tempV = VideoWriter(aviName);
    open(tempV);

    for k = 1:length(mov)
        h = axes;
        set(h,'position',[0 0 1 1])
        image(mov(k).cdata)
        yticks([])
        xticks([])
        colormap('gray')

        frame = getframe(gcf);
        writeVideo(tempV,frame);
    end

    close(tempV);
    close all
    disp(['Video of ',num2str(ei), ' done!'])

end








end