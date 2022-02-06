v = VideoWriter('peaks1.avi');
open(v);

for k = 1:length(mov)
    h = axes;
    set(h,'position',[0 0 1 1])
    image(mov(k).cdata)
    yticks([])
    xticks([])
    colormap('gray')

    frame = getframe(gcf);
    writeVideo(v,frame);
end

close(v);
%%

% Initialize image
app.Figure = figure('Name','Main Video GUI');
app.VidAx = axes(app1.Figure,'XColor','none','YColor','none',...
    'XLim',[0 1],'YLim',[0 1]);
% Can't remember, but here if your video is upside down:
% app1.VidAx.YDir = 'reverse';
app.VFR= vision.VideoFileReader('peaks.avi');
[h,w] = app.VFR.VideoSize;
app.VideoImage = imagesc(app.VidAx,linspace(0,1,w),linspace(0,1,h),zeros(h,w));
...
    % Rest of app constructor
function ButtonPushed(app, event)
% Reset the video first
reset(app.VFR);
% Play video. Every call to the step method reads another frame.
while ~isDone(app.VFR)
    % Update image frame
    app.VideoImage.CData = step(app.VFR);
    drawnow limitrate; % Limits to 20 fps

end
end