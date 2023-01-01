function [] = generateIndivEpochChangePlot(finCHANGE,colorMAP, subNUM)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

scatXax_s3 = reshape(repmat(transpose(1:3),1,6),numel(repmat(transpose(1:3),1,6)),1);
scatYax_s3 = reshape(finCHANGE,numel(finCHANGE),1);
scatcolor_s3 = zeros(18,3);
scatcolor_s3(scatXax_s3 == 1,:) = colorMAP;
scatcolor_s3(scatXax_s3 == 2,:) = colorMAP;
scatcolor_s3(scatXax_s3 == 3,:) = colorMAP;

scatter(scatXax_s3,scatYax_s3,100,scatcolor_s3,'filled')
hold on
linEdata = transpose(finCHANGE);
% Line 
for pi_s3 = 1:6

    line([1 2;2 3], [linEdata(pi_s3,1) linEdata(pi_s3,2); linEdata(pi_s3,2) linEdata(pi_s3,3)],...
        'Color',colorMAP(pi_s3,:))

end

yline(0,'LineWidth',1,'LineStyle','-')
xlim([0.5 3.5])
xticks([1 2 3])
xticklabels({'Night 1' , 'Night 2', 'Night 3'})
ylabel('Epoch # change between Initial and Final review')

subNUMsub = ['Subject ' num2str(subNUM)];
subtitle(subNUMsub)
axf3b = gca;
axf3b.TitleHorizontalAlignment = 'left';
axis square



end