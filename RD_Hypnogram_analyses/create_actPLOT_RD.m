function [] = create_actPLOT_RD(subjectID)


cd('D:\VRAWLFP\ActigraphyProcessALL\PRE_ACT_ALL')

cMAP = cividis;
% [reMAP] = reMapCmap(tmData,cMAP,smSunful,1,'timeBased');
lightCM = cMAP(246,:);
darkCM = cMAP(10,:);

[apData] = getPatDat(subjectID , 'ActALL');

% Activity
nActraw = apData.Activity;
nActrawN = cellfun(@(x) str2double(x), nActraw, 'UniformOutput',true);
sMnRmACTr = smoothdata(nActrawN,'gaussian',300,'omitnan');
nRmACTr = normalize(sMnRmACTr, 'range');
lraw = plot(nRmACTr,'LineWidth',1);
lraw.Color = 'k';
lraw.LineStyle = "-";

hold on
% Circadian Fit
nFitraw = apData.Circadian;
sMnRmFITr = smoothdata(nFitraw,'gaussian',40,'omitnan');
nRmFITr = normalize(sMnRmFITr, 'range');
lcosin = plot(nRmFITr,'LineWidth',2);
lcosin.Color = [0.5 0.5 0.5];
lcosin.LineStyle = "-.";

% Sleep Wake State
% Ronnenberg: 1 = Sleep
% Crespo: 1 = Wake
% Invert Ronneberg
reonUF = apData.Roneenberg(:);
nonNanInd1 = ~isnan(reonUF);
reonUnfurli = ~reonUF(nonNanInd1);
reonUFi = reonUF;
reonUFi(nonNanInd1) = reonUnfurli;
% Get Crespo
cresUF = apData.Crespo(:);
% Find agreement
pairRC = [reonUFi , cresUF];
nanInd2 = isnan(pairRC(:,1));
% Find nonNans
pairMatch = pairRC(:,1) == pairRC(:,2);
swFinMat = pairRC(:,1);
swFinMat(pairMatch) = pairRC(pairMatch,1);
swFinMat(~pairMatch) = nan;
swFinMat(nanInd2) = nan;

xAxisAct = 1:length(swFinMat);
p1 = plot(xAxisAct(swFinMat == 1),ones(size(xAxisAct(swFinMat == 1)))+0.05);
p1.LineStyle = "none";
p1.Color = lightCM;
p1.Marker = 'o';
p2 = plot(xAxisAct(swFinMat == 0),zeros(size(xAxisAct(swFinMat == 0)))-0.05);
p2.Color = darkCM;
p2.Marker = 'o';
p2.LineStyle = "none";
ylim([-0.1 1.1])
leg2 = legend('Raw activity','Cosinor','Wake','Sleep');
leg2.Position = [0.0490 0.6406 0.1126 0.0713];

dayStarts = round(linspace(1,length(swFinMat)-2880,length(swFinMat)/2880));
xticks(dayStarts);
xticklabels(1:length(swFinMat)/2880)
xlim([1, length(swFinMat)])
xlabel('Days of recording')
set(gca,'TickLength',[0 .001])
ylabel('Scaled activity')
yticks(linspace(0,1,3))



end







function [patDATA] = getPatDat(cID , tID)

matDir = dir('*.mat');
matNames = {matDir.name};
matEls = split(matNames,'_');

cIDs1a = matEls(:,:,1);
cIDs1b = matEls(:,:,2);
cIDs2 = cellfun(@(x,y) [x , '_' , y], cIDs1a, cIDs1b, 'UniformOutput', false);

tIDs1 = matEls(:,:,3);
tIDs2 = extractBefore(tIDs1,'.');

matLog = matches(cIDs2,cID) & matches(tIDs2,tID);
load(matNames{matLog},'rawActSlWk')
patDATA = rawActSlWk;

end


