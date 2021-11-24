cd('I:\01_Coding_Datasets\TEST_FILE_sleep');
load("TEST_Sleep.mat",'TT');
scoreD = TT(:,{'LW','ST','CK','MS'});
%%

disLog = false(height(scoreD),1);
for s = 1:height(scoreD)

    tmpRow = table2cell(scoreD(s,:));
    % Determine if there are any empty cells
    empCk = any(cellfun(@(x) isempty(x), tmpRow, 'UniformOutput',true));
    % Get unique vals
    uVs = numel(unique(tmpRow));
    if empCk
        disLog(s) = true;
    elseif uVs > 1
        disLog(s) = true;
    end
end

%%
% Order of sleep IDs
sleePids = {'W','N1','N2','N3','R'};
plotNUMs = nan(height(scoreD),5);
scoreVals = table2cell(scoreD(:,1));
allScores = table2cell(scoreD);
for pi = 1:size(plotNUMs,2)
    for si = 1:length(sleePids)

        switch pi
            case 1

                logLoc = ismember(scoreVals,sleePids{si}) & ~disLog;
                plotNUMs(logLoc,pi) = si;

            case {2,3,4,5}

                logLoc = ismember(allScores(:,pi-1),sleePids{si}) & disLog;
                plotNUMs(logLoc,pi) = si;

        end

    end
end

%%
xAx = 1:height(scoreD);

%%
% Plot agreed upon points
scatter(xAx, plotNUMs(:,1), 20, 'MarkerFaceColor','none','MarkerEdgeColor',...
    [0.5 0.5 0.5],'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
ylim([0 6])


epochof = 25;
idInds = {'LW','ST','CK','MS'};
xline(epochof,'LineStyle','-','LineWidth',3,'Color','g')

epocVals = plotNUMs(epochof,2:5);
uniScores = unique(epocVals);

for u = 1:length(uniScores)

    numIndsS = idInds(epocVals == uniScores(u));

    if length(numIndsS) == 1
        text(epochof,uniScores(u),numIndsS{uiu})

    elseif length(numIndsS) == 2
        for uiu = 1:length(numIndsS)
            if uiu == 1
                text(epochof,uniScores(u),numIndsS{uiu},...
                    'HorizontalAlignment','left');
            else
                text(epochof,uniScores(u),numIndsS{uiu},...
                    'HorizontalAlignment','right');
            end
        end
    else

    end

end



%%

% Plot discordant points
hold on
scatter(xAx, plotNUMs(:,2), 20, 'MarkerFaceColor','none','MarkerEdgeColor',...
    [1 0 0],'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)

scatter(xAx, plotNUMs(:,3), 20, 'MarkerFaceColor','none','MarkerEdgeColor',...
    [0 1 0],'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)

scatter(xAx, plotNUMs(:,4), 20, 'MarkerFaceColor','none','MarkerEdgeColor',...
    [0 0 1],'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)

scatter(xAx, plotNUMs(:,5), 20, 'MarkerFaceColor','none','MarkerEdgeColor',...
    [1 0.5 0.5],'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)



%% Line plot
sleePids = {'W','N1','N2','N3','R'};
plotNUMsL = nan(height(scoreD),4);
allScores = table2cell(scoreD);
for pi = 1:size(plotNUMsL,2)
    for si = 1:length(sleePids)

        logLoc = ismember(allScores(:,pi),sleePids{si});
        plotNUMsL(logLoc,pi) = si;

    end
end

cOlOrs = 'rgbk';
hold on
for iL = 1:4
    plot(xAx,plotNUMsL(:,iL),"Color",cOlOrs(iL))
end
