%preprocess data for each patient
%timetable structure is LFPTT.(var_num){epoch,1}
%var_num in [1,4] is channel
%var_num 5 is sleep label

mainRawLoc = 'E:\Dropbox\Publications_Meta\InProgress\Joel_Z_NREM_update_2022\Raw_MinnData';
mainSaveLoc = 'E:\Dropbox\Publications_Meta\InProgress\Joel_Z_NREM_update_2022\JZ_MinnMats';

best_bipolar = [2 2 3 3 1 1 2 3 2 1]; %from Dulce

numberID = 1;
hertzID = 250;
folderLOC = transpose({'250Hz'});

infoTABLE = table(numberID , hertzID , folderLOC,'VariableNames',...
    {'NumID','HerTZ','FoldLoc'});

dsiis = 250;

tabRow = infoTABLE(1,:);
cd([mainRawLoc , filesep , tabRow.FoldLoc{1}])

LFPPow = cell(1,10);
for ptnum = 1:10
    cd([mainRawLoc , filesep , tabRow.FoldLoc{1}])

    try
        fname = strcat(num2str(ptnum),'_UMin_1_LFPraw.mat');
        load(fname)
    catch
        fname = strcat(num2str(ptnum),'_UMin_1_LFP.mat');
        load(fname)
    end

    TT = LFPTT;

    num_epochs = length(TT.(1)(:,1));
    epochLeng = numel(TT.(1){1});
    %do the power spectral analysis for each epoch

    band_pow = zeros(3,num_epochs,8);
    for epoch=1:num_epochs
        bp = zeros(3,epochLeng);
        for bipolar = 1:3
            %subtract neighboring channels to get bipolar derivations
            bp(bipolar,:) = TT.(bipolar){epoch,1}(1:end) - TT.(bipolar+1){epoch,1}(1:end);
        end
        %calculate pwelch
        numFFTbins = 2^nextpow2(dsiis) + 1;
        % powBINS = floor(numFFTbins/2) + 1;
        pow2use = zeros(3,numFFTbins);
        for bipolar = 1:3
            [pow2use(bipolar,:) , ffreq] = pwelch(bp(bipolar,:),dsiis*2,[],[],dsiis); %2s window is 500 time points at 250 Hz
        end

        deltaIND  = ffreq <= 3;
        thetaIND  = ffreq <= 7  & ffreq >= 3;
        alphaIND  = ffreq <= 13 & ffreq >= 7;
        lbetaIND  = ffreq <= 21 & ffreq >= 13;
        hbetaIND  = ffreq <= 30 & ffreq >= 21;
        lgammaIND = ffreq <= 50 & ffreq >= 30;
        hgammaIND = ffreq <= 90 & ffreq >= 50;
        hfoIND    = ffreq >= 90;

        freqBANDinds = {deltaIND , thetaIND , alphaIND , lbetaIND,...
            hbetaIND , lgammaIND , hgammaIND , hfoIND};

        %remove DC and 60 Hz
        % Find index up to 1Hz
        dcOUT = ffreq <= 1;
        pow2use(:,dcOUT) = nan; %filter out DC, up to 1Hz = 2 cycles / 2s (freqs output in pwelch are 0,1,2,... cycles/2s)
        lineNOISE = ffreq >= 59 & ffreq <= 61;
        pow2use(:,lineNOISE) = nan; %filter around 60Hz

        for bipolar = 1:3
            for bandI = 1:8
                logBANDind = freqBANDinds{bandI};
                band_pow(bipolar,epoch,bandI) = mean(pow2use(bipolar,logBANDind),'omitnan');
            end
        end
        %print so I can see how it's progressing
        if mod(epoch,25)==1
            disp(epoch)
        end
    end
    %now go through and remove outliers
    for bipolar = 1:3
        for band = 1:8
            powsALL = band_pow(bipolar,:,band);
            outliers = find(powsALL > 5*median(powsALL,'omitnan'));
            powsALL(outliers) = 5*median(powsALL,'omitnan');
            band_pow(bipolar,:,band) = powsALL;
        end
    end

    %so output LFPPpow{ptnum} = squeeze(band_pow(best_bipolar,:,:));
    tmpReduceMat = squeeze(band_pow(best_bipolar(ptnum),:,:));

    if isnan(tmpReduceMat(1,8))
        finalMAT = tmpReduceMat(:,1:7);
    else
        finalMAT = tmpReduceMat;
    end

    LFPPow{ptnum} = finalMAT;
    disp(ptnum)
end
% save LFPPow
% CD to correct location
cd([mainSaveLoc , filesep , tabRow.FoldLoc{1}])

% Use patient save name
saveNAME = ['LFPPow_',tabRow.FoldLoc{1},'.mat'];

% Use LFPPow for all
save(saveNAME , 'LFPPow');











