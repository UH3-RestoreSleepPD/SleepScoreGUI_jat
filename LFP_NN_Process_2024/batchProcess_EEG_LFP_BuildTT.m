function [] = batchProcess_EEG_LFP_BuildTT(caseID,SigNal)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% time/date csv
csvLOC = 'D:\';
cd(csvLOC)
timeDATEtab = readtable('SleepDateInfo.csv');
subTab = timeDATEtab(matches(timeDATEtab.sub,caseID),:); % CHANGE TO caseID
% caseID = 'STAN_1'
mainDIR = 'D:\VRAWLFP';
dayMatLOC = [mainDIR, filesep , 'LFP', filesep, caseID]; % CHANGE TO caseID
cd(dayMatLOC)
nightLIST = getDIR(1);

for ni = 1:length(nightLIST)
    tmpNl = [dayMatLOC , filesep , nightLIST{ni}];
    nightTab = subTab(matches(subTab.night,nightLIST{ni}),:);
    outPUTS = fileLOAD(tmpNl,subTab.sub{1});
    outPUTS.eegMODsf = 1375;
    outPUTS.lfpMODsf = 1375;
    outPUTS.recSTART = datetime(nightTab.year,nightTab.month,nightTab.day,...
        nightTab.hour,nightTab.minute,nightTab.second);

    saveLOC = [mainDIR , filesep, 'ProcLFP', filesep , caseID];
    if ~exist(saveLOC,'dir')
        mkdir(saveLOC)
    end

    subPARts = split(caseID,'_');

    outPUTS.subject = subPARts{2};
    outPUTS.Instit = subPARts{1};
    outPUTS.Night = nightLIST{ni}(end);

    selectSignalProcess(SigNal,outPUTS,saveLOC)
    disp([caseID , ' ' nightLIST{ni}, ' DONE!'])
end








% selectSignalProcess(app, 'psg')
% selectSignalProcess(app, 'lfpNF')

















end








function [outLIST] = getDIR(optionS)

switch optionS
    case 1 % folder
        dir1 = dir();
        dir2 = {dir1.name};
        dir3 = dir2(~ismember(dir2,{'.','..'}));
        outLIST = dir3;
    case 2 % mat
        dir1 = dir('*.mat');
        dir2 = {dir1.name};
        outLIST = dir2;
    case 3 % txt
        dir1 = dir('*.txt');
        dir2 = {dir1.name};
        outLIST = dir2;
end



end



function [outPUTS] = fileLOAD(folderDIR , SUBid)

cd(folderDIR)

outPUTS.fileLIST = getDIR(2);
matob = matfile(outPUTS.fileLIST{1});
outPUTS.varlist = who(matob);

if matches(SUBid,'UNMC_1')
    chanLab = {'EEG','EEG','EEG','EEG','EEG','EEG','EEG',...
        'EOG','EOG','EMG','EMG','EMG'};
    chanID = {'F3','Fz','A1','C3','Cz','A2','O1',...
        'EOG1','EOG2','Chin1','Chin2','ChinZ'};
    aoID = {'F3','Fz','A1','C3','Cz','A2','O1',...
        'CEOG_1','CEOG_2','CChin_1','CChin_2','CChin_Z'};
else

    % DETERIME ID for EEG , EOG , EMG -----------------------------
    chanLab = {'EEG','EEG','EEG','EEG','EEG','EEG','EEG','EEG','EEG',...
        'EEG','EEG','EOG','EOG','EMG','EMG','EMG'};
    chanID = {'F3','Fz','F4','A1','C3','Cz','C4','A2','O1','Oz','O2',...
        'EOG1','EOG2','Chin1','Chin2','ChinZ'};
    aoID = {'F3','Fz','F4','A1','C3','Cz','C4','A2','O1','Oz','O2',...
        'CEOG_1','CEOG_2','CChin_1','CChin_2','CChin_Z'};

end

outPUTS.labelTab = table(transpose(aoID),transpose(chanLab),...
    transpose(chanID),'VariableNames',{'AO_ID',...
    'Label','ID'});

% Convert to Time tables with appropriate sample frequencies
outPUTS.chansleep = {'EEG','EOG','EMG'};
outPUTS.chanNUMb = zeros(1,3);
for cNi = 1:length(outPUTS.chansleep)
    outPUTS.chanNUMb(cNi) = sum(ismember(outPUTS.labelTab.Label,outPUTS.chansleep{cNi}));
end

% DETERMINE ID for LFP ----------------------------------------
chanLab_LFP = {'DBS','DBS','DBS','DBS'};
chanID_LFP = {'0','1','2','3'};
aoID_LFP = {'CDBS_0','CDBS_1','CDBS_2','CDBS_3'};

outPUTS.labelTab_LFP = table(transpose(aoID_LFP),transpose(chanLab_LFP),...
    transpose(chanID_LFP),'VariableNames',{'AO_ID',...
    'Label','ID'});

% Convert to Time tables with appropriate sample frequencies
outPUTS.chansleep_LFP = {'DBS'};
outPUTS.chanNUMb_LFP = zeros(1,1);
for cNi = 1:length(outPUTS.chansleep_LFP)
    outPUTS.chanNUMb_LFP(cNi) = sum(ismember(outPUTS.labelTab_LFP.Label,outPUTS.chansleep_LFP{cNi}));
end

disp("File LOADED")


end



function selectSignalProcess(signAl,INPUTS,saveLOC)

switch signAl
    case 'psg'
        tmpEEGarray = zeros(60*60*16*INPUTS.eegMODsf,INPUTS.chanNUMb(1)); % 16 Hours
        tmpEOGarray = zeros(60*60*16*INPUTS.eegMODsf,INPUTS.chanNUMb(2)); % 16 Hours
        tmpEMGarray = zeros(60*60*16*INPUTS.eegMODsf,INPUTS.chanNUMb(3)); % 16 Hours

        chanTabLabs = cell(1,3);
        aoTabLabs = cell(1,3);

        for fi = 1:length(INPUTS.fileLIST)

            % Update progress, report current estimate
            disp(num2str(fi/length(INPUTS.fileLIST)));
            sprintf('%1.0f',fi);

            % file name
            curFile = INPUTS.fileLIST{fi};
            % Load all channels [MAY NEED TO UPDATE with MICRODRIVE]
            load(curFile,INPUTS.varlist{:});

            for ci = 1:length(INPUTS.chansleep)

                chanChoose = INPUTS.chansleep{ci};
                chanTabLabs{ci} = INPUTS.labelTab.ID(ismember(INPUTS.labelTab.Label,chanChoose));
                aoTabLabs{ci} = INPUTS.labelTab.AO_ID(ismember(INPUTS.labelTab.Label,chanChoose));

                % Extract variables relevant to channel
                % Raw data and sampling frequency
                % Get channel list
                for chi = 1:length(aoTabLabs{ci})

                    tmpChan = aoTabLabs{ci}(ismember(aoTabLabs{ci},aoTabLabs{ci}(chi)));
                    tmpRows = INPUTS.varlist(contains(INPUTS.varlist,tmpChan));
                    rawRow = tmpRows{1};

                    chanNums = 1:length(chanTabLabs{ci});

                    dat1 = transpose(eval(rawRow));
                    switch chanChoose
                        case 'EEG'
                            disp('EEG Data...');

                            lastLOC = find(tmpEEGarray(:,length(chanNums)),1,'last');

                            if isempty(lastLOC)
                                lastLOC = 0 + 1;
                            else
                                lastLOC = lastLOC + 1;
                            end

                            % trim off sample points prior to TTL
                            if fi == 1

                                timeOFFset = CDIG_IN_1_TimeBegin - CDBS_0_TimeBegin;
                                ttLSamp = round(timeOFFset*1375);
                                dat1 = dat1(ttLSamp:numel(dat1));

                            end

                            curLeng = length(dat1)-1;
                            tmpEEGarray(lastLOC:lastLOC+curLeng ,chi) = dat1;

                        case 'EOG'
                            disp('EOG Data...');

                            lastLOC = find(tmpEOGarray(:,length(chanNums)),1,'last');

                            if isempty(lastLOC)
                                lastLOC = 0 + 1;
                            else
                                lastLOC = lastLOC + 1;
                            end

                            % trim off sample points prior to TTL
                            if fi == 1

                                timeOFFset = CDIG_IN_1_TimeBegin - CDBS_0_TimeBegin;
                                ttLSamp = round(timeOFFset*1375);
                                dat1 = dat1(ttLSamp:numel(dat1));

                            end

                            curLeng = length(dat1)-1;
                            tmpEOGarray(lastLOC:lastLOC+curLeng ,chi) = dat1;

                        case 'EMG'
                            disp('EMG Data...');

                            lastLOC = find(tmpEMGarray(:,length(chanNums)),1,'last');

                            if isempty(lastLOC)
                                lastLOC = 0 + 1;
                            else
                                lastLOC = lastLOC + 1;
                            end

                            % trim off sample points prior to TTL
                            if fi == 1

                                timeOFFset = CDIG_IN_1_TimeBegin - CDBS_0_TimeBegin;
                                ttLSamp = round(timeOFFset*1375);
                                dat1 = dat1(ttLSamp:numel(dat1));

                            end

                            curLeng = length(dat1)-1;
                            tmpEMGarray(lastLOC:lastLOC+curLeng ,chi) = dat1;

                    end
                end
            end
        end
        % Close the dialog box

        % Clean up array EEG
        lastZEEG = find(tmpEEGarray(:,1),1,'last') + 1;
        tmpEEGarray = tmpEEGarray(1:lastZEEG,:);

        % Clean up array EOG
        lastZEOG = find(tmpEOGarray(:,1),1,'last') + 1;
        tmpEOGarray = tmpEOGarray(1:lastZEOG,:);

        % Clean up array EMG
        lastZEMG = find(tmpEMGarray(:,1),1,'last') + 1;
        tmpEMGarray = tmpEMGarray(1:lastZEMG,:);

        % Create Timetables -------------------------------------------

        % EEG ********************************************************
        disp('Processing EEG');

        tic;
        procEEG = struct;
        tmpEEGtab = array2table(tmpEEGarray);
        tmpEEGtab.Properties.VariableNames = chanTabLabs{1};
        procEEG.eegTTraw = table2timetable(tmpEEGtab,'SampleRate',INPUTS.eegMODsf);
        procEEG.eegProc = procEEG.eegTTraw;

        % Normalize-----DON'T normalize-----------------------------------
        % disp('Normalize');
        % procEEG.eegTTrawN = normalize(procEEG.eegTTraw);
        % Band Stop notch
        disp('Notch Filter');
        procEEG.eegNOTCH = procEEG.eegProc;
        for li = 1:size(procEEG.eegProc,2)
            tmpEEGraw = table2array(procEEG.eegProc(:,li));
            tmpEEGnotch = spectrumInterpolation(tmpEEGraw, INPUTS.eegMODsf, 60, 3, 1);
            tmpEEGntab = array2table(tmpEEGnotch);
            procEEG.eegNOTCH(:,li) = tmpEEGntab;
        end

        % High Filter
        disp('High Pass Filter');

        procEEG.eegTTppHp = procEEG.eegNOTCH;
        for li = 1:size(procEEG.eegTTppHp,2)
            eegTTpTmp = table2array(procEEG.eegTTppHp(:,li));
            procEEG.eegTTppHp{:,li} = highpass(eegTTpTmp,0.2,INPUTS.eegMODsf,'ImpulseResponse','iir','Steepness',0.8);
        end
        % Low Filter
        disp('Low Pass Filter');
        procEEG.eegTTppLp = lowpass(procEEG.eegTTppHp,100); % Only for RC+S
        procEEG.eegTTppSm = procEEG.eegTTppLp;

        % Downsample
        disp('Downsample 100Hz');
        procEEG.eegTTppDs = retime(procEEG.eegTTppSm,'regular','mean','SampleRate',100);

        % Preprocess Signal
        eegTTpp = procEEG.eegTTppDs;

        eegtoc = toc;
        disp(eegtoc)

        % EOG ********************************************************
        disp('Processing EOG');

        procEOG = struct;
        tmpEOGtab = array2table(tmpEOGarray);
        tmpEOGtab.Properties.VariableNames = chanTabLabs{2};
        procEOG.eogTTraw = table2timetable(tmpEOGtab,'SampleRate',INPUTS.eegMODsf);

        % % Normalize
        % disp('Normalize');
        % procEOG.eogTTrawN = normalize(procEOG.eogTTraw);
        % Band Stop notch
        disp('Notch Filter');

        procEOG.eogNOTCH = procEOG.eogTTraw;
        for li = 1:size(procEOG.eogTTraw,2)
            tmpEOGraw = table2array(procEOG.eogTTraw(:,li));
            tmpEOGnotch = spectrumInterpolation(tmpEOGraw, INPUTS.eegMODsf, 60, 3, 1);
            tmpEOGntab = array2table(tmpEOGnotch);
            procEOG.eogNOTCH(:,li) = tmpEOGntab;
        end

        % Low Filter
        disp('Low Pass Filter');
        procEOG.eogTTppLp = lowpass(procEOG.eogNOTCH,35);

        % High Filter
        disp('High Pass Filter');

        procEOG.eogTTppHp = procEOG.eogTTppLp;
        for li = 1:size(procEOG.eogTTppHp,2)
            eogTTpTmp = table2array(procEOG.eogTTppHp(:,li));
            procEOG.eogTTppHp{:,li} = highpass(eogTTpTmp,0.2,INPUTS.eegMODsf,'ImpulseResponse','iir','Steepness',0.8);
        end
        procEOG.eogTTppSm = procEOG.eogTTppLp;

        % Downsample
        disp('Downsample 100Hz');
        procEOG.eogTTppDs = retime(procEOG.eogTTppSm,'regular','mean','SampleRate',100);
        % Preprocess Signal
        eogTTpp = procEOG.eogTTppDs;

        eogtoc = toc;
        disp(eogtoc)

        % EMG ********************************************************
        disp('Processing EMG');

        procEMG = struct;
        tmpEMGtab = array2table(tmpEMGarray);
        tmpEMGtab.Properties.VariableNames = chanTabLabs{3};
        procEMG.emgTTraw = table2timetable(tmpEMGtab,'SampleRate',INPUTS.eegMODsf);
        % Normalize
        % disp('Normalize');
        % procEMG.emgTTrawN = normalize(procEMG.emgTTraw);
        % Band Stop notch
        disp('Notch Filter');
        procEMG.emgNOTCH = procEMG.emgTTraw;
        for li = 1:size(procEMG.emgTTraw,2)
            tmpEMGraw = table2array(procEMG.emgTTraw(:,li));
            tmpEMGnotch = spectrumInterpolation(tmpEMGraw, INPUTS.eegMODsf, 60, 3, 1);
            tmpEMGntab = array2table(tmpEMGnotch);
            procEMG.emgNOTCH(:,li) = tmpEMGntab;
        end

        % Low Filter
        disp('Low Pass Filter');
        procEMG.emgTTppLp = lowpass(procEMG.emgNOTCH,100);
        % High Filter
        disp('High Pass Filter');

        procEMG.emgTTppHp = procEMG.emgTTppLp;
        for li = 1:size(procEMG.emgTTppHp,2)
            emgTTpTmp = table2array(procEMG.emgTTppHp(:,li));
            procEMG.emgTTppHp{:,li} = highpass(emgTTpTmp,10,INPUTS.eegMODsf,'ImpulseResponse','iir','Steepness',0.8);
        end
        procEMG.emgTTppSm = procEMG.emgTTppHp;

        % Downsample
        disp('Downsample 100Hz');
        procEMG.emgTTppDs = retime(procEMG.emgTTppSm,'regular','mean','SampleRate',100);
        % Preprocess Signal
        emgTTpp = procEMG.emgTTppDs;

        emgtoc = toc;
        disp(emgtoc)

        % Create FINAL Timetables -------------------------------------
        eegTABf = createBINS(eegTTpp,'eeg',chanTabLabs{1},100);
        eogTABf = createBINS(eogTTpp,'eog',chanTabLabs{2},100);
        emgTABf = createBINS(emgTTpp,'emg',chanTabLabs{3},100);

        % Create FINAL EEG/EOG/EMG Timetables -------------------------
        TT = [eegTABf,eogTABf,emgTABf];
        % Add Sleep Score raters
        TT.LW = repmat({''},height(TT),1);
        TT.CK = repmat({''},height(TT),1);
        TT.ST = repmat({''},height(TT),1);
        TT.MS = repmat({''},height(TT),1);

        if isempty(INPUTS.recSTART)
            addTval = 0;
        else
            addTval = INPUTS.recSTART;
        end

        actimeVec = transpose(0:seconds(30):seconds((height(TT)-1)*30)) + addTval;
        TT.ActTime = actimeVec;
        % SAVE data
        cd(saveLOC)
        fname = [INPUTS.subject,'_',INPUTS.Instit,'_',INPUTS.Night];
        fnameMAT = [fname,'_PSG.mat'];
        save(fnameMAT,'TT')

    case 'lfp'
        % CD back to raw data
        cd(app.fileLOC)


        app.tasK.Text = 'FULL';

        progBAR = uiprogressdlg(app.UIFigure,'Title','Combining Files...',...
            'Message','1');
        drawnow

        % LFP
        tmpLFParray = zeros(60*60*16*app.lfpMODsf,app.chanNUMb_LFP(1)); % 14 Hours
        chanTabLabs_LFP = cell(1,1);
        aoTabLabs_LFP = cell(1,1);

        % LFP Processing
        for fi = 1:length(app.fileLIST)
            % Update progress, report current estimate
            progBAR.Value = fi/length(app.fileLIST);
            progBAR.Message = sprintf('%1.0f',fi);
            progBAR.Title = 'LFP Data...';
            % file name
            curFile = app.fileLIST{fi};
            % Load all channels
            load(curFile,app.varlist{:});

            for ci = 1:length(app.chansleep_LFP)

                chanChooselfp = app.chansleep_LFP{ci};
                chanTabLabs_LFP{ci} = app.labelTab_LFP.ID(ismember(app.labelTab_LFP.Label,chanChooselfp));
                aoTabLabs_LFP{ci} = app.labelTab_LFP.AO_ID(ismember(app.labelTab_LFP.Label,chanChooselfp));

                for chi = 1:length(aoTabLabs_LFP{ci})

                    tmpChanlfp = aoTabLabs_LFP{ci}(ismember(aoTabLabs_LFP{ci},aoTabLabs_LFP{ci}(chi)));
                    tmpRowslfp = app.varlist(contains(app.varlist,tmpChanlfp));
                    rawRowlfp = tmpRowslfp{1};

                    chanNumslfp = 1:length(chanTabLabs_LFP{ci});

                    lfptemp = transpose(eval(rawRowlfp));

                    lastLOClfp = find(tmpLFParray(:,length(chanNumslfp)),1,'last');

                    if isempty(lastLOClfp)
                        lastLOClfp = 0 + 1;
                    else
                        lastLOClfp = lastLOClfp + 1;
                    end

                    % trim off sample points prior to TTL
                    if fi == 1

                        timeOFFset = CDIG_IN_1_TimeBegin - CDBS_0_TimeBegin;
                        ttLSamp = round(timeOFFset*1375);
                        lfptemp = lfptemp(ttLSamp:numel(lfptemp));

                    end

                    curLeng = length(lfptemp)-1;
                    tmpLFParray(lastLOClfp:lastLOClfp+curLeng ,chi) = lfptemp;
                end
            end
        end

        % Close the dialog box
        close(progBAR)

        % Clean up array LFP
        lastZLFP = find(tmpLFParray(:,1),1,'last') + 1;
        tmpLFParray = tmpLFParray(1:lastZLFP,:);

        % LFP ********************************************************
        lfpPB = uiprogressdlg(app.UIFigure,'Title','Processing LFP',...
            'Indeterminate','on');
        drawnow

        tmpLFPtab = array2table(tmpLFParray);
        tmpLFPtab.Properties.VariableNames = chanTabLabs_LFP{1};
        app.lfpTTraw = table2timetable(tmpLFPtab,'SampleRate',app.lfpMODsf);

        % Normalize
        lfpPB.Message = sprintf('%s','Normalize');
        lfpTTrawN = normalize(app.lfpTTraw);
        % Band Stop notch
        lfpPB.Message = sprintf('%s','Notch Filter');
        lfpNOTCH = lfpTTrawN;
        for li = 1:size(lfpTTrawN,2)
            tmpLFPraw = table2array(lfpTTrawN(:,li));
            tmpLFPnotch = spectrumInterpolation(tmpLFPraw, app.lfpMODsf, 60, 3, 1);
            tmpLFPntab = array2table(tmpLFPnotch);
            lfpNOTCH(:,li) = tmpLFPntab;
        end

        % High Filter
        lfpPB.Message = sprintf('%s','High Pass Filter');
        %                     lfpTTppHp = highpass(lfpNOTCH,0.2); % Only for RC+S

        lfpTTppHp = lfpNOTCH;
        for li = 1:size(lfpTTppHp,2)
            lfpTTpTmp = table2array(lfpTTppHp(:,li));
            lfpTTppHp{:,li} = highpass(lfpTTpTmp,0.2,1375,'ImpulseResponse','iir','Steepness',0.8);
        end

        % Low Filter
        lfpPB.Message = sprintf('%s','Low Pass Filter');
        lfpTTppLp = lowpass(lfpTTppHp,250); % Only for RC+S
        lfpTTppSm = lfpTTppLp;

        % Downsample
        % CREATE RC+S version
        % CREATE RAW version - no downsampling
        lfpPB.Message = sprintf('%s','Downsample 250Hz');
        lfpTTppDs = retime(lfpTTppSm,'regular','mean','SampleRate',250);
        % Preprocess Signal
        app.lfpTTpp = lfpTTppDs;
        app.lfpTTNdnS = lfpTTppHp;

        close(lfpPB)

        lfpTABf = createBINS(app,'lfp',chanTabLabs_LFP{1},250);
        lfpTABfND = createBINS(app,'lfpr',chanTabLabs_LFP{1},1375);

        % Add ActTime % USE AO CONVERTER TO GET START TIME OF FIRST
        % FILE

        % Create FINAL LFP Timetables -------------------------
        LFPTT = lfpTABf;

        actimeVec = transpose(0:seconds(30):seconds((height(LFPTT)-1)*30)) + app.recSTART;

        LFPTT.ActTime = actimeVec;
        LFPTTRaw = lfpTABfND;
        LFPTTRaw.ActTime = actimeVec;
        % SAVE OUT
        cd(app.saveLOC)

        app.fname = [app.outP{1},'_',app.outP{2},'_',app.outP{3}];
        fnameMATlfp = [app.fname,'_LFP.mat'];
        fnameMATlfpR = [app.fname,'_LFPraw.mat'];

        save(fnameMATlfp,'LFPTT')
        save(fnameMATlfpR,'LFPTTRaw','-v7.3')

        app.INFOLabel.Text = "DONE";
        app.tasK.Text = 'DONE';

    case 'lfpNF'
        % CD back to raw data

        disp('Combining Files');

        % LFP
        tmpLFParray = zeros(60*60*16*INPUTS.lfpMODsf,INPUTS.chanNUMb_LFP(1)); % 14 Hours
        chanTabLabs_LFP = cell(1,1);
        aoTabLabs_LFP = cell(1,1);

        % LFP Processing
        for fi = 1:length(INPUTS.fileLIST)
            % Update progress, report current estimate
            disp('LFP Data...');
            % file name
            curFile = INPUTS.fileLIST{fi};
            % Load all channels
            load(curFile,INPUTS.varlist{:});

            for ci = 1:length(INPUTS.chansleep_LFP)

                chanChooselfp = INPUTS.chansleep_LFP{ci};
                chanTabLabs_LFP{ci} = INPUTS.labelTab_LFP.ID(ismember(INPUTS.labelTab_LFP.Label,chanChooselfp));
                aoTabLabs_LFP{ci} = INPUTS.labelTab_LFP.AO_ID(ismember(INPUTS.labelTab_LFP.Label,chanChooselfp));

                for chi = 1:length(aoTabLabs_LFP{ci})

                    tmpChanlfp = aoTabLabs_LFP{ci}(ismember(aoTabLabs_LFP{ci},aoTabLabs_LFP{ci}(chi)));
                    tmpRowslfp = INPUTS.varlist(contains(INPUTS.varlist,tmpChanlfp));
                    rawRowlfp = tmpRowslfp{1};

                    chanNumslfp = 1:length(chanTabLabs_LFP{ci});

                    lfptemp = transpose(eval(rawRowlfp));

                    lastLOClfp = find(tmpLFParray(:,length(chanNumslfp)),1,'last');

                    if isempty(lastLOClfp)
                        lastLOClfp = 0 + 1;
                    else
                        lastLOClfp = lastLOClfp + 1;
                    end

                    % trim off sample points prior to TTL
                    if fi == 1

                        timeOFFset = CDIG_IN_1_TimeBegin - CDBS_0_TimeBegin;
                        ttLSamp = round(timeOFFset*1375);
                        lfptemp = lfptemp(ttLSamp:numel(lfptemp));

                    end

                    curLeng = length(lfptemp)-1;
                    tmpLFParray(lastLOClfp:lastLOClfp+curLeng ,chi) = lfptemp;
                end
            end
        end

        % Clean up array LFP
        lastZLFP = find(tmpLFParray(:,1),1,'last') + 1;
        tmpLFParray = tmpLFParray(1:lastZLFP,:);

        % LFP ********************************************************
        disp('Processing LFPnf');

        procLFPnf = struct;
        tmpLFPtab = array2table(tmpLFParray);
        tmpLFPtab.Properties.VariableNames = chanTabLabs_LFP{1};
        procLFPnf.lfpTTrawNF = table2timetable(tmpLFPtab,'SampleRate',INPUTS.lfpMODsf);
        procLFPnf.lfpTTNdnS = procLFPnf.lfpTTrawNF;

        lfpNFTABf = procLFPnf.lfpTTNdnS;

        lfpTABfND = createBINS(lfpNFTABf,'lfpr',chanTabLabs_LFP{1},1375);
        LFPTTRaw = lfpTABfND;

        if isempty(INPUTS.recSTART)
            addTval = 0;
        else
            addTval = INPUTS.recSTART;
        end

        % actimeVec = transpose(0:seconds(30):seconds((height(TT)-1)*30)) + addTval;

        actimeVec = transpose(0:seconds(30):seconds((height(LFPTTRaw)-1)*30)) + addTval;
        LFPTTRaw.ActTime = actimeVec;

        % SAVE data
        cd(saveLOC)

        fname = [INPUTS.subject,'_',INPUTS.Instit,'_',INPUTS.Night];
        fnameMATlfpRnf = [fname,'_LFPrawNF.mat'];

        save(fnameMATlfpRnf,'LFPTTRaw','-v7.3')
        disp("DONE!");




end

end




function outTTab = createBINS(rawDAT,chan,chanLABS,SF)

% switch chan
%     case 'emg'
%         rawDAT = app.emgTTpp;
%     case 'eog'
%         rawDAT = app.eogTTpp;
%     case 'eeg'
%         rawDAT = app.eegTTpp;
%     case 'lfp'
%         rawDAT = app.lfpTTpp;
%     case 'lfpr'
%         rawDAT = app.lfpTTNdnS;
% 
% end

disp(['Creating Bins ',chan]);

num30bins = floor(height(rawDAT)/SF/30);
tstart = seconds(0);
tend = seconds(30);
outBINS = cell(num30bins,size(rawDAT,2));
for si = 1:num30bins

    timeLimits = timerange(tstart,tend);

    tmpT = rawDAT(timeLimits,:);
    tt2t = timetable2table(tmpT);
    t2c = table2cell(tt2t);
    cleanCel = cell2mat(t2c(:,2:end));

    for cchi = 1:size(rawDAT,2)
        outBINS{si,cchi} = cleanCel(:,cchi);
    end

    tstart = tend;
    tend = tstart + seconds(30);

    disp([num2str(si) ,' out of ' , num2str(num30bins)])

end

% Create new TimeTable
% Time
timeBINS = transpose(seconds(0:30:30*num30bins-1));
% Data
finalTT = cell2table(outBINS);
% Variable names
varNAMES = chanLABS;
% Create Table
finalTT.Properties.VariableNames = varNAMES;
finalTT.Time = timeBINS;

outTTab = table2timetable(finalTT);

end


