function [] = fixUNMC_1(datDIR,newDIR)

% 1 – F3
% 2 – Fz
% 3 – Depth 0 (black)
% 4 – A1 (M1)
% 5 – C3
% 6 – Cz
% 7 – Depth 1 (red)
% 8 – A2  (M2)
% 9 – O1
% 10 – Depth 2 (yellow)
% 11 – EOG1
% 12 – EOG2
% 13 – Chin1
% 14 – Chin2
% 15 - ChinZ
% 16 – Depth 3 (yellow)
cd(datDIR)
matDir1 = dir('*.mat');
matDir2 = {matDir1.name};

for mi = 1:length(matDir2)
    cd(datDIR)

    % get file names
    tmMat = matDir2{mi};
    tmMatf = matfile(tmMat);
    whoTm = who(tmMatf);

    % Get EEG rows
    load(tmMat)

    cases = {'CDBS','EOG','EMG'};
    % EMG = CChin_1 , CChin_2, CChin_Z
    % CChin_2_TimeEnd
    % CDBS = CDBS_0, CDBS_1, CDBS_2, CDBS_3
    % CDBS_0_TimeEnd
    % EOG = CEOG_1, CEOG_2
    % CEOG_1_KHz
    for ci = 1:length(cases)
        tmpCASE = cases{ci};
        switch tmpCASE
            case 'CDBS'
                % 3 7 10 16

                eegAll = whoTm(contains(whoTm,'CEEG'));
                extNums = extractAfter(extractBetween(eegAll,'__','__'),'_');
                conNums = cellfun(@(x) str2double(x), extNums,'UniformOutput',true);
                dbsNums = conNums(ismember(conNums,[3 7 10 16]));
                dbsRows = eegAll(ismember(conNums,[3 7 10 16]));

                for dbi = 1:length(dbsRows)

                    tmpNum = dbsNums(dbi);
                    tmpName = dbsRows{dbi};

                    switch tmpNum
                        case 3
                            if length(tmpName) == 24
                                newName = 'CDBS_0';
                                localfcn(newName,eval(tmpName))
                                clearvars(tmpName)
                            else
                                newName = ['CDBS_0',extractAfter(tmpName,24)];
                                localfcn(newName,eval(tmpName))
                                clearvars(tmpName)
                            end

                        case 7
                            if length(tmpName) == 24
                                newName = 'CDBS_1';
                                localfcn(newName,eval(tmpName))
                                clearvars(tmpName)
                            else
                                newName = ['CDBS_1',extractAfter(tmpName,24)];
                                localfcn(newName,eval(tmpName))
                                clearvars(tmpName)
                            end

                        case 10
                            if length(tmpName) == 24
                                newName = 'CDBS_2';
                                localfcn(newName,eval(tmpName))
                                clearvars(tmpName)
                            else
                                newName = ['CDBS_2',extractAfter(tmpName,24)];
                                localfcn(newName,eval(tmpName))
                                clearvars(tmpName)
                            end


                        case 16
                            if length(tmpName) == 24
                                newName = 'CDBS_3';
                                localfcn(newName,eval(tmpName))
                                clearvars(tmpName)
                            else
                                newName = ['CDBS_3',extractAfter(tmpName,24)];
                                localfcn(newName,eval(tmpName))
                                clearvars(tmpName)
                            end

                    end



                end

            case 'EOG'

                eegAll = whoTm(contains(whoTm,'CEEG'));
                extNums = extractAfter(extractBetween(eegAll,'__','__'),'_');
                conNums = cellfun(@(x) str2double(x), extNums,'UniformOutput',true);
                eogNums = conNums(ismember(conNums,[11 12]));
                eogRows = eegAll(ismember(conNums,[11 12]));

                for eogi = 1:length(eogRows)

                    tmpNum = eogNums(eogi);
                    tmpName = eogRows{eogi};

                    switch tmpNum
                        case 11
                            if length(tmpName) == 24
                                newName = 'CEOG_1';
                                localfcn(newName,eval(tmpName))
                                clearvars(tmpName)
                            else
                                newName = ['CEOG_1',extractAfter(tmpName,24)];
                                localfcn(newName,eval(tmpName))
                                clearvars(tmpName)
                            end

                        case 12
                            if length(tmpName) == 24
                                newName = 'CEOG_2';
                                localfcn(newName,eval(tmpName))
                                clearvars(tmpName)
                            else
                                newName = ['CEOG_2',extractAfter(tmpName,24)];
                                localfcn(newName,eval(tmpName))
                                clearvars(tmpName)
                            end


                    end


                end


            case 'EMG'

                eegAll = whoTm(contains(whoTm,'CEEG'));
                extNums = extractAfter(extractBetween(eegAll,'__','__'),'_');
                conNums = cellfun(@(x) str2double(x), extNums,'UniformOutput',true);
                emgNums = conNums(ismember(conNums,[13 14 15]));
                emgRows = eegAll(ismember(conNums,[13 14 15]));


                for emgi = 1:length(emgRows)

                    tmpNum = emgNums(emgi);
                    tmpName = emgRows{emgi};

                    switch tmpNum
                        case 13
                            if length(tmpName) == 24
                                newName = 'CChin_1';
                                localfcn(newName,eval(tmpName))
                                clearvars(tmpName)
                            else
                                newName = ['CChin_1',extractAfter(tmpName,24)];
                                localfcn(newName,eval(tmpName))
                                clearvars(tmpName)
                            end

                        case 14
                            if length(tmpName) == 24
                                newName = 'CChin_2';
                                localfcn(newName,eval(tmpName))
                                clearvars(tmpName)
                            else
                                newName = ['CChin_2',extractAfter(tmpName,24)];
                                localfcn(newName,eval(tmpName))
                                clearvars(tmpName)
                            end

                        case 15
                            if length(tmpName) == 24
                                newName = 'CChin_Z';
                                localfcn(newName,eval(tmpName))
                                clearvars(tmpName)
                            else
                                newName = ['CChin_Z',extractAfter(tmpName,24)];
                                localfcn(newName,eval(tmpName))
                                clearvars(tmpName)
                            end


                    end


                end



        end

    end

    newWHO = whos;
    newWHO2 = {newWHO.name};
    toSAVE = newWHO2(contains(newWHO2,{'CChin','CDBS','CDIG','CEEG','CEOG'}));
    cd(newDIR);
    save(tmMat , toSAVE{:});
    clearvars(toSAVE{:})

end



end





function localfcn(in,out)
assignin('caller',in,out)
end