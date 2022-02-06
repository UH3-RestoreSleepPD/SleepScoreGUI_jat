function [] = convet_oldTT_IndvTT(matLOC)
%% Duplicate table
cd(matLOC)

% Get mat list
dirMAT = dir('*.mat');
dirMATu = {dirMAT.name};

for ti = 1:length(dirMATu)

    load(dirMATu{ti},'TT')

    ttN = TT;
    ttN.CK = ttN.STNF;
    ttN.LW = ttN.STNF;
    ttN.ST = ttN.STNF;
    ttN.MS = ttN.STNF;
    % Overwrite column names
    ttN = removevars(ttN,{'STNF','UNMC'});
    % Create new columns
    ttN = movevars(ttN,{'CK','LW','ST','MS'},'After','ChinZ');
    TT = ttN;

    %% Save out
    save(dirMATu{ti},"TT");
    disp(['File ',num2str(ti), ' Done!'])

end