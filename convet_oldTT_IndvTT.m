function [] = convet_oldTT_IndvTT(matLOC)
%% Duplicate table
cd(matLOC)

% Get mat list
dirMAT = dir('*.mat');
dirMATu = {dirMAT.name};

for ti = 1:length(dirMATu)

    load(dirMATu{ti},'TT')

    ttN = TT;

    tmpEmp = repmat({''}, height(TT),1);

    ttN.CK = tmpEmp;
    ttN.LW = tmpEmp;
    ttN.ST = tmpEmp;
    ttN.MS = tmpEmp;

    % Overwrite column names
    tabNames = ttN.Properties.VariableNames;
    if matches('STNF',tabNames)
        ttN = removevars(ttN,{'STNF','UNMC'});
        % Create new columns
        ttN = movevars(ttN,{'CK','LW','ST','MS'},'After','ChinZ');
        TT = ttN;
    else
        TT = ttN;
    end

    %% Save out
    save(dirMATu{ti},"TT");
    disp(['File ',num2str(ti), ' Done!'])

end