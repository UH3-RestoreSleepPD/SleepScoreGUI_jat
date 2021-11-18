%% Duplicate table
ttN = TT;
% Overwrite column names
ttNr = renamevars(ttN,["STNF","UNMC"], ["LW","ST"]);
% Create new columns
CK = ttNr.LW;
MS = ttNr.ST;
% Modify some of the values
uniVals = unique(MS);
uniVals2 = uniVals(cellfun(@(x) ~isempty(x), uniVals, 'UniformOutput', true));
twnRv1 = randperm(height(MS),15);
twnRv2 = randperm(height(MS),15);
seq1 = datasample(uniVals2,15);
seq2 = datasample(uniVals2,15);

CK(twnRv1) = seq1;
MS(twnRv2) = seq2;

%% Add new columns
TT = addvars(ttNr,CK,'After','ST');
TT = addvars(TT,MS,'After','CK');

%% Save out
save("TEST_Sleep.mat","TT");