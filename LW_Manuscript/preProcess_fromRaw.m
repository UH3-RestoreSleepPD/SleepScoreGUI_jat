clear
%% Github loc
gitloc = 'C:\Users\Admin\Documents\Github\SleepScoreGUI_jat\LW_Manuscript';
cd(gitloc)

%% initial data loc
todo = 'f';
switch todo
    case 'i'

        curLoc = 'H:\SleepStudy2_Aim1\vrawLFP_initialScore';
        curFile = 'H:\SleepStudy2_Aim1\vrawLFP_initialScore\InitialReview.mat';

        % final data loc
    case 'f'
        curLoc = 'H:\SleepStudy2_Aim1\vrawLFP_finalScore';
        curFile = 'H:\SleepStudy2_Aim1\vrawLFP_finalScore\finalReview.mat';
end

load(curFile)
%% process initial current loc

cd(curLoc)

folddir = dir();
folddir1 = {folddir.name};
folddir2 = folddir1([folddir.isdir]);
folddir3 = folddir2(~ismember(folddir2,{'..','.'}));

startIND = length(finalNlist) + 1;

for fi = 1:length(folddir3)

    tmpLoc = [curLoc , filesep , folddir3{fi}];
    cd(tmpLoc)

    matdir = dir('*.mat');
    matdir1 = {matdir.name};

    for mi = 1:length(matdir1)

        mitemp = matdir1{mi};
        load(mitemp,'TT')

        tmpTT = TT(:,17:20);

        finalConfL{1,startIND} = tmpTT;
        finalNlist{1,startIND} = matdir1{1};

        startIND = startIND + 1;


    end

end

%% process final current loc 


cd(curLoc)

folddir = dir();
folddir1 = {folddir.name};
folddir2 = folddir1([folddir.isdir]);
folddir3 = folddir2(~ismember(folddir2,{'..','.'}));

startIND = length(finalNlist) + 1;

for fi = 1:length(folddir3)

    tmpLoc = [curLoc , filesep , folddir3{fi}];
    cd(tmpLoc)

    folddirB = dir();
    folddirB1 = {folddirB.name};
    folddirB3 = folddirB1(~ismember(folddirB1,{'..','.'}));

    for ffi = 1:length(folddirB3)
        tmpLoc2 = [tmpLoc , filesep , folddirB3{ffi}];
        cd(tmpLoc2)

        matdir = dir('*.mat');
        matdir1 = {matdir.name};

        for mi = 1:length(matdir1)

            mitemp = matdir1{mi};
            load(mitemp,'TT')

            tmpTT = TT.FINALSCORE;
            tmpSpli = split(folddirB3{ffi},'_');
            scFname = tmpSpli(1);

    

        end

        finalConfL{1,startIND} = tmpTT;
        finalNlist{1,startIND} = matdir1{1};

        startIND = startIND + 1;



    end





end

%%
% Day ID = 'UPEN_1_1
rowID = 39;
finalNlist(rowID)
%%

userID = 'MS';

fxCol = TT.FINALSCORE;

sum(cellfun(@(x) isempty(x), fxCol))

finalConfL{1,rowID}.(userID) = fxCol;

%%

up3_1 = TT(:,17:20);
up3_1.LW = [];
up3_1.CK = [];
up3_1.MS = TT.FINALSCORE;

%%

up3_2 = TT(:,17:20);
up3_2.LW = [];
up3_2.CK = [];
up3_2.MS = TT.FINALSCORE;

%%

up3_3 = TT(:,17:20);
up3_3.LW = [];
up3_3.CK = [];
up3_3.MS = TT.FINALSCORE;

%%
up3_1.ST = TT.FINALSCORE;

%%
up3_2.ST = TT.FINALSCORE;

%%
up3_3.ST = TT.FINALSCORE;

%% 

finalConfL{43} = up3_1;
finalConfL{44} = up3_2;
finalConfL{45} = up3_3;
finalNlist{43} = 'UPEN_3_1';
finalNlist{44} = 'UPEN_3_2';
finalNlist{45} = 'UPEN_3_3';


%% 

todo = 'i';
switch todo

    case 'i'
        datAOI = initialDat;
        datEOI = allLIST;

    case 'f'
        datAOI = finalConfL;
        datEOI = finalNlist;
end

cAll = nan(200,1);
cCount = 1;
cName = cell(200,1);
scName = cell(200,1);
for fii = 1:length(datAOI)

    tmpC = datAOI{fii};

    numC = width(tmpC);

    for ni = 1:numC
        tmpN = table2cell(tmpC(:,ni));
        fracM = (sum(cellfun(@(x) isempty(x), tmpN)) / length(tmpN))*100;
        cAll(cCount) = fracM;
        cName{cCount} = datEOI{fii};
        tmpNN = tmpC(:,ni);
        tmpNN1 = tmpNN.Properties.VariableNames;
        scName{cCount} = tmpNN1{1};
        cCount = cCount + 1;
    end
end

cAll = cAll(~isnan(cAll));
cName = cName(cellfun(@(x) ~isempty(x), cName));
scName = scName(cellfun(@(x) ~isempty(x), scName));

outTable = table(cName, scName, cAll, 'VariableNames',{'IS_Night','Scorer','PercMiss'});

writetable(outTable,'PercentMissingIN.csv');








