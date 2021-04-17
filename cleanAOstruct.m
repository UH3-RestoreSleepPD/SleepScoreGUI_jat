%% Load
cd('G:\01 Coding Datasets\SLEEP_SCORE_EDF\RAW_UNMC_02_N3\old')

raw = 'G:\01 Coding Datasets\SLEEP_SCORE_EDF\RAW_UNMC_02_N3\old';
new = 'G:\01 Coding Datasets\SLEEP_SCORE_EDF\RAW_UNMC_02_N3\new';

mdir = dir('*.mat');
mdir2 = {mdir.name};
fileLIST = mdir2;


%% Load



for fi = 1:length(fileLIST)

    matob = matfile(fileLIST{fi});
    varlist = who(matob);
    chanlist = varlist(~contains(varlist,{'SF','Ports','CADD','Channel','CANALOG'}));
    
    cd(raw)
    load(fileLIST{fi},chanlist{:});
    cd(new)
    save(fileLIST{fi},chanlist{:});
    
    
end


%% Check ttl

new = 'C:\Users\John\Desktop\sleepQCcheck\samplechck2';
cd(new)
mdir = dir('*.mat');
mdir2 = {mdir.name};
fileLIST = mdir2;

ttlC = cell(1,length(fileLIST))
ttlA = cell(1,length(fileLIST))
for fi = 1:length(fileLIST)
    
    
    matob = matfile(fileLIST{fi});
    varlist = who(matob);
    chanlist = varlist(contains(varlist,{'CDIG'}));
    
    if isempty(chanlist)
       
        continue
        
    else
        
        load(fileLIST{fi},'CDIG_IN_1_Up');
        
        ttlA{fi} = CDIG_IN_1_Up;
        
        dt  = diff(CDIG_IN_1_Up / 44000);
        
        
        ttlC{fi} = dt;
        
        clearvars CDIG_IN_1_Up
    end
    
    
    
    
    
    
    
end

