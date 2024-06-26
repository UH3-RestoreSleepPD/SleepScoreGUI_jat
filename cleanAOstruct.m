%% Load
cd('C:\Users\Admin\Downloads\oldMW2')

raw = 'C:\Users\Admin\Downloads\oldMW2';
new = 'C:\Users\Admin\Downloads\newMW2';

mdir = dir('*.mat');
mdir2 = {mdir.name};
fileLIST = mdir2;


% Load



for fi = 1:length(fileLIST)
    cd(raw)
    matob = matfile(fileLIST{fi});
    varlist = who(matob);
    chanlist = varlist(~contains(varlist,{'SF','Ports','CADD','Channel','CANALOG','new','raw','fileLIST','CMacro','CRAW','CSPK','CLFP'}));
    
    cd(raw)
    load(fileLIST{fi},chanlist{:});
    cd(new)
    save(fileLIST{fi},chanlist{:});
    
    
end


%% Check ttl

% new = 'C:\Users\John\Desktop\sleepQCcheck\samplechck2';
% cd(new)
% mdir = dir('*.mat');
% mdir2 = {mdir.name};
% fileLIST = mdir2;
% 
% ttlC = cell(1,length(fileLIST))
% ttlA = cell(1,length(fileLIST))
% for fi = 1:length(fileLIST)
%     
%     
%     matob = matfile(fileLIST{fi});
%     varlist = who(matob);
%     chanlist = varlist(contains(varlist,{'CDIG'}));
%     
%     if isempty(chanlist)
%        
%         continue
%         
%     else
%         
%         load(fileLIST{fi},'CDIG_IN_1_Up');
%         
%         ttlA{fi} = CDIG_IN_1_Up;
%         
%         dt  = diff(CDIG_IN_1_Up / 44000);
%         
%         
%         ttlC{fi} = dt;
%         
%         clearvars CDIG_IN_1_Up
%     end
%     
%     
%     
%     
%     
%     
%     
% end

%% Single mat file

cd('C:\Users\Admin\Downloads\oldMW2')

raw = 'C:\Users\Admin\Downloads\oldMW2';
new = 'C:\Users\Admin\Downloads\newMW2';

mdir = dir('*.mat');
mdir2 = {mdir.name};
fileLIST = mdir2;


% Load



for fi = 1:length(fileLIST)
    cd(raw)
    matob = matfile(fileLIST{fi});
    varlist = who(matob);
    chanlist = varlist(~contains(varlist,{'SF','Ports','CADD','Channel','CANALOG','new','raw','fileLIST','CMacro','CRAW','CSPK','CLFP'}));
    
    cd(raw)
    load(fileLIST{fi},chanlist{:});
    cd(new)
    save(fileLIST{fi},chanlist{:});
    
    
end