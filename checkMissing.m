function [] = checkMissing()
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

[fileName,fileLocation] = uigetfile();

cd(fileLocation);
load(fileName,'TT');

tmpN = TT.FINALSCORE;
fracM = (sum(cellfun(@(x) isempty(x), tmpN)) / length(tmpN))*100;
ronFracM = round(fracM,1);

clc
disp(['The total fraction of missing values is '  num2str(ronFracM) ' %']);
disp('..')
disp('..')

if ronFracM < 1
    disp('Done! Please upload!!!');
else
    disp('Please RECHECK still too many missing values!')
end

end