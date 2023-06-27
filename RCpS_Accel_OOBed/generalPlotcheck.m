function [] = generalPlotcheck(recDataLoc)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here




cd(recDataLoc)

matDirall = dir('*.mat');
matDirnames = {matDirall.name};
dataFileName = matDirnames{contains(matDirnames,'Data')};


csvDirall = dir('*.csv');
csvAll = {csvDirall.name};
tempTab = readtable(csvAll{1});



% Load data file
load(dataFileName,'GeneratedData');


genXtemp = GeneratedData.Accel_XSamples;
genYtemp = GeneratedData.Accel_YSamples;
genZtemp = GeneratedData.Accel_ZSamples;

xAxis = 1:length(genXtemp);

plot(xAxis, genXtemp)
hold on; 
plot(xAxis, genZtemp)
plot(xAxis, genYtemp)

xline(tempTab.StartInd,'r-')





end