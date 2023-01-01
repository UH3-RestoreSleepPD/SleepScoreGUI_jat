function [iniBar , finBar , iniCount , finCount] = getIndivBarData(subNUM , inLIST, iniDATA , finDATA , fintSSu)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

subLISTitems = split(inLIST,{'_','.'});
subIDsall = subLISTitems(:,:,1);
instiIDsall = subLISTitems(:,:,2);
% nightIDsall = subLISTitems(:,:,3);
sortCaseall = cellfun(@(x,y) [x , '_' y], subIDsall , instiIDsall, 'UniformOutput',false);
uniSubList = unique(sortCaseall);
subraw = uniSubList{subNUM};
inSumSsub = iniDATA(matches(sortCaseall,subraw));
finSumSsub = finDATA(matches(sortCaseall,subraw));

iniBar = zeros(3,6);
finBar = zeros(3,6);
iniCount = zeros(3,6);
finCount = zeros(3,6);
% order of sleep states : fintSSu
for si = 1:length(inSumSsub)
    % Get night
    tmpNight_si = inSumSsub{si};
    tmpNight_sf = finSumSsub{si};
    for ei = 1:6
        tNsEi = sum(matches(tmpNight_si,fintSSu{ei})) / length(tmpNight_si);
        tNsEip = round(tNsEi,3);
        iniCount(si,ei) = sum(matches(tmpNight_si,fintSSu{ei}));
        iniBar(si,ei) = tNsEip;

        tNsEf = sum(matches(tmpNight_sf,fintSSu{ei})) / length(tmpNight_sf);
        tNsEfp = round(tNsEf,2);
        finCount(si,ei) = sum(matches(tmpNight_sf,fintSSu{ei}));
        finBar(si,ei) = tNsEfp;
    end
end

end



