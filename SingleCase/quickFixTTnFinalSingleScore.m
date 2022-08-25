%% Quick fix TTn

% SubNUM
subNUM = '2';
institut = 'STAN';
saveLOC = 'F:\SleepStudy2_Aim1';

% Creat final cell
ttnNights = fieldnames(TTn);

finalCS = struct;
for ti = 1:length(ttnNights)

    finalCS.(ttnNights{ti}) = TTn.(ttnNights{ti}).FINAL;

end

%%

saveID = [subNUM, '_' , institut ,'_Fscores.mat'];

%%
cd(saveLOC)
save(saveID,"finalCS")