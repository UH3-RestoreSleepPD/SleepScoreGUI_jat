subplot(4,1,1)
plot(FileMade.BiData01{1})
subplot(4,1,2)
plot(FileMade.BiData23{1})
subplot(4,1,3)
plot(FileMade.corrBD01{1})
subplot(4,1,4)
plot(FileMade.corrBD23{1})

%%
figure;
subplot(4,1,1)
plot(FileMade.TD_key0{1})
subplot(4,1,2)
plot(FileMade.TD_key1{1})
subplot(4,1,3)
plot(FileMade.TD_key2{1})
subplot(4,1,4)
plot(FileMade.TD_key3{1})


%%
figure;
subplot(6,1,1)
plot(FileMade.Klean_left{1,1}.bipolar)
subplot(6,1,2)
plot(FileMade.Klean_left{1,1}.ChannelOne)
subplot(6,1,3)
plot(FileMade.Klean_left{1,1}.ChannelTwo)
subplot(6,1,4)
plot(FileMade.Klean_right{1,1}.bipolar)
subplot(6,1,5)
plot(FileMade.Klean_right{1,1}.ChannelOne)
subplot(6,1,6)
plot(FileMade.Klean_right{1,1}.ChannelTwo)

% isequal(FileMade.Klean_left{1,1}.ChannelTwo,FileMade.Klean_right{1,1}.ChannelOne)
