 
mainDIR = 'F:\SleepStudy2_Aim1\VRAWLFP\LFP';

updir = {'UPEN_1','UPEN_2','UPEN_3','UPEN_4'};

nightS = {'Night 1', 'Night 2', 'Night 3'};

for upD = 1:length(updir)

    tmpD1 = [mainDIR , filesep , updir{upD}];
    cd(tmpD1);

    for ni = 1:length(nightS)
        tmpD2 = [tmpD1 , filesep , nightS{ni}];
        newDIR1 = [tmpD1 , filesep , [nightS{ni} 'b']];
        newDIR2 = [tmpD1 , filesep , [nightS{ni} 'c']];

        if ~exist(newDIR1,'dir')
            mkdir(newDIR1)
        end

        if ~exist(newDIR2,'dir')
            mkdir(newDIR2)
        end

        fixPENN(1,tmpD2,newDIR1)

        if upD == 1
            continue
        else
            fixPENN(2,newDIR1,newDIR2)
        end
    end

end



