function [] = fixPENN(problem,datDIR,newDIR)



switch problem

    case 1
        cd(datDIR)
        matdir1 = dir('*.mat');
        matdir2 = {matdir1.name};

        for mi = 1:length(matdir2)

            cd(datDIR)
            tmMat = matdir2{mi};
            tmMatf = matfile(tmMat);
            whoTm = who(tmMatf);

            % Load everything
            load(tmMat)

            e1L = whoTm(contains(whoTm,'E1'));
            e2L = whoTm(contains(whoTm,'E2'));

            ne1L = cell(size(e1L));
            ne2L = cell(size(e1L));

            for ei = 1:length(e2L)

                tmpE1 = e1L{ei};
                %                 load(tmMat,tmpE1);
                tmpE1a = replace(tmpE1,'E1','CEOG_1');
                ne1L{ei} = tmpE1a;
                localfcn(tmpE1a,eval(tmpE1))
                %                 assignin('caller',tmpE1a,eval(tmpE1))
                clearvars(tmpE1)

                tmpE2 = e2L{ei};
                %                 load(tmMat,tmpE2);
                tmpE2a = replace(tmpE2,'E2','CEOG_2');
                ne2L{ei} = tmpE2a;
                localfcn(tmpE2a,eval(tmpE2))
                %                 assignin('caller',tmpE2a,eval(tmpE2))
                clearvars(tmpE2)

            end


            % get new who
            whoTm2a = whoTm(~matches(whoTm,[e1L , e2L]));
            whoTm2b = [whoTm2a ; [ne1L ; ne2L]];
            % Save over file
            cd(newDIR)
            save(tmMat,whoTm2b{:})

        end

    case 2
        cd(datDIR)
        matdir1 = dir('*.mat');
        matdir2 = {matdir1.name};
        for mi = 1:length(matdir2)

            cd(datDIR)
            tmMat = matdir2{mi};
            tmMatf = matfile(tmMat);
            whoTm = who(tmMatf);

            % Load everything
            load(tmMat)

            e1L = whoTm(contains(whoTm,'Chin_1'));
            e2L = whoTm(contains(whoTm,'Chin_2'));
            e3L = whoTm(contains(whoTm,'Chin_z'));

            ne1L = cell(size(e1L));
            ne2L = cell(size(e1L));
            ne3L = cell(size(e1L));

            for ei = 1:length(e2L)

                tmpE1 = e1L{ei};
                tmpE1a = replace(extractAfter(tmpE1,14),'Chin_1','CChin_1');
                ne1L{ei} = tmpE1a;
                localfcn(tmpE1a,eval(tmpE1))
                clearvars(tmpE1)

                tmpE2 = e2L{ei};
                tmpE2a = replace(extractAfter(tmpE2,14),'Chin_2','CChin_2');
                ne2L{ei} = tmpE2a;
                localfcn(tmpE2a,eval(tmpE2))
                clearvars(tmpE2)

                tmpE3 = e3L{ei};
                tmpE3a = replace(extractAfter(tmpE3,14),'Chin_z','CChin_Z');
                ne3L{ei} = tmpE3a;
                localfcn(tmpE3a,eval(tmpE3))
                clearvars(tmpE3)

            end


            % get new who
            whoTm2a = whoTm(~matches(whoTm,[e1L , e2L , e3L]));
            whoTm2b = [whoTm2a ; [ne1L ; ne2L ; ne3L]];
            % Save over file
            cd(newDIR)
            save(tmMat,whoTm2b{:})

        end

end




end





function localfcn(in,out)
assignin('caller',in,out)
end