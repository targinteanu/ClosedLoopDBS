%function AE_test(UserQueue)

Lst = afterEach(UserQueue, @myFunc);

    function myFunc(RecdData)
        L = evalin('base','Lst');
        DT = string(datetime);
        strOut = "Received data <" + RecdData + "> at time " + DT;
        disp(strOut);
        if strcmpi(RecdData, "start")
            Niter = 20;
            iter = 0;
            while iter < Niter
                pause(.5)
                iter = iter+1;
                disp("Process " + num2str(100*iter/Niter) + "% complete...")
            end
        elseif strcmpi(RecdData, "stop")
            delete(L)
        end
    end

%end