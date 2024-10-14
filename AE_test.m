function AE_test(UserQueue)

Lst = afterEach(UserQueue, @myFunc);

    function myFunc(RecdData)
        L = evalin('base','Lst');
        DT = string(datetime);
        strOut = "Received data <" + RecdData + "> at time " + DT;
        disp(strOut);
        if strcmpi(RecdData, "start")
            iter = 0;
            while iter < 100
                pause(.5)
                iter = iter+1;
                disp("Process " + num2str(iter) + "% complete...")
            end
        elseif strcmpi(RecdData, "stop")
            delete(L)
        end
    end

end