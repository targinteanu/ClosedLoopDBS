function AE_test(UserQueue)

afterEach(UserQueue, @myFunc)

    function myFunc(RecdData)
        DT = string(datetime);
        strOut = "Received data <" + RecdData + "> at time " + DT;
        disp(strOut);
    end

end