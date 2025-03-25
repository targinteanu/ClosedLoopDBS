chsel = 10000:10095;
Tcom = .5*rand(100,1); % computer delay time 
Tars = zeros(size(Tcom)); % ARS delay time 
for it = 1:length(Tcom)
    tic
    pause(.1);
    %[~,t1] = AO_GetLatestTimeStamp();
    [~,~,~,t1] = AO_GetAlignedData(chsel);
    pause(Tcom(it)); 
    pause(.1);
    t3 = toc;
    %[~,t2] = AO_GetLatestTimeStamp();
    [~,~,~,t2] = AO_GetAlignedData(chsel);
    Tcom(it) = t3;
    Tars(it) = t2-t1;
end

Tcom = Tcom(Tars > eps); Tars = Tars(Tars > eps);
figure; plot(Tcom, Tars, '.'); grid on; 
xlabel('computer'); ylabel('\alpha RS');
title(num2str(Tcom\Tars));