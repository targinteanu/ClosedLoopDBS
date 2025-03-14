Tcom = .1*rand(100,1); % computer delay time 
Tars = zeros(size(Tcom)); % ARS delay time 
for it = 1:length(Tcom)
    tic
    [~,t1] = AO_GetLatestTimeStamp();
    pause(Tcom(it)); 
    t3 = toc;
    [~,t2] = AO_GetLatestTimeStamp();
    Tcom(it) = t3;
    Tars(it) = t2-t1;
end

figure; plot(Tcom, Tars, '.'); grid on; 
xlabel('computer'); ylabel('\alpha RS');