dT = .01; % s between data requests 
N = 100; % # of data samples to obtain 
channelIndex = 1;

connect_cbmex(); 
t0 = datetime - seconds(cbmex('time'));
pause(dT);

D = initRawData_cbmex(channelIndex, 100, t0);
d = D{1}; 
channelName = d.Properties.VariableNames{1};
figure; plot(d, channelName); hold on; grid on;
pause(1);

for ti = 1:N
    try
    D = getNewRawData_cbmex(channelIndex, t0);
    d = D{1}; 
    plot(d, channelName); 
    catch ME
    end
    disp(['Progress: ',num2str(100*ti/N,3),'%'])
    pause(dT);
end

disconnect_cbmex();