dT = .01; % s between data requests 
N = round(3/dT); % # of data samples to obtain 
channelIndex = 1;

connect_cbmex(); 
t0 = datetime - seconds(cbmex('time'));
pause(dT);

D = initRawData_cbmex(channelIndex, 60000, t0);
d = D{1}; 
d0 = d; d1 = []; d2 = [];
channelName = d.Properties.VariableNames{1};
figure; ax(1) = subplot(2,1,1);
%plot(d, channelName); 
hold on; grid on;
pause(1);

for ti = 1:N
    try
    D = getNewRawData_cbmex(channelIndex, t0);
    d = D{1}; 
    plot(d, channelName); 
    d1 = d2; d2 = d;
    if ~isempty(d1)
        [d0,d2] = bufferAndRetime(d0, d1, d2);
    end
    catch ME
    end
    disp(['Progress: ',num2str(100*ti/N,3),'%'])
    pause(dT);
end

ax(2) = subplot(2,1,2); 
plot(d0, channelName); grid on; 
linkaxes(ax, 'xy');

disconnect_cbmex();