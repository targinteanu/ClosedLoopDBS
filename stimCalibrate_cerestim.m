%% connect to hardware and determine stim trigger channel
% much of this should be replaced into the user interface 
buffSize = 10000;
[~,contData,~,chanInfos,initTic] = initRawData_cbmex([], buffSize);
chanInfos = cell2mat(chanInfos);
chanNames = {chanInfos.Name};
chanIDs = [chanInfos.IDnumber];
[chanInd, chanSel] = listdlg("PromptString", "Select Stim Trigger Channel", ...
    "SelectionMode","single", "ListString",chanNames);
if ~chanSel
    warning('Stim Trigger Channel must be selected for calibration.')
    chanInd = find(contains(lower(chanNames), 'ainp1'));
    if isempty(chanInd)
        error('Stim Trigger Channel not found.')
    else
        warning(['Using Channel ',chanNames{chanInd}])
    end
end
chanName = chanNames{chanInd}; chanID = chanIDs(chanInd);
chanInfo = {chanInfos(chanInd)};

time0 = datetime - seconds(getTime_cbmex(initTic)); % time-of-day when cbmex time was 0

stimulator = stimSetup_cerestim(struct(...
    'channel1', 1, ...
    'channel2', 2, ...
    'amp1', 1000, ...
    'amp2', 1000, ...
    'width1', 100, ...
    'width2', 100, ...
    'interphase', 100, ...
    'frequency', 100, ...
    'pulses', 1));

%% perform the test 
disp('Beginning calibration'); tic
maxTestTime = 20; % s; at least 100 will not exceed the buffer 
numTests = 500; StimSentTime = nan(numTests, 1);
maxTestDur = maxTestTime/numTests;
contData = getNewRawData_cbmex(chanID); % clear out buffer 
for n = 1:numTests
    pause(.01); 
    testDur = rand*maxTestDur; % s 
    testDur = ceil(testDur*1000)/1000; % round up to ms 
    pause(testDur);
    stimtime1 = getTime_cbmex(initTic);
    stimulator = stimPulse_cerestim(stimulator);
    stimtime2 = getTime_cbmex(initTic);
    StimSentTime(n) = mean([stimtime1, stimtime2]);
end
pause(1)
[contData] = getNewRawData_cbmex(chanID);
toc
disp('Calibration complete.')

%% evaluate the test results 
dataTbl = data2timetable(contData, chanInfo, time0);
dataTbl = dataTbl{1};
StimTrainRec = [false; diff(dataTbl.Variables) > 2000]; % rising edge 
StimTimeRec = dataTbl.Time(StimTrainRec);
StimTimeSent = seconds(StimSentTime) + time0;

figure; subplot(3,1,1); 
stem(StimTimeSent, 2*ones(size(StimTimeSent)));
hold on; stem(StimTimeRec, ones(size(StimTimeRec)));
title('Stimulus Pulses');
legend('Sent', 'Received', 'Location','eastoutside');

subplot(3,1,2); histogram(seconds(diff(StimTimeSent))); 
hold on; histogram(seconds(diff(StimTimeRec)));
title('Inter-Stimulus Time'); xlabel('Seconds');
legend('Sent', 'Received', 'Location','eastoutside');
StimDelay = StimTimeRec - StimTimeSent;
subplot(3,1,3); histogram(seconds(StimDelay));
title('Received-Sent Delay'); xlabel('Seconds')

%% disconnect from hardware 
disconnect_cbmex();
stimShutdown_cerestim(stimulator);