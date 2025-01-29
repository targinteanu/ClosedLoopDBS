function [avgLag, stdLag] = stimCalibrate_cerestim(UserArgs, showplot, ...
    connectRecording, disconnectRecording, getTime, initRawData, getNewRawData, ...
    stimSetup, stimShutdown, stimPulse)

if nargin < 2
    showplot = true;
    if nargin < 1
        UserArgs = [];
    end
end

% functions 
if nargin < 10
    if isempty(UserArgs)
        % use BlackRock CereStim 
        stimSetup = @stimSetup_cerestim;
        stimShutdown = @stimShutdown_cerestim;
        stimPulse = @stimPulse_cerestim;
    else
        % use UserArgs 
        stimSetup = UserArgs.HardwareFuncs.SetupStimulator; 
        stimShutdown = UserArgs.HardwareFuncs.ShutdownStimulator; 
        stimPulse = UserArgs.HardwareFuncs.PulseStimulator;
    end
end
if nargin < 7
    if isempty(UserArgs)
        % use BlackRock cbmex 
        connectRecording = @connect_cbmex;
        disconnectRecording = @disconnect_cbmex;
        getTime = @getTime_cbmex;
        initRawData = @initRawData_cbmex;
        getNewRawData = @getNewRawData_cbmex;
    else
        % use UserArgs 
        connectRecording = UserArgs.HardwareFuncs.SetupRecording;
        disconnectRecording = UserArgs.HardwareFuncs.ShutdownRecording;
        getTime = UserArgs.HardwareFuncs.GetTime;
        initRawData = UserArgs.HardwareFuncs.InitRawData;
        getNewRawData = UserArgs.HardwareFuncs.GetNewRawData;
    end
end

if isempty(UserArgs)
    buffSize = 10000;
    [~,~,~,chanInfos,initTic] = initRawData([], buffSize);
    StimSetupArgs = struct(...
        'channel1', 1, ...
        'channel2', 2, ...
        'amp1', 1000, ...
        'amp2', 1000, ...
        'width1', 100, ...
        'width2', 100, ...
        'interphase', 100, ...
        'frequency', 100, ...
        'pulses', 1);
else
    initTic = tic;
    chanInfos = UserArgs.allChannelInfo;
    connectRecording();
    StimSetupArgs = UserData.StimSetupArgs;
end

%% connect to hardware and determine stim trigger channel
% much of this should be replaced into the user interface 
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

time0 = datetime - seconds(getTime(initTic)); % time-of-day when cbmex time was 0

stimulator = stimSetup(StimSetupArgs);

%% perform the test 
disp('Beginning calibration'); tic
maxTestTime = 20; % s; at least 100 will not exceed the buffer 
numTests = 500; StimSentTime = nan(numTests, 1);
maxTestDur = maxTestTime/numTests;
contData = getNewRawData(chanID); % clear out buffer 
for n = 1:numTests
    pause(.01); 
    testDur = rand*maxTestDur; % s 
    testDur = ceil(testDur*1000)/1000; % round up to ms 
    pause(testDur);
    stimtime1 = getTime(initTic);
    stimulator = stimPulse(stimulator);
    stimtime2 = getTime(initTic);
    StimSentTime(n) = mean([stimtime1, stimtime2]);
end
pause(1)
[contData] = getNewRawData(chanID);
toc
disp('Calibration complete.')

%% evaluate the test results 
dataTbl = data2timetable(contData, chanInfo, time0);
dataTbl = dataTbl{1};
StimTrainRec = [false; diff(dataTbl.Variables) > 2000]; % rising edge 
StimTimeRec = dataTbl.Time(StimTrainRec);
StimTimeSent = seconds(StimSentTime) + time0;
StimDelay = StimTimeRec - StimTimeSent;

avgLag = mean(seconds(stimDelay)); stdLag = std(seconds(StimDelay));

if showplot
    figure; subplot(3,1,1); 
    stem(StimTimeSent, 2*ones(size(StimTimeSent)));
    hold on; stem(StimTimeRec, ones(size(StimTimeRec)));
    title('Stimulus Pulses');
    legend('Sent', 'Received', 'Location','eastoutside');
    
    subplot(3,1,2); histogram(seconds(diff(StimTimeSent))); 
    hold on; histogram(seconds(diff(StimTimeRec)));
    title('Inter-Stimulus Time'); xlabel('Seconds');
    legend('Sent', 'Received', 'Location','eastoutside');
    subplot(3,1,3); histogram(seconds(StimDelay));
    title('Received-Sent Delay'); xlabel('Seconds')
end

    disp(['Stim Received - Sent Delay = avg. ', ...
        num2str(avgLag),', SD ',num2str(stdLag),' seconds.'])

%% disconnect from hardware 
disconnectRecording();
stimShutdown(stimulator);

end