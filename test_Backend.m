function test_Backend(DQ)

try

%% set filter 
minfac         = 2;    % this many (lo)cutoff-freq cycles in filter
min_filtorder  = 15;   % minimum filter length
srate = 1000; 
loco = 13; hico = 30; % Hz 
filtorder = minfac*fix(srate/loco);
TimeShiftFIR = filtorder/(2*srate); % seconds

% build filter 
% usage i.e.: 
% >> filteredSignal = filter(filtwts, 1, unfilteredSignal) 
% -- OR --
% >> filteredSignal = filtfilt(filtwts, 1, unfilteredSignal)
filtwts = fir1(filtorder, [loco, hico]./(srate/2));

%% load AR model 
load(fullfile('Saved Data Test', 'ARmdl.mat'));

%% init 

forecastwin = 1000; % # samples ahead to forecast
buffSize = 20000; % samples
chInd = 65;

%connect_cbmex();
connect_AO();
%[~,~,~,allChannelInfo] = initRawData_cbmex([], ceil(.02*buffSize));
[~,~,~,allChannelInfo] = initRawData_AO([], ceil(.02*buffSize));
%disconnect_cbmex();
disconnect_AO();
UserArgs.allChannelInfo = allChannelInfo;

    selRaw2Art = chInd;
        selRaw2Flt = chInd;
            selFlt2For = 1;
            selFor2Art = 1;
    selRaw2For = [];
    UserArgs.selInds = struct(...
        'selRaw2Flt', selRaw2Flt, ...
        'selFlt2For', selFlt2For, ...
        'selFor2Art', selFor2Art, ...
        'selRaw2Art', selRaw2Art, ...
        'selRaw2For', selRaw2For);

    %{
    [rawD, artD, fltD, forD, timeBuffs, initTic] = ...
        InitializeRecording_cbmex(buffSize, filtorder, forecastwin, ...
        [], selRaw2Art, selRaw2Flt, selRaw2For, selFlt2For);
    %}
    [rawD, artD, fltD, forD, timeBuffs, initTic] = ...
        InitializeRecording_AO(buffSize, filtorder, forecastwin, ...
        [], selRaw2Art, selRaw2Flt, selRaw2For, selFlt2For);
    rawD1 = rawD(1,:); rawD4 = rawD(4,:);
    artD1 = artD(1,:); artD4 = artD(4,:);
    fltD1 = fltD(1,:); fltD4 = fltD(4,:);
    forD1 = forD(1,:); forD4 = forD(4,:);
    timeBuff = timeBuffs{chInd};
    buffSize2 = (buffSize / 1000) * .5 * 50;
    buffSize2 = ceil(buffSize2);
    forBuff = nan(buffSize2,2);
    tSt = nan(buffSize2,1);
    recDataStructs.forBuffs = {forBuff}; recDataStructs.stimBuff = tSt;
    for v = ["rawD", "artD", "fltD", "forD", "timeBuffs", "initTic"]
        eval("recDataStructs."+v+" = "+v+";");
    end
    UserArgs.recDataStructs = recDataStructs;

%% setup structure 
UserArgs.StimSetupArgs = [];
UserArgs.StimTriggerMode = true;
UserArgs.DAQstatus = true; UserArgs.RunMainLoop = true; 
UserArgs.FilterSetUp = true; UserArgs.MdlSetUp = true;
UserArgs.check_artifact_Value = true;
UserArgs.ControllerResult = true;
UserArgs.StimulatorLagTime = .01; 
UserArgs.FilterOrder = filtorder; UserArgs.BPF = filtwts; 
UserArgs.hicutoff = hico; UserArgs.locutoff = loco; 
UserArgs.Mdl = ARmdl; 
UserArgs.channelIndex = chInd; 
UserArgs.allChannelIDs = [];
UserArgs.PDSwin1 = forecastwin; UserArgs.PDSwin2 = ceil(.02*forecastwin);
UserArgs.bufferSize = buffSize; UserArgs.bufferSizeGrid = ceil(.02*buffSize);
UserArgs.PhaseOfInterest = [0 pi];
UserArgs.StimActive = true; UserArgs.stimMaxFreq = 50;
UserArgs.check_artifact.Value = true;
UserArgs.SerialArgs = struct('UserData',struct('ReceivedData',''), 'NoSerial',true);

%% loop 
%{
bg_PhaseDetect(UserArgs, DQ, [], ...
    @connect_cbmex, @disconnect_cbmex, ...
    @stimSetup_cerestim, @stimShutdown_cerestim, ...
    ...@stimPulse_cerestim, @(~) [], ...
    @stimPulse_cpod, @srlSetup_cpod, ...
    @getNewRawData_cbmex, @getTime_cbmex);
%}
%bg_PhaseDetect_BlackRock(UserArgs, DQ, @Controller_PDS_PD);
bg_PhaseDetect(UserArgs, DQ, [], ...
    @connect_AO, @disconnect_AO, ...
    @(~) 0, @(~,~) 0, ...
    @(~,~) 0, @(~) [], ...
    @getNewRawData_AO, @getTime_AO);

catch ME
    getReport(ME)
    send(DQ, ME);
end

end