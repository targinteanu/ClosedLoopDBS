function test_Backend(DQ)

try

%% set filter 
minfac         = 1;    % this many (lo)cutoff-freq cycles in filter
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
load("Saved Data Test\PD24N007001_ARmdl.mat","ARmdl");

%% init 

forecastwin = 1000; % # samples ahead to forecast
buffSize = 20000; % samples
chInd = 1;

%% setup structure 
UserArgs.DAQstatus = true; UserArgs.RunMainLoop = true; 
UserArgs.FilterSetUp = true; UserArgs.MdlSetUp = true;
UserArgs.filtorder = filtorder; UserArgs.BPF = filtwts; 
UserArgs.hicutoff = hico; UserArgs.locutoff = loco; 
UserArgs.Mdl = ARmdl; 
UserArgs.channelIndex = chInd; 
UserArgs.PDSwin1 = forecastwin; UserArgs.PDSwin2 = ceil(.02*forecastwin);
UserArgs.bufferSize = buffSize; UserArgs.bufferSizeGrid = ceil(.02*buffSize);
UserArgs.PhaseOfInterest = [0 pi];

%% loop 
bg_PhaseDetect(UserArgs, DQ, [], ...
    @InitializeRecording_cbmex, @disconnect_cbmex, ...
    @initRawData_cbmex, @getNewRawData_cbmex, []);

catch ME
    getReport(ME)
    send(DQ, ME);
end

end