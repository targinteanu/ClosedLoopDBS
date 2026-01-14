function handles = helperGUIv0_OpeningInitialize(handles, ud, svloc, srlfunc)
% 
% Automate & standardize the common opening steps in GUIv0 to be run in the
% GUIDE opening function. 
% 
% handles: GUIDE hadles
% ud: default UserData field for serialport object 
% svloc: save location filepath 
% srlfunc: serial receiver callback function that takes in
%          (serialport_handle, event)
% 

% start receiver serial communication from paradigm computer
handles.textSrl.String = 'attempting to start serial com here ...';
thisportname = FindMySerialPort();
thisportconnected = sum(strcmpi(serialportlist, thisportname));
if ~thisportconnected
    resp = questdlg('Proceed without serial connection?', ...
        'Serial not connected.');
    if strcmp(resp, 'Cancel')
        error('Serial not connected; startup aborted.')
    elseif strcmp(resp, 'Yes')
        thisportname = '';
    end
end
noSerialSetup = isempty(thisportname);
if ~noSerialSetup
    receiverSerial = serialport(thisportname, 9600);
end
receiverSerial.UserData = ud;
if ~noSerialSetup
configureCallback(receiverSerial,"terminator",srlfunc);
end
handles.srl = receiverSerial; 
handles.srlLastMsg  = ud.ReceivedData;

% initialize data collection flags
handles.cbmexStatus = false;
handles.StimActive = false;
handles.RunMainLoop = false;
handles.FilterSetUp = false;
handles.MdlSetUp    = false;
handles.showElecGrid = false;

% tracking stimulations and artifact
handles.lastSampleProcTime = -inf;
handles.stimLastTime = -inf;
handles.stimNewTime = -inf;
handles.stimind = -1;
handles.artReplaceRemaining = [];

% defaults: these should be settable by user input
handles.ArtifactDuration = .04; % set artifact duration (seconds) 
handles.ArtifactStartBefore = .01; % artifact start uncertainty (seconds)
handles.bpthresh = 1000; % min band power cutoff; orig at 1000
handles.StimulatorLagTime = 0.03; % stimulate this many seconds early
handles.PhaseOfInterest = [0, pi, nan, nan, nan, nan];
handles.PhaseOfInterestName = ["Peak", "Trough", "Phase3", "Phase4", "Phase5", "Phase6"];
handles.ARlearnrate = 0;

% serial log storage
emptyStorage = nan(100000,1);
handles.pkStorage1 = emptyStorage; handles.pkP1 = 1;
handles.pkStorage2 = emptyStorage; 
handles.trStorage1 = emptyStorage; handles.trP1 = 1;
handles.trStorage2 = emptyStorage; 
handles.stStorage1 = emptyStorage; handles.stP1 = 1;
ud.TimeStamp = nan;
handles.udBlank = ud;
handles.srlStorage1 = repmat(ud,[1000,1]);
handles.srlP1 = 1; 

% empty objects that will be filled in later 
handles.peakDataBuffer = [];
handles.trouDataBuffer = [];
handles.h_peakTrace = [];
handles.h_trouTrace = [];

% data saving
handles.SaveFileLoc = svloc;
handles.SaveFileN = 1;

end