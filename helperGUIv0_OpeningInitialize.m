function handles = helperGUIv0_OpeningInitialize(handles, ud, svloc)

% start receiver serial communication from paradigm computer
handles.textSrl.String = 'attempting to start serial com here ...';
thisportname = FindMySerialPort();
noSerialSetup = isempty(thisportname);
if ~noSerialSetup
    receiverSerial = serialport(thisportname, 9600);
end
receiverSerial.UserData = ud;
if ~noSerialSetup
configureCallback(receiverSerial,"terminator",...
    @(hsrl,evt)CharSerialCallbackReceiver_Memory_v0(hsrl,evt, ...
                    handles.textSrl, handles.ParadigmInfoTable));
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
handles.stimLastTime = -inf;
handles.stimNewTime = -inf;
handles.stimind = -1;
handles.artReplaceRemaining = [];

% defaults: these should be settable by user input
handles.ArtifactDuration = .04; % set artifact duration (seconds) 
handles.ArtifactStartBefore = .01; % artifact start uncertainty (seconds)
handles.bpthresh = 1000; % min band power cutoff; orig at 1000
handles.StimulatorLagTime = 0.03; % stimulate this many seconds early

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

% data saving
handles.SaveFileLoc = svloc;
handles.SaveFileN = 1;

end