function varargout = GUI_Memory_v1(varargin)
% GUI_MEMORY_V1 MATLAB code for GUI_Memory_v1.fig
%      GUI_MEMORY_V1, by itself, creates a new GUI_MEMORY_V1 or raises the existing
%      singleton*.
%
%      H = GUI_MEMORY_V1 returns the handle to a new GUI_MEMORY_V1 or the handle to
%      the existing singleton*.
%
%      GUI_MEMORY_V1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_MEMORY_V1.M with the given input arguments.
%
%      GUI_MEMORY_V1('Property','Value',...) creates a new GUI_MEMORY_V1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_Memory_v1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_Memory_v1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_Memory_v1

% Last Modified by GUIDE v2.5 28-Mar-2025 17:45:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_Memory_v1_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_Memory_v1_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before GUI_Memory_v1 is made visible.
function GUI_Memory_v1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_Memory_v1 (see VARARGIN)

% Choose default command line output for GUI_Memory_v1
handles.output = hObject;

% --- Begin My Code ---
wb = waitbar(0, 'Starting up GUI...');

% change ax_polar to polar axes 
waitbar(0, wb, 'Setting up polar plot...')
axes(handles.ax_polar);
polarhistogram(pi/2);
handles.ax_polar = gca;

% Create a timer object that will be used to grab data and refresh
% analysis/plotting
handles.timer = timer(...
    'ExecutionMode', 'fixedSpacing', ...       % Run timer repeatedly
    'Period', 0.01, ...                      % Initial period is 100 ms
    'TimerFcn', {@updateDisplay,hObject}, ... % callback function.  Pass the figure handle
    'StartFcn', {@StartMainLoop,hObject}, ...
    'StopFcn',  {@StopMainLoop,hObject}, ...     % callback to execute when timer starts
    'ErrorFcn', {@TimerError,hObject}, ...
    'BusyMode', 'error'); 

% start receiver serial communication from paradigm computer
waitbar(.05, wb, 'Setting up serial com...')
handles.textSrl.String = 'attempting to start serial com here ...';
thisportname = FindMySerialPort();
noSerialSetup = isempty(thisportname);
ud = struct('ReceivedData', '', ...
            'TrialNumber', -1, ...
            'StimOn', false, ...
            'ParadigmPhase', 'WAIT', ...
            'ImageVisible', false);
handles.SerialArgs = struct('UserData', ud, ...
                         'CallbackFcn', @CharSerialCallbackReceiver_Memory_v0, ...
                         'PortName', thisportname, ...
                         'NoSerial', noSerialSetup);
handles.srlHere = false;
try
    handles = connectSerial(handles);
catch ME1
    getReport(ME1)
    % try delete(instrfind) ??
    keyboard
end

% serial saving 
ud.TimeStamp = nan;
handles.udBlank = ud;
handles.srlStorage1 = repmat(ud,[1000,1]);
handles.srlP1 = 1; 

% initiate other vars ...
waitbar(.1, wb, 'Initializing variables...')
handles.DAQstatus = false;
handles.StimActive = false;
handles.RunMainLoop = false;
handles.FilterSetUp = false;
handles.MdlSetUp    = false;
handles.showElecGrid = false;
handles.srlLastMsg  = ud.ReceivedData;
handles.bufferSize     = 10; 
handles.bufferSizeGrid = 10;
handles.allChannelIDs = [];
handles.channelIndex = [];
handles.PhaseOfInterest = [0, pi];

% init other storage buffers
handles.phStorage = nan(100000,width(handles.PhaseOfInterest)); 
handles.phP = 1;
handles.stStorage = nan(100000,1); handles.stP = 1;

% init empty plot handles to avoid errors later 
handles.h_rawDataTrace = [];
handles.h_artDataTrace = [];
handles.h_filtDataTrace = []; 
handles.h_timingTrace = [];
handles.h_timeDispTrace = [];
handles.h_predTrace = [];
handles.h_peakTrace = [];
handles.h_trouTrace = [];
handles.h_stimTrace = [];
handles.h_peakPhase = [];
handles.h_trouPhase = [];
handles.elecGridImg = [];

% file saving 
waitbar(.15, wb, 'Setting file save location...')
svloc = ['Saved Data Memory',filesep,'Saved Data ',...
    datestr(datetime, 'yyyy-mm-dd HH.MM.SS')];
pause(1)
mkdir(svloc); 
handles.SaveFileLoc = svloc;
handles.SaveFileName = [svloc,filesep,'SaveFile'];
handles.SaveFileN = 1;

% default values 
waitbar(.2, wb, 'Setting default values...')
%handles.channelIndex = get(handles.pop_channels,'Value'); 
PDSwin = str2double(get(handles.txt_PDSwin,'String'));
PDSwin = ceil(PDSwin*1000); handles.PDSwin1 = PDSwin;
handles.PDSwin2 = ceil(.02*PDSwin); 
handles.bufferSize = str2double(get(handles.txt_display,'String')) * 1000;
handles.bufferSizeGrid = str2double(get(handles.txt_griddur,'String')) * 1000;
handles.stimMaxFreq = eval(get(handles.txt_MaxStimFreq, 'String'));
handles.StimulatorLagTime = 0.012; % or should it start at 0?
handles.check_artifact_Value = false;
handles.ControllerResult = 0;

% hardware-specific functions 
waitbar(.25, wb, 'Setting up hardware...')
[handles.HardwareFuncs, handles.StimTriggerMode] = helperGUIv1_DefHardwareFuncs();
handles.initTic = tic;

% start parallel pool(s) 
waitbar(.3, wb, 'Starting parallel pool...')
handles.pool = gcp('nocreate');
if isempty(handles.pool)
    handles.pool = parpool;
end

% Set up a data queue(s)
waitbar(.95, wb, 'Setting parallel data queue(s)...')
%handles.userQueue = parallel.pool.PollableDataQueue;
handles.dataQueue = parallel.pool.PollableDataQueue;
handles.stimQueue = parallel.pool.PollableDataQueue;
%handles.modlQueue = parallel.pool.PollableDataQueue;
handles.f_PhaseDetect = [];

% remove these fields when sending handles over the userQueue
handles.rmfieldList = {...
    'dataQueue', 'stimQueue', 'pool', ...
    'f_PhaseDetect', ...
    'udBlank', 'srlStorage1', 'srlP1', 'phStorage', 'phP', 'stStorage', 'stP', ...
    'output', ...
    'timer', 'srl', ...
    'HardwareFuncs', ...
    'h_rawDataTrace', 'h_artDataTrace', 'h_filtDataTrace', 'h_timingTrace', 'h_timeDispTrace', ...
    'h_predTrace', 'h_peakTrace', 'h_trouTrace', 'h_stimTrace', ...
    'h_peakPhase', 'h_trouPhase', ...
    'figure1', 'scribeOverlay', 'output', ...
    'pnl_elecgrid', 'pnl_stim', 'pnl_filt', 'pnl_controls', ...
    'check_polar', 'check_artifact', ...
    'txt_gridmax', 'txt_gridmin', 'txt_griddur', ...
    'txt_MaxStimFreq', 'txt_interphase', 'txt_width2', 'txt_width1', 'txt_amp2', 'txt_amp1', ...
    'txt_AR', 'txt_hico', 'txt_loco', 'txt_PDSwin', 'txt_display', ...
    'lbl_channel', 'lbl_display', ...
    'textSrl', ...
    'text25', 'text26', 'text24', 'text23', 'text22', 'text21', 'text20', ...
    'text19', 'text17', 'text16', 'text15', ...
    'text14', 'text12', 'text13', 'text11', 'text9', 'text8', 'text7', 'text6', ...
    'text27', 'text28', 'text29', ...
    'ax_elecgrid', 'ax_polar', 'ax_timing', 'ax_filt', 'ax_raw', ...
    'elecGridImg', ...
    'pop_elecgrid', 'pop_HoldStim', 'pop_EncodeStim', 'pop_DecodeStim', ...
    'pop_channel5', 'pop_channel4', 'pop_channel3', 'pop_channel2', 'pop_channel1', 'pop_channels', ...
    'tgl_stim', 'tgl_StartStop', ...
    'push_AR', 'push_filter', 'push_remchan', 'push_stimCalibrate', ...
    'cmd_cbmexOpen', 'cmd_cbmexClose'};
try
    h1 = rmfield(handles, handles.rmfieldList);
catch ME2
    getReport(ME2)
    % fix issue with removed field list here
    keyboard
end

% Update handles structure
waitbar(1, wb, 'Updating GUI handles...')
pause(.01)
close(wb)
pause(.01)
guidata(hObject, handles);

clc
%clear all

% UIWAIT makes GUI_Memory_v1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_Memory_v1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% ----------------------------------------------------------------------- %
% ----                                                                --- %
% ----        Figure Objects Create and Callback Functions            --- %
% ----                                                                --- %
% ----------------------------------------------------------------------- %

% --- Executes on button press in cmd_cbmexOpen.
function cmd_cbmexOpen_Callback(hObject, eventdata, handles)
% hObject    handle to cmd_cbmexOpen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try

handles = helperGUIv1_DAQopen(handles);

guidata(hObject,handles)

catch ME
    getReport(ME)
    errordlg(ME.message, 'Connect DAQ issue')
end

% --- Executes on button press in cmd_cbmexClose.
function cmd_cbmexClose_Callback(hObject, eventdata, handles)
% hObject    handle to cmd_cbmexClose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.DAQstatus = false;
[dataRecd, handles.SaveFileN, ~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~, ...
    handles.phStorage, handles.phP, handles.stStorage, handles.stP] = ...
    pollDataQueue_PhaseDetect_v1(handles.dataQueue, handles.channelIndex, ...
        handles.SaveFileName, handles.SaveFileN, handles.time0, 10, ...
        handles.phStorage, handles.phP, handles.stStorage, handles.stP);
if ~dataRecd
    warning('Polling data queue timed out.')
    keyboard
end
cancel(handles.f_PhaseDetect); 
cancelAll(handles.pool.FevalQueue);

% handles = connectSerial(handles); % ?
guidata(hObject,handles)

function txt_display_Callback(hObject, eventdata, handles)
% hObject    handle to txt_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_display as text
%        str2double(get(hObject,'String')) returns contents of txt_display as a double
settingChange(hObject)

% --- Executes during object creation, after setting all properties.
function txt_display_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function tgl_StartStop_ButtonDownFcn(hObject, eventdata, handles)
keyboard
tgl_StartStop_Callback(hObject, eventdata, handles)
% Unclear why the code sometimes gets here instead of Callback below. 


% --- Executes on button press in tgl_StartStop.
function tgl_StartStop_Callback(hObject, eventdata, handles)
% hObject    handle to tgl_StartStop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tgl_StartStop

handles = guidata(hObject);

% if Start
if get(hObject,'Value') == 1
    
    % Check to make sure cbmex connection is open
    if ~handles.DAQstatus
        errordlg('No DAQ connection.  Open connection before starting','Not Connected')
        return
    end
    
    set(hObject,'String','Stop');

    % This starts the timer and also executes the StartFnc which grabs the
    % data, creates the buffer and plots the first bit of data
    start(handles.timer)
    
% Stop
else
    set(hObject,'String','Start')
    stop(handles.timer)
end

% --- Executes on selection change in pop_channels.
function pop_channels_Callback(hObject, eventdata, handles)
% hObject    handle to pop_channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_channels contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_channels
%hObject.Value
%hObject.String(hObject.Value,:)
settingChange(hObject)

% --- Executes during object creation, after setting all properties.
function pop_channels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

drawnow

try
% stop stim
if handles.StimActive
    handles.StimActive = false;
    guidata(hObject, handles)
    setTgl(hObject, eventdata, handles, handles.tgl_StartStop, 0);
end
catch ME1
    getReport(ME1)
end

try
handles.HardwareFuncs.ShutdownRecording();
handles.DAQstatus = false;
[dataRecd, handles.SaveFileN, ~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~, ...
    handles.phStorage, handles.phP, handles.stStorage, handles.stP] = ...
    pollDataQueue_PhaseDetect_v1(handles.dataQueue, handles.channelIndex, ...
        handles.SaveFileName, handles.SaveFileN, handles.time0, 10, ...
        handles.phStorage, handles.phP, handles.stStorage, handles.stP);
if ~dataRecd
    warning('Polling data queue timed out.')
    %keyboard
end

if ~isempty(handles.f_PhaseDetect)
    cancel(handles.f_PhaseDetect); 
end
cancelAll(handles.pool.FevalQueue);
guidata(hObject, handles)
stop(handles.timer)
delete(handles.timer)
if ~handles.SerialArgs.NoSerial
    delete(handles.srl);
end
catch ME3
    getReport(ME3)
end

try
% save stored data 
% TO DO: if Stim/PeakTrough are all nan or SerialLog has only nan
% timestamps, do not save that variable 
Stim = handles.stStorage; PeakTrough = handles.phStorage;
SerialLog = handles.srlStorage1;
svfn = [handles.SaveFileName,num2str(handles.SaveFileN),'.mat'];
disp(['Saving Remaining Data to ',svfn])
save(svfn, 'SerialLog', 'PeakTrough', 'Stim');
disp('Saving all data...')
ConsolidateSavedData(handles.SaveFileLoc)
catch ME2
    getReport(ME2)
end

% Hint: delete(hObject) closes the figure
delete(hObject);


% ----------------------------------------------------------------------- %
% ----                                                                --- %
% ----          Main Loop, Runs every time the timer fires            --- %
% ----                                                                --- %
% ----------------------------------------------------------------------- %
function updateDisplay(obj, evt, hObject)

handles = guidata(hObject);

try
    
    % ensure connection with hardware 
    handles = guidata(hObject);
    if ~handles.DAQstatus
        stop(handles.timer)
    end

    % timing 
    % pause(.05); % s between displays - is this necessary since timer has
    %                                    fixed spacing period?
    tNow = datetime;
    timeDisp2 = handles.timeDisp0 + toc(handles.timeDisp1);
    handles.timeDispBuff = bufferData(handles.timeDispBuff, timeDisp2);
    guidata(hObject, handles);

    % get data from Central
    if handles.dataQueue.QueueLength > 1000 
        % LIMIT DATA QUEUE LENGTH
        cancel(handles.f_PhaseDetect);
        cancelAll(handles.pool.FevalQueue);
        %{
        handles = connectSerial(handles);
        guidata(hObject, handles);
        %}
        error('Data Queue Overflow');
    end
    [dataRecd, handles.SaveFileN, timeBuff, forBuff, tSt, ...
    tPltRng, rawPlt, fltPlt, forPlt, artPlt, ...
    rawD1, rawD4, fltD1, fltD4, forD1, forD4, artD1, artD4, ...
    handles.phStorage, handles.phP, handles.stStorage, handles.stP] = ...
    pollDataQueue_PhaseDetect_v1(handles.dataQueue, handles.channelIndex, ...
        handles.SaveFileName, handles.SaveFileN, handles.time0, 10, ...
        handles.phStorage, handles.phP, handles.stStorage, handles.stP);
    if ~dataRecd
        error('Data aquisition timed out.')
    end
    lastSampleProcTime = timeBuff(end);
    rawIDs = cellfun(@(s) s.IDnumber, rawD1);

    % x axes alignment 
    common_xlim = tPltRng - tNow;
    common_xdiff = diff(common_xlim); 
    ext_xdiff = common_xdiff * handles.ax_filt.InnerPosition(3) / ...
        handles.ax_raw.InnerPosition(3); 
    ext_xlim = [0, ext_xdiff] + common_xlim(1); % align left

    % update serial log 
    % TO DO: there should be a better way to do this; serial callback
    % should trigger an event or listener that logs the info 
    if handles.srlHere
        if handles.FilterSetUp
            if numel(fltPlt) % should this be necessary when above is met?
                ControllerLastResult = handles.ControllerResult;
                ControllerResult = Controller_PDS_Memory( handles.srl, ...
                    rmfield(handles, handles.rmfieldList), ...
                    fltPlt{(end-handles.PDSwin1+1):end, : } );
                if ControllerResult ~= ControllerLastResult
                    handles.ControllerResult = ControllerResult;
                    guidata(hObject, handles);
                    requeryPhaseDetect(hObject, 1);
                end
            end
        end
    ReceivedData = handles.srl.UserData.ReceivedData; 
    if ~strcmp(ReceivedData, handles.srlLastMsg)
        ud = handles.srl.UserData; 
        ud.TimeStamp = lastSampleProcTime;
        if handles.srlP1 <= length(handles.srlStorage1)
            handles.srlStorage1(handles.srlP1) = ud;
            handles.srlP1 = handles.srlP1+1;
        else
            ud = handles.udBlank;
            % storage full; save
            SerialLog = handles.srlStorage1;
            svfn = [handles.SaveFileLoc,filesep,'SaveFile',num2str(handles.SaveFileN),'.mat'];
            disp(['Saving Serial to ',svfn])
            save(svfn,'SerialLog');
            handles.SaveFileN = handles.SaveFileN + 1;
            handles.srlP1 = 1;
            handles.srlStorage1 = repmat(ud, size(handles.srlStorage1));
        end
    end
    handles.srlLastMsg = ReceivedData;
    end

    % update electrode grid
    if handles.showElecGrid
        try
            handles.elecGridImg.CData = helperGUIv1_ElectrodeGridUpdate(...
                handles.elecGridImg, handles.elecGridFunc, ...
                handles.channelIDlist, rawIDs, rawD4, ...
                handles.bufferSizeGrid, handles.fSamples);
        catch ME4
            getReport(ME4)
            errordlg(ME4.message, 'Electrode Grid Issue');
            handles.showElecGrid = false;
            pause(.01);
        end
    end

    if handles.FilterSetUp
        try
        % update filtered data plot
        handles = helperGUIv1_plotFlt(handles, tNow, fltPlt, ext_xlim);

        if handles.MdlSetUp
            try
                handles = helperGUIv1_plotMdl(handles, tNow, fltPlt, forPlt, forBuff, tSt, artPlt);
            catch ME2
                getReport(ME2)
                errordlg(ME2.message, 'Model Prediction Issue');
                handles.MdlSetUp = false;
                guidata(hObject, handles);
                pause(.01);
            end
        end

        catch ME1
            getReport(ME1)
            errordlg(ME1.message, 'Filtering Issue');
            handles.FilterSetUp = false;
            handles.MdlSetUp = false;
            guidata(hObject, handles);
            pause(.01);
        end
    end

    % update raw data and timing plots
    handles = helperGUIv1_plotRaw(handles, tNow, rawPlt, timeBuff, common_xlim);

    guidata(hObject,handles)

catch ME 
    getReport(ME)
    guidata(hObject, handles)
    stop(handles.timer)
    keyboard
end
    
% ----------------------------------------------------------------------- %
% ----                                                                --- %
% ----                      Timer Start Function                      --- %
% ----                                                                --- %
% ----------------------------------------------------------------------- %

% Runs once when timer starts.  It initializes plot and buffer and
% accommodates any new selection by user.
function  StartMainLoop(obj, evt, hObject)
% Put the whole function in a try-catch block.  This makes debugging much
% easier because it captures the error and displays a report to the Matlab
% Command Window

handles = guidata(hObject);

try

    tNow = datetime;

    % setup data and initiate raw and timing plots
    [handles, fltPlt, forPlt, forBuff, tSt, common_xlim, unitname] = ...
        helperGUIv1_plotSetupRaw(handles, tNow);

    % initiate filtered data plot
    if handles.FilterSetUp
        try
            handles = helperGUIv1_plotSetupFltMdl(handles, tNow, ...
                fltPlt, forPlt, forBuff, tSt, common_xlim, unitname);
        catch ME1
            getReport(ME1)
            errordlg(ME1.message, 'Filtering Issue');
            handles.FilterSetUp = false;
            handles.MdlSetUp = false;
            guidata(hObject, handles);
            requeryPhaseDetect(hObject, 1);
            handles = guidata(hObject);
            pause(.01);
        end
    end

    % handles = disconnectSerial(handles);
    handles.RunMainLoop = true; 
    guidata(hObject,handles)
    requeryPhaseDetect(hObject, 1);
    
catch ME
    getReport(ME)
    handles.RunMainLoop = false; 
    guidata(hObject, handles);
    requeryPhaseDetect(hObject, 1);
    stop(handles.timer);
    %{
    handles = connectSerial(handles);
    guidata(hObject, handles);
    %}
    % keyboard
end

% ----------------------------------------------------------------------- %
% ----                                                                --- %
% ----                      Timer Stop Function                       --- %
% ----                                                                --- %
% ----------------------------------------------------------------------- %

% Runs once when timer stops.
function  StopMainLoop(obj, evt, hObject)
% Put the whole function in a try-catch block.  This makes debugging much
% easier because it captures the error and displays a report to the Matlab
% Command Window

handles = guidata(hObject);

try
    handles.RunMainLoop = false;
    guidata(hObject, handles)
    requeryPhaseDetect(hObject, 1);
    % handles = connectSerial(handles);
    guidata(hObject, handles);
catch ME
    getReport(ME)
    keyboard
end


function TimerError(obj, evt, hObject)
% timer has encountered an error! Why was it not caught?
handles = guidata(hObject);
keyboard

% ----------------------------------------------------------------------- %
% ----                                                                --- %
% ----                      Helper Functions                          --- %
% ----                                                                --- %
% ----------------------------------------------------------------------- %

function settingChange(hObject)
handles = guidata(hObject);

% if the timer is running, stop it and restart it (which will use the newly
% selected channel.  If the timer isn't running, don't do anything.
if strcmp(handles.timer.Running,'on')
    stop(handles.timer)
    start(handles.timer)
end


% --- New Helpers ---


function requeryPhaseDetect(hObject, timeoutdur)
handles = guidata(hObject);
if timeoutdur >= 0
[dataRecd, handles.SaveFileN, ~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~, ...
    handles.phStorage, handles.phP, handles.stStorage, handles.stP] = ...
    pollDataQueue_PhaseDetect_v1(handles.dataQueue, handles.channelIndex, ...
        handles.SaveFileName, handles.SaveFileN, handles.time0, timeoutdur, ...
        handles.phStorage, handles.phP, handles.stStorage, handles.stP);
if ~dataRecd
    warning('Polling data queue timed out.')
    %hObject
    %eventdata
    %keyboard
end
end
try
if ~isempty(handles.f_PhaseDetect)
cancel(handles.f_PhaseDetect); 
% For some reason, after push_filter, the above does not immediately cancel
% all RunningFutures, but waiting at least 3 seconds will make them empty
if handles.FilterSetUp
    pause(6)
end
% cancellAll may be necessary in case anything is still running, but if the
% FevalQueue is not empty, it causes annoying problems like extra GUI
% windows trying to open or opening. 
cancelAll(handles.pool.FevalQueue);
end
% handles = disconnectSerial(handles);
handles.f_PhaseDetect = parfeval(handles.pool, @bg_PhaseDetect, 1, ...
    rmfield(handles, handles.rmfieldList), ...
    handles.dataQueue, handles.stimQueue, ...
    handles.HardwareFuncs.SetupRecording, handles.HardwareFuncs.ShutdownRecording, ...
    handles.HardwareFuncs.SetupStimulator, handles.HardwareFuncs.ShutdownStimulator, handles.HardwareFuncs.PulseStimulator, handles.HardwareFuncs.SetupStimTTL, ...
    handles.HardwareFuncs.GetNewRawData, handles.HardwareFuncs.GetTime);
catch ME2
    warning(ME2.message);
    % if there is a problem here, consider stopping everything 
    keyboard
end
guidata(hObject, handles);


function handles = connectSerial(handles)
if handles.srlHere
    % should already be connected? 
    keyboard
else
try
ud = handles.SerialArgs.UserData; 
CBFn = handles.SerialArgs.CallbackFcn;
thisportname = handles.SerialArgs.PortName; 
noSerialSetup = handles.SerialArgs.NoSerial;
if ~noSerialSetup
    try
        pause(.1);
        receiverSerial = serialport(thisportname, 9600);
    catch ME1 
        if contains(ME1.message, ...
                'Verify that a device is connected to the port, the port is not in use')
            % the port may need some time (how much?) before it is cleared 
            warning('Failed to connect serial port on first try; trying again soon...')
            pause(3);
            receiverSerial = serialport(thisportname, 9600);
        else
            rethrow(ME1)
        end
    end
end
receiverSerial.UserData = ud;
if ~noSerialSetup
configureCallback(receiverSerial,"terminator",...
    @(hsrl,evt)CBFn(hsrl,evt, ...
                    handles.textSrl, handles.ParadigmInfoTable)); 
end
handles.textSrl.String = 'Serial is connected on user thread.';
handles.srlHere = true;
handles.srl = receiverSerial;
catch ME
    getReport(ME)
    handles.textSrl.String = ME.message;
end
end

function handles = disconnectSerial(handles)
if ~handles.SerialArgs.NoSerial
    delete(handles.srl);
end
handles.srlHere = false;
handles.textSrl.String = 'Serial disconnected from user thread.'; 


function setTgl(hObject, eventdata, handles, hTgl, newValue)
% set a toggle button to a desired Value and activate its callback if it is
% not currently at that value. 
curValue = hTgl.Value; 
hTgl.Value = newValue; guidata(hObject, handles);
if ~(curValue == newValue)
    hTgl.Callback(hTgl, eventdata);
    guidata(hObject, handles);
end
 

% ----------------------------------------------------------------------- %
% ----                                                                --- %
% ----                        New Callbacks                           --- %
% ----                                                                --- %
% ----------------------------------------------------------------------- %


% --- Executes on selection change in pop_EncodeStim.
function pop_EncodeStim_Callback(hObject, eventdata, handles)
% hObject    handle to pop_EncodeStim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_EncodeStim contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_EncodeStim
contents = cellstr(get(hObject,'String'));
EncodeStim = contents{get(hObject,'Value')};
handles.StimMode.ENCODE = EncodeStim; 
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function pop_EncodeStim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_EncodeStim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
pop_EncodeStim_Callback(hObject, eventdata, handles)


% --- Executes on selection change in pop_DecodeStim.
function pop_DecodeStim_Callback(hObject, eventdata, handles)
% hObject    handle to pop_DecodeStim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_DecodeStim contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_DecodeStim
contents = cellstr(get(hObject,'String'));
DecodeStim = contents{get(hObject,'Value')};
handles.StimMode.DECODE = DecodeStim; 
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function pop_DecodeStim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_DecodeStim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
pop_DecodeStim_Callback(hObject, eventdata, handles)


% --- Executes on selection change in pop_HoldStim.
function pop_HoldStim_Callback(hObject, eventdata, handles)
% hObject    handle to pop_HoldStim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_HoldStim contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_HoldStim
contents = cellstr(get(hObject,'String'));
HoldStim = contents{get(hObject,'Value')};
handles.StimMode.HOLD = HoldStim; 
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function pop_HoldStim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_HoldStim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
pop_HoldStim_Callback(hObject, eventdata, handles)



function txt_loco_Callback(hObject, eventdata, handles)
% hObject    handle to txt_loco (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_loco as text
%        str2double(get(hObject,'String')) returns contents of txt_loco as a double


% --- Executes during object creation, after setting all properties.
function txt_loco_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_loco (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_hico_Callback(hObject, eventdata, handles)
% hObject    handle to txt_hico (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_hico as text
%        str2double(get(hObject,'String')) returns contents of txt_hico as a double


% --- Executes during object creation, after setting all properties.
function txt_hico_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_hico (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_filter.
% setup the filter, as in Myeegfilt
function push_filter_Callback(hObject, eventdata, handles)
% hObject    handle to push_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

stop(handles.timer); % to avoid dataQueue overflow 

% details from data 
srate = handles.fSample;

% filtering bound rules 
minfac         = 2;    % this many (lo)cutoff-freq cycles in filter
min_filtorder  = 15;   % minimum filter length

% access desired parameters
hicutoff = str2double(handles.txt_hico.String); 
locutoff = str2double(handles.txt_loco.String); 

% build filter
[filtwts, filtorder] = buildFIRBPF(srate, locutoff, hicutoff, minfac, min_filtorder);

handles.hicutoff = hicutoff; handles.locutoff = locutoff; 
handles.FilterOrder = filtorder;
handles.TimeShiftFIR = filtorder/(2*srate); % seconds 
handles.IndShiftFIR = ceil(filtorder/2); % samples ???

handles.BPF = filtwts; 

% restart timer and plots
pause(.01);
handles.FilterSetUp = true;
guidata(hObject, handles)
pause(.01);
start(handles.timer);



function txt_PDSwin_Callback(hObject, eventdata, handles)
% hObject    handle to txt_PDSwin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_PDSwin as text
%        str2double(get(hObject,'String')) returns contents of txt_PDSwin as a double


% --- Executes during object creation, after setting all properties.
function txt_PDSwin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_PDSwin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_AR_Callback(hObject, eventdata, handles)
% hObject    handle to txt_AR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_AR as text
%        str2double(get(hObject,'String')) returns contents of txt_AR as a double


% --- Executes during object creation, after setting all properties.
function txt_AR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_AR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_AR.
% estimate the A.R. model
function push_AR_Callback(hObject, eventdata, handles)
% hObject    handle to push_AR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

stop(handles.timer); % to avoid dataQueue overflow 

n = str2double(get(handles.txt_AR,'String'));
N = str2double(get(handles.txt_PDSwin,'String'));
PDSwin = ceil(N*handles.fSample); handles.PDSwin1 = PDSwin;
handles.PDSwin2 = ceil(.02*PDSwin); 

try

    % catch mistakes 
    if PDSwin > handles.bufferSize
        error('Phase estimation window cannot be larger than display window.')
    end
    if handles.PDSwin1 <= handles.IndShiftFIR
        error('Forecast length does not overcome filter delay.')
    end

    [dataRecd, handles.SaveFileN, timeBuff, forBuff, stimBuff, ...
    tPltRng, rawPlt, fltPlt, forPlt, artPlt, ...
    rawD1, rawD4, fltD1, fltD4, forD1, forD4, artD1, artD4, ...
    handles.phStorage, handles.phP, handles.stStorage, handles.stP] = ...
    pollDataQueue_PhaseDetect_v1(handles.dataQueue, handles.channelIndex, ...
        handles.SaveFileName, handles.SaveFileN, handles.time0, 10, ...
        handles.phStorage, handles.phP, handles.stStorage, handles.stP);
    if ~dataRecd
        error('Data aquisition timed out.')
    end

    y = fltD4{1}(:,2); 
    L = min(length(y), 3*PDSwin) - 1;
    y = y((end-L):end);
    y = iddata(y,[],1/handles.fSample);
    ARmdl = ar(y,n,'yw');
    
    handles.Mdl = ARmdl; 
    handles.MdlSetUp = true;
    
% restart timer and plots
pause(.01);
guidata(hObject, handles)
pause(.01);
start(handles.timer);

catch ME
    getReport(ME)
    errordlg(ME.message, 'Model ID Issue')
end


% --- Executes on button press in check_polar.
function check_polar_Callback(hObject, eventdata, handles)
% hObject    handle to check_polar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_polar


% --- Executes on selection change in pop_channel1.
function pop_channel1_Callback(hObject, eventdata, handles)
% hObject    handle to pop_channel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_channel1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_channel1
setTgl(hObject, eventdata, handles, handles.tgl_stim, 0);


% --- Executes during object creation, after setting all properties.
function pop_channel1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_channel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pop_channel2.
function pop_channel2_Callback(hObject, eventdata, handles)
% hObject    handle to pop_channel2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_channel2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_channel2
setTgl(hObject, eventdata, handles, handles.tgl_stim, 0);


% --- Executes during object creation, after setting all properties.
function pop_channel2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_channel2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_amp1_Callback(hObject, eventdata, handles)
% hObject    handle to txt_amp1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_amp1 as text
%        str2double(get(hObject,'String')) returns contents of txt_amp1 as a double
setTgl(hObject, eventdata, handles, handles.tgl_stim, 0);


% --- Executes during object creation, after setting all properties.
function txt_amp1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_amp1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_amp2_Callback(hObject, eventdata, handles)
% hObject    handle to txt_amp2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_amp2 as text
%        str2double(get(hObject,'String')) returns contents of txt_amp2 as a double
setTgl(hObject, eventdata, handles, handles.tgl_stim, 0);


% --- Executes during object creation, after setting all properties.
function txt_amp2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_amp2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_width1_Callback(hObject, eventdata, handles)
% hObject    handle to txt_width1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_width1 as text
%        str2double(get(hObject,'String')) returns contents of txt_width1 as a double
setTgl(hObject, eventdata, handles, handles.tgl_stim, 0);


% --- Executes during object creation, after setting all properties.
function txt_width1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_width1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_width2_Callback(hObject, eventdata, handles)
% hObject    handle to txt_width2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_width2 as text
%        str2double(get(hObject,'String')) returns contents of txt_width2 as a double
setTgl(hObject, eventdata, handles, handles.tgl_stim, 0);


% --- Executes during object creation, after setting all properties.
function txt_width2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_width2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_interphase_Callback(hObject, eventdata, handles)
% hObject    handle to txt_interphase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_interphase as text
%        str2double(get(hObject,'String')) returns contents of txt_interphase as a double
setTgl(hObject, eventdata, handles, handles.tgl_stim, 0);


% --- Executes during object creation, after setting all properties.
function txt_interphase_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_interphase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_MaxStimFreq_Callback(hObject, eventdata, handles)
% hObject    handle to txt_MaxStimFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_MaxStimFreq as text
%        str2double(get(hObject,'String')) returns contents of txt_MaxStimFreq as a double
handles.stimMaxFreq = eval(get(hObject, 'String'));
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function txt_MaxStimFreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_MaxStimFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tgl_stim.
function tgl_stim_Callback(hObject, eventdata, handles)
% hObject    handle to tgl_stim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tgl_stim

if get(hObject, 'Value') == 1
    % start stimulus 

    try

    handles.HardwareFuncs.CheckConnectionStimulator();

    handles.stimMaxFreq = eval(get(handles.txt_MaxStimFreq, 'String'));

    handles.StimSetupArgs = stimGetSetupArgs(handles);

    if handles.StimTriggerMode
        stimulator = handles.HardwareFuncs.SetStimTriggerMode([], handles.StimSetupArgs);
    end

    handles.StimActive = true;
    set(hObject, 'String', 'Stim On'); 

    catch ME
        hObject.Value = 0;
        handles.StimActive = false;
        guidata(hObject, handles);
        getReport(ME)
        errordlg(ME.message, 'Stim Setup Issue');
    end

else
    % stop stimulus 
    if handles.StimTriggerMode
        stimulator = handles.HardwareFuncs.ShutdownStimulator(handles.StimSetupArgs);
    end
    handles.StimActive = false;
    set(hObject, 'String', 'Stim Off');
end

guidata(hObject, handles)
settingChange(hObject);


% --- Executes during object creation, after setting all properties.
function tgl_stim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tgl_stim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
hObject.Value = 0; 
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function tgl_StartStop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tgl_StartStop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
hObject.Value = 0; 
guidata(hObject, handles)


% --- Executes on button press in check_artifact.
function check_artifact_Callback(hObject, eventdata, handles)
% hObject    handle to check_artifact (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_artifact
handles.check_artifact_Value = get(hObject, 'Value'); 
guidata(hObject, handles);
settingChange(hObject)


% --- Executes on selection change in pop_elecgrid.
function pop_elecgrid_Callback(hObject, eventdata, handles)
% hObject    handle to pop_elecgrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_elecgrid contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_elecgrid
contents = cellstr(get(hObject,'String'));
sel = contents{get(hObject,'Value')};
try
    [handles.showElecGrid, handles.elecGridFunc] = ...
        helperGUIv1_ElectrodeGridFuncSelect(sel, ...
        handles.FilterSetUp, handles.locutoff, handles.hicutoff);
catch ME4 
    getReport(ME4)
    errordlg(ME4.message, 'Electrode Grid Selection Issue');
    handles.showElecGrid = false;
    guidata(hObject, handles);
    pause(.01);
end
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function pop_elecgrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_elecgrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_gridmin_Callback(hObject, eventdata, handles)
% hObject    handle to txt_gridmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_gridmin as text
%        str2double(get(hObject,'String')) returns contents of txt_gridmin as a double
handles.ax_elecgrid.CLim(1) = eval(get(hObject, 'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txt_gridmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_gridmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_gridmax_Callback(hObject, eventdata, handles)
% hObject    handle to txt_gridmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_gridmax as text
%        str2double(get(hObject,'String')) returns contents of txt_gridmax as a double
handles.ax_elecgrid.CLim(2) = eval(get(hObject, 'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txt_gridmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_gridmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_griddur_Callback(hObject, eventdata, handles)
% hObject    handle to txt_griddur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_griddur as text
%        str2double(get(hObject,'String')) returns contents of txt_griddur as a double
settingChange(hObject)

% --- Executes during object creation, after setting all properties.
function txt_griddur_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_griddur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_remchan.
function push_remchan_Callback(hObject, eventdata, handles)
% hObject    handle to push_remchan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[remind, ok] = listdlg("ListString",handles.channelList, ...
    "PromptString", 'REMOVE these channels:');
if ok
    channelIDlist = handles.channelIDlist;
    remtf = false(size(channelIDlist)); remtf(remind) = true;
    handles.allChannelIDs = channelIDlist(~remtf);
    guidata(hObject, handles);
    settingChange(hObject);
end


% --- Executes on button press in push_stimCalibrate.
function push_stimCalibrate_Callback(hObject, eventdata, handles)
% hObject    handle to push_stimCalibrate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.RunMainLoop
    error('Cannot calibrate while running data aquisition. Press Stop to proceed.')
end
if handles.StimActive
    error('Stimulator is already active. Press Stim Off to proceed.')
end
handles.HardwareFuncs.CheckConnectionStimulator();
handles.StimSetupArgs = stimGetSetupArgs(handles);
handles.StimulatorLagTime = handles.HardwareFuncs.CalibrateStimulator(...
    rmfield(handles, handles.rmfieldList), false, ...
    handles.HardwareFuncs.SetupRecording, handles.HardwareFuncs.ShutdownRecording, ...
    handles.HardwareFuncs.GetTime, handles.HardwareFuncs.InitRawData, handles.HardwareFuncs.GetNewRawData, ...
    handles.HardwareFuncs.SetupStimulator, handles.HardwareFuncs.ShutdownStimulator, ...
    handles.HardwareFuncs.PulseStimulator);
guidata(hObject, handles);


% --- Executes on selection change in pop_channel3.
function pop_channel3_Callback(hObject, eventdata, handles)
% hObject    handle to pop_channel3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_channel3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_channel3
setTgl(hObject, eventdata, handles, handles.tgl_stim, 0);


% --- Executes during object creation, after setting all properties.
function pop_channel3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_channel3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pop_channel4.
function pop_channel4_Callback(hObject, eventdata, handles)
% hObject    handle to pop_channel4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_channel4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_channel4
setTgl(hObject, eventdata, handles, handles.tgl_stim, 0);


% --- Executes during object creation, after setting all properties.
function pop_channel4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_channel4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pop_channel5.
function pop_channel5_Callback(hObject, eventdata, handles)
% hObject    handle to pop_channel5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_channel5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_channel5
setTgl(hObject, eventdata, handles, handles.tgl_stim, 0);


% --- Executes during object creation, after setting all properties.
function pop_channel5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_channel5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
