function varargout = GUI_PD_v1b(varargin)
% GUI_PD_V1B MATLAB code for GUI_PD_v1b.fig
%      GUI_PD_V1B, by itself, creates a new GUI_PD_V1B or raises the existing
%      singleton*.
%
%      H = GUI_PD_V1B returns the handle to a new GUI_PD_V1B or the handle to
%      the existing singleton*.
%
%      GUI_PD_V1B('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_PD_V1B.M with the given input arguments.
%
%      GUI_PD_V1B('Property','Value',...) creates a new GUI_PD_V1B or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_PD_v1b_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_PD_v1b_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_PD_v1b

% Last Modified by GUIDE v2.5 05-Dec-2025 10:40:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_PD_v1b_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_PD_v1b_OutputFcn, ...
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


% --- Executes just before GUI_PD_v1b is made visible.
function GUI_PD_v1b_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_PD_v1b (see VARARGIN)

% Choose default command line output for GUI_PD_v1b
handles.output = hObject;

% --- Begin My Code ---

% change ax_polar to polar axes 
axes(handles.ax_polar);
polarhistogram(pi/2);
handles.ax_polar = gca;

% Create a timer object that will be used to grab data and refresh
% analysis/plotting
%{
handles.timer = timer(...
    'ExecutionMode', 'fixedSpacing', ...       % Run timer repeatedly
    'Period', 0.01, ...                      
    'TimerFcn', {@updateDisplay,hObject}, ... % callback function.  Pass the figure handle
    'StartFcn', {@StartMainLoop,hObject}, ... % callback to execute when timer starts
    'StopFcn',  {@StopMainLoop,hObject}, ...
    'ErrorFcn', {@TimerError,hObject}, ...
    'BusyMode', 'error'); 
%}

% serial user data
ud = struct('ReceivedData', '', ...
            'ParadigmPhase', 'Stopped');
RecSrlCallback = @(hsrl,evt) CharSerialCallbackReceiver_PD_v0(hsrl,evt, ...
                    handles.textSrl, handles.txt_Status);

% save location
svloc = ['Saved Data PD',filesep,'Saved Data ',...
    datestr(datetime, 'yyyy-mm-dd HH.MM.SS')];
pause(1)
mkdir(svloc); 

handles = helperGUIv0_OpeningInitialize(handles, ud, svloc, RecSrlCallback);

% initiate other vars ...
handles.DAQstatus = handles.cbmexStatus;
handles.bufferSize     = 10; 
handles.bufferSizeGrid = 10;
handles.allChannelIDs = [];
handles.channelIndex = [];

% additional phase tracking buffers & objects
emptyStorage = nan(100000,1);
handles.redStorage1 = emptyStorage; handles.redP1 = 1;
handles.redStorage2 = emptyStorage; 
handles.yelStorage1 = emptyStorage; handles.yelP1 = 1;
handles.yelStorage2 = emptyStorage; 
handles.grnStorage1 = emptyStorage; handles.grnP1 = 1;
handles.grnStorage2 = emptyStorage; 
handles.stpStorage1 = emptyStorage; handles.stpP1 = 1;
handles.stpStorage2 = emptyStorage; 
handles.redDataBuffer = []; handles.h_redTrace = [];
handles.yelDataBuffer = []; handles.h_yelTrace = [];
handles.grnDataBuffer = []; handles.h_grnTrace = [];
handles.stpDataBuffer = []; handles.h_stpTrace = [];

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

% default values 
PDSwin = str2double(get(handles.txt_PDSwin,'String'));
PDSwin = ceil(PDSwin*1000); handles.PDSwin1 = PDSwin;
handles.PDSwin2 = ceil(.02*PDSwin); 
handles.bufferSize = str2double(get(handles.txt_display,'String')) * 1000;
handles.bufferSizeGrid = str2double(get(handles.txt_griddur,'String')) * 1000;
handles.stimMaxFreq = eval(get(handles.txt_MaxStimFreq, 'String'));
%handles.check_artifact_Value = false;
handles.ControllerResult = 0;

handles.QueuedStim = timer(...
    'StartDelay', 10, ...
    'TimerFcn',   {@myPULSE, hObject}, ...
    'StopFcn',    {@finishPULSE, hObject}, ...
    'StartFcn',   {@schedulePULSE, hObject}, ...
    'UserData',   Inf);

% turn off graphs' interactivity 
disableDefaultInteractivity(handles.ax_polar);
disableDefaultInteractivity(handles.ax_raw);
disableDefaultInteractivity(handles.ax_filt);
disableDefaultInteractivity(handles.ax_timing);
disableDefaultInteractivity(handles.ax_elecgrid);

% hardware-specific functions 
[handles.HardwareFuncs, handles.StimTriggerMode] = helperGUIv1_DefHardwareFuncs();
handles.initTic = tic;

% Update handles structure
guidata(hObject, handles);

clc
%clear all

% UIWAIT makes GUI_PD_v1b wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_PD_v1b_OutputFcn(hObject, eventdata, handles) 
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

handles = helperGUIv1b_DAQopen(handles);

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
handles.HardwareFuncs.ShutdownRecording();
handles.DAQstatus = false;
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
% This has been added for debugging. Most likely occurs when right click
% instead of left. 


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
        errordlg('No DAQ connection. Open connection before starting','Not Connected')
        return
    end
    
    set(hObject,'String','Stop');

    % This starts the timer and also executes the StartFnc which grabs the
    % data, creates the buffer and plots the first bit of data
    % start(handles.timer)
    StartMainLoop(hObject, eventdata, handles);
    
% Stop
else
    set(hObject,'String','Start')
    % stop(handles.timer)
    StopMainLoop(hObject, eventdata, handles);
end

% --- Executes on selection change in pop_channels.
function pop_channels_Callback(hObject, eventdata, handles)
% hObject    handle to pop_channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_channels contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_channels
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
    stop(handles.QueuedStim)
    handles.HardwareFuncs.ShutdownStimulator(handles.stimulator);
    handles.StimActive = false;
    guidata(hObject, handles)
    setTgl(hObject, eventdata, handles, handles.tgl_StartStop, 0);
end
catch ME1
    getReport(ME1)
end

try
% save stored data 
phStorage = {handles.pkStorage1, handles.trStorage1, ...
    handles.redStorage1, handles.yelStorage1, handles.grnStorage1, handles.stpStorage1};
for iph = 1:length(handles.PhaseOfInterest)
    phname = handles.PhaseOfInterestName(iph);
    phname = phname+"Time";
    phStorage1 = phStorage{iph};
    if numel(phStorage1)
        if sum(~isnan(phStorage1))
            %assignin("caller",phname,phStorage1);
            eval(phname+" = phStorage1;");
            svfn = [handles.SaveFileLoc,filesep,'SaveFile',num2str(handles.SaveFileN),'.mat'];
            disp("Saving "+phname+" to "+svfn)
            save(svfn,phname);
            handles.SaveFileN = handles.SaveFileN + 1;
        end
    end
end
StimTime = handles.stStorage1;
if numel(StimTime)
    if sum(~isnan(StimTime))
        svfn = [handles.SaveFileLoc,filesep,'SaveFile',num2str(handles.SaveFileN),'.mat'];
        disp(['Saving Stimulus to ',svfn])
        save(svfn,'StimTime');
        handles.SaveFileN = handles.SaveFileN + 1;
    end
end
SerialLog = handles.srlStorage1;
svfn = [handles.SaveFileLoc,filesep,'SaveFile',num2str(handles.SaveFileN),'.mat'];
disp(['Saving Serial to ',svfn])
save(svfn,'SerialLog');
disp('Saving all data...')
ConsolidateSavedData(handles.SaveFileLoc)
catch ME2
    getReport(ME2)
end

try
handles.HardwareFuncs.ShutdownRecording();
handles.DAQstatus = false;
guidata(hObject, handles)
%stop(handles.timer)
StopMainLoop(hObject,eventdata,handles)
%delete(handles.timer)
if ~handles.SerialArgs.NoSerial
    delete(handles.srl);
end
catch ME3
    getReport(ME3)
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
while handles.RunMainLoop

try

    % ensure connection with hardware 
    handles = guidata(hObject);
    if ~handles.DAQstatus
        %stop(handles.timer)
        StopMainLoop(hObject,eventdata,handles)
    end

    % timing 
    tNow = datetime;
    %timeDisp2 = handles.timeDisp0 + toc(handles.timeDisp1);
    %handles.timeDispBuff = bufferData(handles.timeDispBuff, timeDisp2);
    %guidata(hObject, handles);

    % main iteration ================================================
    selRaw2Flt = handles.selInds.selRaw2Flt;
    selFlt2For = handles.selInds.selFlt2For;
    selFor2Art = handles.selInds.selFor2Art;
    selRaw2Art = handles.selInds.selRaw2Art;
    selRaw2For = handles.selInds.selRaw2For;

    rawData = handles.recDataStructs.rawD;
    artRemData = handles.recDataStructs.artD;
    fltData = handles.recDataStructs.fltD;
    forData = handles.recDataStructs.forD;
    timeBuffs = handles.recDataStructs.timeBuffs;
    initTic = handles.recDataStructs.initTic;
    forBuffs = handles.recDataStructs.forBuffs;

    % enable defined filtering/forecasting funcs only if ready
    if handles.FilterSetUp
        filtfun = @filtFun;
    else
        filtfun = [];
    end
    if handles.MdlSetUp
        forefun = @foreFun;
    else
        forefun = [];
    end
    if handles.check_artifact.Value && handles.MdlSetUp
        artremfun = @artRemFun;
    else
        artremfun = [];
    end

    [...
    timeBuffs, rawData, ...
    artRemData, handles.artRemArgs, ...
    fltData, handles.filtArgs, ...
    forBuffs, forData, handles.foreArgs] = ...
    iterReadBrain(...
        timeBuffs, rawData, @() handles.HardwareFuncs.GetNewRawData(handles.allChannelIDs), ...
        selRaw2Art, selFor2Art, selRaw2Flt, selRaw2For, selFlt2For, ...
        artRemData, artremfun, handles.artRemArgs, ...
        fltData, filtfun, handles.filtArgs, ...
        forBuffs, forData, forefun, handles.foreArgs);

    handles.recDataStructs.rawD = rawData;
    handles.recDataStructs.artD = artRemData;
    handles.recDataStructs.fltD = fltData;
    handles.recDataStructs.forD = forData;
    handles.recDataStructs.timeBuffs = timeBuffs;
    handles.recDataStructs.forBuffs = forBuffs;
    % ===============================================================

    lastSampleProcTime = handles.recDataStructs.timeBuffs{handles.channelIndex}(end);
    selFlt = handles.selInds.selFlt2For;
    if isempty(selFlt) || isnan(selFlt)
        selFlt = 1;
    end
    fltPlt = handles.recDataStructs.fltData{4,selFlt}{:,2};

    % stimulation ===================================================
    Stim2Q = false;
    forBuff = forBuffs{1}; 
    forBuffNew = [max(forBuff(:,1)), max(forBuff(:,2))]; 
    forBuffNew = forBuffNew - timeBuffs{chInd}(end,:); % [t2p, t2t]
    StimController = handles.ControllerResult;
    doStim = ((~isempty(StimController)) && handles.StimActive) && (handles.FilterSetUp && handles.MdlSetUp);
    if doStim && (StimController > 0)
        Stim2Q = true;
        t2Q = forBuffNew(:,StimController);
    end
    if Stim2Q && (t2Q >= 0)
        t2Q = .001*floor(1000*t2Q); % round to nearest 1ms 
        if t2Q < (100/handles.locutoff + handles.TimeShiftFIR)
            t2Qabs = t2Q + lastSampleProcTime; % in NSP time "absolute"
            Dt2Q = t2Qabs - handles.stimLastTime; 
            if 1/Dt2Q <= handles.stimMaxFreq
                stim2Q_proceed = true;
                if strcmp(handles.QueuedStim.Running, 'on')
                    % last queued stim has not yet fired 
                    if t2Qabs > handles.QueuedStim.UserData
                        % new requested point is later than current timer
                        stim2Q_proceed = t2Q > handles.StimulatorLagTime; 
                            % is there enough time to make a change
                        stim2Q_proceed = stim2Q_proceed && ...
                            (t2Qabs - handles.QueuedStim.UserData) > handles.StimulatorLagTime; 
                            % is the change outside margin of error
                        stim2Q_proceed = stim2Q_proceed && ...
                            (t2Qabs - handles.QueuedStim.UserData) < 1/handles.hicutoff; 
                            % is it trying to target the next cycle
                    end
                end
                if stim2Q_proceed
                    if strcmp(handles.QueuedStim.Running, 'on')
                        stop(handles.QueuedStim);
                    end
                    handles.QueuedStim.StartDelay = t2Q;
                    handles.QueuedStim.UserData = t2Qabs;
                    start(handles.QueuedStim);
                end
            end
        end
    else
        if strcmp(handles.QueuedStim.Running, 'on')
            stop(handles.QueuedStim);
        end
    end
    % ===============================================================
    
    if handles.check_polar.Value
    % x axes alignment 
    xTbl = data2timetable(rawPlt, [], []); % FIX THIS 
    tPltRng = xTbl.Time;
    common_xlim = tPltRng - tNow;
    common_xdiff = diff(common_xlim); 
    ext_xdiff = common_xdiff * handles.ax_filt.InnerPosition(3) / ...
        handles.ax_raw.InnerPosition(3); 
    ext_xlim = [0, ext_xdiff] + common_xlim(1); % align left 
    end

    % update serial log 
    % TO DO: there should be a better way to do this; serial callback
    % should trigger an event or listener that logs the info 
    [handles, newsrl] = helperGUIv0_UpdateSerialLog(handles);
    if newsrl
        if handles.FilterSetUp
            if numel(fltPlt) % should this be necessary when above is met?
                ControllerLastResult = handles.ControllerResult;
                ControllerResult = Controller_PDS_PD( handles.srl, handles, ...
                    fltPlt((end-handles.PDSwin1+1):end) );
                if ControllerResult ~= ControllerLastResult
                    handles.ControllerResult = ControllerResult;
                    guidata(hObject, handles);
                end
            end
        end
    end

    % update electrode grid 
    if handles.showElecGrid
        try
            handles.elecGridImg.CData = helperGUIv1_ElectrodeGridUpdate(...
                handles.elecGridImg, handles.elecGridFunc, ...
                handles.channelIDlist, handles.allChannelIDs, rawData(4,:), ...
                handles.bufferSizeGrid, handles.fSamples);
        catch ME4 
            getReport(ME4)
            errordlg(ME4.message, 'Electrode Grid Issue');
            handles.showElecGrid = false;
            guidata(hObject, handles);
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
    set(handles.h_rawDataTrace,'YData',rawPlt);
    set(handles.h_timingTrace,'YData',[nan; diff(handles.recDataStructs.timeBuffs{1})]);
    if handles.check_polar.Value
        % also update time (x) axis
    end
    guidata(hObject,handles)

catch ME 
    getReport(ME)
    guidata(hObject, handles)
    %stop(handles.timer)
    %keyboard
end

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
        helperGUIv1b_plotSetupRaw(handles, tNow);

    % initiate filtered data plot
    if handles.FilterSetUp
        try        
            handles = helperGUIv1_plotSetupFltMdl(handles, tNow, ...
                fltPlt, forPlt, forBuff, tSt, common_xlim, unitname);
            % setup additional phase-tracking buffers & plots
            if handles.MdlSetUp
                tR = forBuff(:,3); tY = forBuff(:,4); tG = forBuff(:,5); tS = forBuff(:,6);
                axes(handles.ax_filt)
                handles.h_redTrace = plot(handles.time0 + seconds(tR) - tNow, zeros(size(tR)), ...
                    's', 'Color',"#A2142F");
                handles.h_yelTrace = plot(handles.time0 + seconds(tY) - tNow, zeros(size(tY)), ...
                    's', 'Color',"#EDB120");
                handles.h_grnTrace = plot(handles.time0 + seconds(tG) - tNow, zeros(size(tG)), ...
                    's', 'Color',"#77AC30");
                handles.h_stpTrace = plot(handles.time0 + seconds(tS) - tNow, zeros(size(tS)), ...
                    's', 'Color',"#7E2F8E");
            end
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

    % define input args for filtering/forecasting funcs 
    selRaw2Flt = handles.selInds.selRaw2Flt;
    selFlt2For = handles.selInds.selFlt2For;
    selFor2Art = handles.selInds.selFor2Art;
    selRaw2Art = handles.selInds.selRaw2Art;
    selRaw2For = handles.selInds.selRaw2For;
    Fs = handles.fSamples;
    if FilterSetUp
        fIC = arrayfun(@(ord) zeros(ord,1), handles.FilterOrder, 'UniformOutput',false);
        filtArgs.fltInit = fIC; filtArgs.fltObj = {1; handles.BPF};
        filtArgs.TimeShift = handles.FilterOrder(1)/(2*Fs(selRaw2Flt));
        if MdlSetUp
            foreArgs.K = handles.PDSwin1; foreArgs.k = handles.PDSwin2;
            foreArgs.TimeStart = nan(size([selRaw2For, selFlt2For]));
            foreArgs.TimeShift = [zeros(size(selRaw2For)), filtArgs.TimeShift(selFlt2For)];
            foreArgs.ARmdls = {handles.Mdl};
            foreArgs.SampleRates = Fs(selRaw2Flt); % !!! needs improvement
            foreArgs.FreqRange = [handles.locutoff, handles.hicutoff];
            foreArgs.PhaseOfInterest = handles.PhaseOfInterest;
            foreArgs.StimulatorLagTime = handles.StimulatorLagTime;
            artRemArgs.SampleRates = Fs(selRaw2Art);
            artRemArgs.StimDur = .11; % seconds !!! MAKE ADJUSTABLE
            artRemArgs.StimTimes = cell(size(artRemArgs.SampleRates));
            artRemArgs.nOverlap = zeros(size(artRemArgs.SampleRates));
        else
            foreArgs = [];
            artRemArgs = [];
        end
    else
        filtArgs = [];
        foreArgs = [];
        artRemArgs = [];
    end
    handles.filtArgs = filtArgs; handles.foreArgs = foreArgs; handles.artRemArgs = artRemArgs;

    % handles = disconnectSerial(handles);
    handles.RunMainLoop = true; 
    guidata(hObject,handles)
    updateDisplay(hObject,eventdata)
    
catch ME
    getReport(ME)
    handles.RunMainLoop = false; 
    guidata(hObject, handles);
    %stop(handles.timer);
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
    if isfield(handles, 'StimActive')
        if handles.StimActive
            stop(handles.QueuedStim)
        end
    end
    if isfield(handles, 'RunMainLoop')
        handles.RunMainLoop = false;
    end
    guidata(hObject, handles)
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
%{
if strcmp(handles.timer.Running,'on')
    stop(handles.timer)
    start(handles.timer)
end
%}
if handles.RunMainLoop
    StopMainLoop(hObject,[],handles)
    StartMainLoop(hObject,[],handles)
end


% --- GUI Helpers ---

% function handles = initRec(handles)
% TO DO: (re-)init RecDataStructs here only 
% can be called by settingchange and only called by startmainloop at the
% beginning 

function setTgl(hObject, eventdata, handles, hTgl, newValue)
% set a toggle button to a desired Value and activate its callback if it is
% not currently at that value. 
curValue = hTgl.Value; 
hTgl.Value = newValue; guidata(hObject, handles);
if ~(curValue == newValue)
    hTgl.Callback(hTgl, eventdata);
    guidata(hObject, handles);
end

% --- PhaseDetect helpers ---

function [artRemTails, artRemArgs] = ...
        artRemFun(artRemArgs, rawTails, forTails)
% forTails must be the forecast data starting at the same time as rawTails
artRemTails = cell(size(rawTails));
FsArt = artRemArgs.SampleRates;
StimDur = artRemArgs.StimDur; % seconds to remove
StimLen = ceil(StimDur.*FsArt); % #samples to remove 
StimTimesTail = artRemArgs.StimTimes; 
for ch_art = 1:size(rawTails,2)
    tXfor = forTails{ch_art};
    tX = rawTails{ch_art}; % [time, data]
    stimtimes = StimTimesTail{ch_art}; % time to stim (sec)
    stiminds = round(stimtimes * FsArt(ch_art));
    stiminds = stiminds(stiminds > 0);
    for i1 = stiminds'
        i2 = i1 + StimLen(ch_art);
        if i2 > height(tX)
            % artifact will carry over into the next packet 
            nO = i2 - height(tX);
            artRemArgs.nOverlap(ch_art) = nO;
            tX = [tX; nan(nO,2)]; 
        else
            % artifact limited to this packet 
            artRemArgs.nOverlap(ch_art) = 0;
        end
        i2 = min(i2, height(tX)); % why is this necessary??
        tX(i1:i2,2) = tXfor(i1:i2,2);
    end
    artRemTails{ch_art} = tX;
end

function [foreTails, foreBuffsAdd, foreArgs] = foreFun(foreArgs, inData)
% inData should just be the (filtered) channel(s) of interest
% foreTails will be the forecast tails starting at the current time minus
% any filtering delay
if isempty(foreArgs)
    error('Attempted forecast function with empty arguments.')
end
K = foreArgs.K; k = foreArgs.k;
ARmdls = foreArgs.ARmdls;
fs = foreArgs.SampleRates;
Ts = foreArgs.TimeShift;
Fco = foreArgs.FreqRange;
phis = foreArgs.PhaseOfInterest;
TstimLag = foreArgs.StimulatorLagTime;
foreTails = cell(size(inData)); foreBuffsAdd = foreTails; 
for ch_fore = 1:size(inData,2)
    armdl = ARmdls{ch_fore};
    FT = myFastForecastAR(armdl, inData{ch_fore}(:,2), K);
    foreTails{ch_fore} = [nan(height(FT),1), FT];
    foreTails{ch_fore}(1,1) = foreArgs.TimeStart(ch_fore);

    FT = FT(1:k,:); % use limited duration for hilbert padding
    %[t2,i2,phi_inst,f_inst] = blockPDS(...
    t2 = blockPDS(...
        inData{ch_fore}(:,2), FT, fs(ch_fore), phis, ...
        Ts(ch_fore)+TstimLag, Fco(1), Fco(2));
    t2 = t2-Ts(ch_fore); % [t2peak, t2trough]
    % t2 = max(t2,0); % Should this be necessary? Will this cause problems?
    t2(t2 < 0) = nan;
    foreBuffsAdd{ch_fore} = t2; 
end

function [fltTails, fltArgs] = filtFun(fltArgs, rawTails)
filtObj = fltArgs.fltObj;
filtInit = fltArgs.fltInit;
filtFin = cell(size(filtInit));
fltTails = cell(size(rawTails));
for ch_flt = 1:size(rawTails,2)
    a = filtObj{1,ch_flt}; b = filtObj{2,ch_flt};
    [FT,filtFin{ch_flt}] = filter(b,a,rawTails{ch_flt}(:,2:end),filtInit{ch_flt});
    fltTails{ch_flt} = [rawTails{ch_flt}(:,1) - fltArgs.TimeShift(ch_flt), FT];
end
fltArgs.fltInit = filtFin;

% ----------------------------------------------------------------------- %
% ----                                                                --- %
% ----                        New Callbacks                           --- %
% ----                                                                --- %
% ----------------------------------------------------------------------- %


% --- Executes on selection change in pop_RedStim.
function pop_RedStim_Callback(hObject, eventdata, handles)
% hObject    handle to pop_RedStim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_RedStim contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_RedStim
contents = cellstr(get(hObject,'String'));
RedStim = contents{get(hObject,'Value')};
handles.StimMode.red = RedStim; 
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function pop_RedStim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_RedStim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
pop_RedStim_Callback(hObject, eventdata, handles)


% --- Executes on selection change in pop_YellowStim.
function pop_YellowStim_Callback(hObject, eventdata, handles)
% hObject    handle to pop_YellowStim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_YellowStim contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_YellowStim
contents = cellstr(get(hObject,'String'));
YellowStim = contents{get(hObject,'Value')};
handles.StimMode.yellow = YellowStim; 
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function pop_YellowStim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_YellowStim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
pop_YellowStim_Callback(hObject, eventdata, handles)


% --- Executes on selection change in pop_GreenStim.
function pop_GreenStim_Callback(hObject, eventdata, handles)
% hObject    handle to pop_GreenStim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_GreenStim contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_GreenStim
contents = cellstr(get(hObject,'String'));
GreenStim = contents{get(hObject,'Value')};
handles.StimMode.green = GreenStim; 
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function pop_GreenStim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_GreenStim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pop_StopStim.
function pop_StopStim_Callback(hObject, eventdata, handles)
% hObject    handle to pop_StopStim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_StopStim contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_StopStim
contents = cellstr(get(hObject,'String'));
StopStim = contents{get(hObject,'Value')};
handles.StimMode.Stopped = StopStim; 
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function pop_StopStim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_StopStim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
pop_StopStim_Callback(hObject, eventdata, handles)



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

%stop(handles.timer); % to avoid dataQueue overflow 

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
StopMainLoop(hObject,eventdata,handles)
guidata(hObject, handles)
pause(.01);
StartMainLoop(hObject,eventdata,handles)
%start(handles.timer);



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

%stop(handles.timer); % to avoid dataQueue overflow 

n = str2double(get(handles.txt_AR,'String'));
N = str2double(get(handles.txt_PDSwin,'String'));
PDSwin = ceil(N*handles.fSample); handles.PDSwin1 = PDSwin;
handles.PDSwin2 = ceil(.02*PDSwin); 

try

    y = handles.recDataStructs.fltD{4,1}{:,2};
    handles = helperGUIv0_pushAR(handles, PDSwin, n, y);
    
% restart timer and plots
pause(.01);
StopMainLoop(hObject,eventdata,handles)
guidata(hObject, handles)
pause(.01);
StartMainLoop(hObject,eventdata,handles)
%start(handles.timer);

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
 
    if handles.StimActive
        % there is a problem if StimActive is already enabled 
        error('Issue last time turning off stim; please try again.')
    end

    handles.HardwareFuncs.CheckConnectionStimulator();

    handles.stimMaxFreq = eval(get(handles.txt_MaxStimFreq, 'String'));

    handles.StimSetupArgs = stimGetSetupArgs(handles);

    if handles.StimTriggerMode
        stimulator = handles.HardwareFuncs.SetupStimulator(handles.StimSetupArgs);
        stimulator = handles.HardwareFuncs.SetStimTriggerMode(stimulator);
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
%handles.check_artifact_Value = get(hObject, 'Value'); 
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
    handles.HardwareFuncs.PulseStimulator, ...
    handles.HardwareFuncs.SetupStimTTL);
guidata(hObject, handles);
