function varargout = GUI_PD_v0(varargin)
% GUI_PD_V0 MATLAB code for GUI_PD_v0.fig
%      GUI_PD_V0, by itself, creates a new GUI_PD_V0 or raises the existing
%      singleton*.
%
%      H = GUI_PD_V0 returns the handle to a new GUI_PD_V0 or the handle to
%      the existing singleton*.
%
%      GUI_PD_V0('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_PD_V0.M with the given input arguments.
%
%      GUI_PD_V0('Property','Value',...) creates a new GUI_PD_V0 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_PD_v0_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_PD_v0_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_PD_v0

% Last Modified by GUIDE v2.5 10-Jul-2025 01:39:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_PD_v0_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_PD_v0_OutputFcn, ...
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


% --- Executes just before GUI_PD_v0 is made visible.
function GUI_PD_v0_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_PD_v0 (see VARARGIN)

% Choose default command line output for GUI_PD_v0
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
    'Period', 0.01, ...                      % Initial period is 100 ms
    'TimerFcn', {@updateDisplay,hObject}, ... % callback function.  Pass the figure handle
    'StartFcn', {@startTimer,hObject}, ...
    'StopFcn',  {@stopTimer,hObject});     % callback to execute when timer starts
%}

% serial user data
ud = struct('ReceivedData', '', ...
            'ParadigmPhase', 'Stopped');

% save location
svloc = ['Saved Data PD',filesep,'Saved Data ',...
    datestr(datetime, 'yyyy-mm-dd HH.MM.SS')];
pause(1)
mkdir(svloc); 

handles = helperGUIv0_OpeningInitialize(handles, ud, svloc);

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

% Update handles structure
guidata(hObject, handles);

clc
%clear all

% UIWAIT makes GUI_PD_v0 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_PD_v0_OutputFcn(hObject, eventdata, handles) 
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

% Use a TRY-CATCH in case cbmex is already open.  If you try to open it
% when it's already open, Matlab throws a 'MATLAB:unassignedOutputs'
% MException.
try
    cbmex('open');
catch ME
    if strcmp(ME.identifier,'MATLAB:unassignedOutputs')
        % Dont need to do anything because cbmex is already open and it
        % already sends a message stating that
    else
        disp(ME)
    end
end
handles.cbmexStatus = true;

cbmex('trialconfig',1,'absolute','double')
pause(0.1)

% Acquire some data to get channel information.  Determine which channels
% are enabled
[spikeEvents, time, continuousData] = cbmex('trialdata',1);
handles.channelList = [continuousData{:,1}];
chL = [handles.channelList, nan];
% set channel popup meno to hold channels
set(handles.pop_channels, 'String', handles.channelList);
for pop_ = [handles.pop_channel1, ...
            handles.pop_channel2, ...
            handles.pop_channel3, ...
            handles.pop_channel4, ...
            handles.pop_channel5]
    set(pop_, 'String', chL);
end
% electrode grid 
axes(handles.ax_elecgrid)
ncol = 3; % # of columns
nrow = 21;
img = nan(nrow, ncol);
[X,Y] = meshgrid(1:ncol, 1:nrow);
handles.elecGridImg = imagesc(img); colormap('parula'); colorbar; hold on;
for ch = 1:min(length(handles.channelList),63)
    chname = handles.channelList(ch); chname = num2str(chname);
    text(X(ch),Y(ch), chname, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'FontWeight', 'bold', ...
        'Color',[1 0 0]);
end

% Set the Start/Stop toggle button to stopped state (String is 'Start' and
% Value is 1)
set(handles.tgl_StartStop,'String','Start', 'Value',0)

guidata(hObject,handles)

% --- Executes on button press in cmd_cbmexClose.
function cmd_cbmexClose_Callback(hObject, eventdata, handles)
% hObject    handle to cmd_cbmexClose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cbmex('close')
handles.cbmexStatus = false;
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
% This has been added for debugging. Does the code ever get here?


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
    if ~handles.cbmexStatus
        errordlg('No cbmex connection.  Open connection before starting','Not Connected')
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

try
if handles.StimActive
    stop(handles.QueuedStim)
    handles.StimActive = false;
    guidata(hObject, handles)
    setTgl(hObject, eventdata, handles, handles.tgl_StartStop, 0);
end
catch ME1
    getReport(ME1)
end

try
% save stored data 
PeakTime = handles.pkStorage1;
svfn = [handles.SaveFileLoc,filesep,'SaveFile',num2str(handles.SaveFileN),'.mat'];
disp(['Saving Peaks to ',svfn])
save(svfn,'PeakTime');
handles.SaveFileN = handles.SaveFileN + 1;
TroughTime = handles.trStorage1;
svfn = [handles.SaveFileLoc,filesep,'SaveFile',num2str(handles.SaveFileN),'.mat'];
disp(['Saving Troughs to ',svfn])
save(svfn,'TroughTime');
handles.SaveFileN = handles.SaveFileN + 1;
StimTime = handles.stStorage1;
svfn = [handles.SaveFileLoc,filesep,'SaveFile',num2str(handles.SaveFileN),'.mat'];
disp(['Saving Stimulus to ',svfn])
save(svfn,'StimTime');
handles.SaveFileN = handles.SaveFileN + 1;
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
cbmex('close')
handles.cbmexStatus = false;
guidata(hObject, handles)
%stop(handles.timer)
StopMainLoop(hObject,eventdata,handles)
%delete(handles.timer)
delete(handles.srl);
if handles.StimActive
    stop(handles.QueuedStim)
    handles.stimulator.stop();
    handles.stimulator.disconnect;
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
function updateDisplay(hObject, eventdata)

handles = guidata(hObject);
while handles.RunMainLoop

try
    
    handles = guidata(hObject);
    if ~handles.cbmexStatus
        %stop(handles.timer)
        StopMainLoop(hObject,eventdata,handles)
    end
    
    [events, time, continuousData] = cbmex('trialdata',1);

    if ~isempty(continuousData)
    if size(continuousData,1) >= handles.channelIndex
        % !!!!! sometimes CD does not have all channels?????
        % also add a check whether channel name at channelIndex is the same
        % as specified in UI *****
        
    newContinuousData = continuousData{handles.channelIndex,3}; 
    N = length(newContinuousData);

    [handles, newContinuousData] = ...
        helperGUIv0_MainLoopPrepareNewData(handles, newContinuousData, time);

    guidata(hObject,handles)

    % update serial log 
    handles = helperGUIv0_UpdateSerialLog(handles);

    % update electrode grid 
    if handles.showElecGrid
        try
            elecimg = handles.elecGridImg.CData;
            for ch = 1:min(size(continuousData,1),63)
                newContinuousData_ch = continuousData{ch,3};
                fSample_ch = continuousData{ch,2};
                elecimg(ch) = handles.elecGridFunc(newContinuousData_ch, fSample_ch);
            end
            handles.elecGridImg.CData = elecimg;
        catch ME4 
            getReport(ME4);
            errordlg(ME4.message, 'Electrode Grid Issue');
            handles.showElecGrid = false;
            pause(.01);
        end
    end

    [handles, ~, phBuffers, phStorage] = ...
        helperGUIv0_MainProcessAndPlot(handles, ...
        newContinuousData, @Controller_PDS_PD, ...
        {handles.h_peakTrace, handles.h_trouTrace}, ...
        {handles.peakDataBuffer, handles.trouDataBuffer}, ...
        {handles.pkStorage1, handles.trStorage1; ...
         handles.pkP1,       handles.trP1; ...
         handles.pkStorage2, handles.trStorage2});

    % additional phase tracking buffers & plots

    handles.peakDataBuffer = phBuffers{1}; handles.trouDataBuffer = phBuffers{2}; 

    handles.pkStorage1   = phStorage{1,1}; handles.pkP1   = phStorage{2,1}; 
    handles.pkStorage2   = phStorage{3,1}; 
    handles.trStorage1   = phStorage{1,2}; handles.trP1   = phStorage{2,2}; 
    handles.trStorage2   = phStorage{3,2}; 

    % update YData of ax_rawData
    if ~(length(handles.rawDataBuffer) == handles.bufferSize)
        % malfunctioning artifact removal can cause this 
        % there will be a problem if it gets here 
        error('Raw data buffer changed size. Check artifact removal code.')
        % keyboard
    end
    set(handles.h_rawDataTrace,'YData',handles.rawDataBuffer)
    set(handles.h_timingTrace,'YData',handles.diffSampleProcTime)

    guidata(hObject,handles)
    pause(.001)
    drawnow 
    pause(.001)

    end
    end

catch ME 
    getReport(ME)
    keyboard
end

end
    
% ----------------------------------------------------------------------- %
% ----                                                                --- %
% ----                      Timer Start Function                      --- %
% ----                                                                --- %
% ----------------------------------------------------------------------- %

% Runs once when timer starts.  It initializes plot and buffer and
% accommodates any new selection by user.
function  StartMainLoop(hObject, eventdata, handles)
% Put the whole function in a try-catch block.  This makes debugging much
% easier because it captures the error and displays a report to the Matlab
% Command Window

try

    handles = helperGUIv0_StartMainLoop(handles);

    handles.RunMainLoop = true; 
    guidata(hObject,handles)
    updateDisplay(hObject,eventdata)
    
catch ME
    getReport(ME)
    keyboard
end

% ----------------------------------------------------------------------- %
% ----                                                                --- %
% ----                      Timer Stop Function                       --- %
% ----                                                                --- %
% ----------------------------------------------------------------------- %

% Runs once when timer stops.
function  StopMainLoop(hObject, eventdata, handles)
% Put the whole function in a try-catch block.  This makes debugging much
% easier because it captures the error and displays a report to the Matlab
% Command Window

try
    if handles.StimActive
        stop(handles.QueuedStim)
    end
    handles.RunMainLoop = false;
    guidata(hObject, handles)
catch ME
    getReport(ME)
    keyboard
end

% ----------------------------------------------------------------------- %
% ----                                                                --- %
% ----               Stimulus Pulse Timer Functions                   --- %
% ----                                                                --- %
% ----------------------------------------------------------------------- %

function myPULSE(hTimer,eventdata,hFigure)
handles = guidata(hFigure);
stimulator = handles.stimulator;
if ~stimulator.isConnected()
    warning('Stimulator is not connected.')
end
if stimulator.isLocked()
    warning('Stimulator is locked.')
end
stimstatus = stimulator.getSequenceStatus;
if stimstatus == 2
    warning('Stimulator is already playing.')
    % if this warning shows up, consider wrapping the rest of the function
    % in a conditional that stimstatus == 0 [stopped] or 1 [paused] (?), or
    % try using stimulator.stop()
end
% if stimstatus == 0
stimtime1 = cbmex('time');
stimulator.play(1); % consider changing to groupStimulus or manualStim to save time
stimtime2 = cbmex('time');
dstimtime = stimtime2 - stimtime1; 
% disp time of pulse using eventdata
eventTime = datestr(eventdata.Data.time);
stimtime = .5*(stimtime1 + stimtime2);
stimschedtime = hTimer.UserData; 
disp(['Stimulus pulsed at ',eventTime,' within ',num2str(dstimtime),'s, ',...
      num2str(stimtime - stimschedtime),' s late'])
handles.stimLastTime = stimtime; handles.stimNewTime = stimtime; 
if handles.stP1 <= length(handles.stStorage1)
    handles.stStorage1(handles.stP1) = stimtime; 
    handles.stP1 = handles.stP1 + 1;
else
    % storage full; save
    StimTime = handles.stStorage1;
    svfn = [handles.SaveFileLoc,filesep,'SaveFile',num2str(handles.SaveFileN),'.mat'];
    disp(['Saving Stimulus to ',svfn])
    save(svfn,'StimTime');
    handles.SaveFileN = handles.SaveFileN + 1;
    handles.stP1 = 1;
    handles.stStorage1 = nan(size(handles.stStorage1));
end
% end
guidata(hFigure,handles);

function schedulePULSE(hTimer,eventdata,hFigure)
%{
% announce when stimulus will go off
eventTime = datestr(eventdata.Data.time);
disp(['at ',eventTime,...
    ' stimulus pulse scheduled for NSP time ',...
    num2str(hTimer.UserData),...
    ' in ',num2str(hTimer.StartDelay),' s'])
%}

function finishPULSE(hTimer,eventdata,hFigure)
%{
% announce that this stim has completed or been aborted.
eventTime = datestr(eventdata.Data.time);
disp(['at ',eventTime,...
    ' stimulus scheduled for NSP time ',...
    num2str(hTimer.UserData),...
    ' has been completed or aborted.'])
%}

% ----------------------------------------------------------------------- %
% ----                                                                --- %
% ----                      Helper Functions                          --- %
% ----                                                                --- %
% ----------------------------------------------------------------------- %

function newBuffer = cycleBuffer(oldBuffer, newData)
N = length(newData);
if N >= length(oldBuffer)
    newBuffer = newData(end-length(oldBuffer)+1:end);
else
    newBuffer = [oldBuffer(N+1:end); newData];
end


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


% --- New Helpers ---


function plotData = plotLogical(logData)
% take in a logical array and output an array that will plot true as 1 
% and will not plot false
plotData = double(logData); 
plotData(~logData) = nan; 


function [storage1, p1, storage2, p2] = ...
    cycleStorage(storage1, p1, storage2, newData)
N = length(newData); 
if p1+N-1 > length(storage1)
    % storage 1 is now full 
    if N > length(storage2)
        warning('Data overloaded save buffer; some data may not be saved.')
        N = length(storage2);
        newData = newData(1:N);
    end
    p1 = 0; 
    storage2(1:N) = newData; 
    p2 = N+1;
else
    storage1(p1:(p1+N-1)) = newData;
    p1 = p1+N;
    p2 = [];
end


function [newBuffer, lastBuffer] = CombineAndCycle(oldBuffer, newData, N, M)
% ?? does this still work when N is longer than length oldBuffer ??
% M = length(newData); 
newBuffer = false(size(oldBuffer)); 
newBuffer(1:(end-N)) = oldBuffer((N+1):end);
lastBuffer = oldBuffer; 
if N < length(lastBuffer)
    lastBuffer = lastBuffer(1:N);
end
newBuffer((end-M+1):end) = newData((end-M+1):end);


function newBuffer = OverwriteAndCycle(oldBuffer, newData, N)
% N = # of points of data that is actually new 
if N <= length(newData)
    if N >= length(oldBuffer)
        newBuffer = newData(end-length(oldBuffer)+1:end);
    else
        % buffer AND overwrite 
        L = length(newData)-N; % length to overwrite
        newBuffer = [oldBuffer((N+1):(end-L)); newData];
    end
else
    % there is no data to overwrite; in fact, there is not enough new data
    % nan-pad newData to length N and cycle buffer 
    newData = [newData; nan(N-length(newData),1)];
    newBuffer = cycleBuffer(oldBuffer, newData);
end


function [newFiltBuffer, filterFinalCond] = FilterAndCycle(...
    oldFiltBuffer, newUnfilt, filtobj, filterInitCond)
[newFilt,filterFinalCond] = filter(filtobj,1,newUnfilt,filterInitCond);
newFiltBuffer = cycleBuffer(oldFiltBuffer, newFilt);


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


% --- Executes on selection change in pop_RedStim.
function pop_RedStim_Callback(hObject, eventdata, handles)
% hObject    handle to pop_RedStim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_RedStim contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_RedStim
contents = cellstr(get(hObject,'String'));
RedStim = contents{get(hObject,'Value')};
if strcmpi(RedStim, 'Peak')
    RedStim = 1;
elseif strcmpi(RedStim, 'Trough')
    RedStim = 2;
else
    RedStim = 0;
end
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
if strcmpi(YellowStim, 'Peak')
    YellowStim = 1;
elseif strcmpi(YellowStim, 'Trough')
    YellowStim = 2;
else
    YellowStim = 0;
end
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
if strcmpi(GreenStim, 'Peak')
    GreenStim = 1;
elseif strcmpi(GreenStim, 'Trough')
    GreenStim = 2;
else
    GreenStim = 0;
end
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
if strcmpi(StopStim, 'Peak')
    StopStim = 1;
elseif strcmpi(StopStim, 'Trough')
    StopStim = 2;
else
    StopStim = 0;
end
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
handles.FilterSetUp = true;

%stop(handles.timer)
StopMainLoop(hObject,eventdata,handles)
guidata(hObject, handles)
%start(handles.timer) % restart timer and plots 
StartMainLoop(hObject,eventdata,handles)



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

n = str2double(get(handles.txt_AR,'String'));
N = str2double(get(handles.txt_PDSwin,'String'));
PDSwin = ceil(N*handles.fSample); handles.PDSwin1 = PDSwin;

try

    handles = helperGUIv0_pushAR(handles, PDSwin, n);
    
    %stop(handles.timer)
    StopMainLoop(hObject,eventdata,handles)
    pause(.01)
    guidata(hObject, handles)
    %start(handles.timer) % restart timer and plots
    StartMainLoop(hObject,eventdata,handles)
    pause(.001)

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

% if there is a stimulator, stop it 
if handles.StimActive
    try
        handles.stimulator.stop();
        handles.stimulator.disconnect;
        pause(.1)
    catch ME0
        getReport(ME0)
        keyboard
        % should not get here because StimActive should be false
    end
end

if get(hObject, 'Value') == 1
    % start stimulus 

    try

    handles.stimMaxFreq = eval(get(handles.txt_MaxStimFreq, 'String'));

    amp1 = eval(handles.txt_amp1.String); 
    amp2 = eval(handles.txt_amp2.String); 
    width1 = eval(handles.txt_width1.String); 
    width2 = eval(handles.txt_width2.String); 
    interphase = eval(handles.txt_interphase.String); 
    frequency = 100; 
    pulses = 1;

    channel = nan(1,5); p = 1;
    for pop_ = [handles.pop_channel1, ...
                handles.pop_channel2, ...
                handles.pop_channel3, ...
                handles.pop_channel4, ...
                handles.pop_channel5]
        chanind = pop_.Value;
        if chanind <= size(pop_.String,1)
            channel(p) = str2double(pop_.String(chanind,:));
        end
        p = p+1;
    end
    channel1 = channel(1)
    channel2 = channel(2:end)

    handles.stimulator = defineSTIM4(channel1, channel2, amp1, amp2, ...
        width1, width2, interphase, frequency, pulses);

    handles.StimActive = true;
    set(hObject, 'String', 'Stim On'); 

    catch ME
        hObject.Value = 0;
        handles.StimActive = false;
        getReport(ME)
        errordlg(ME.message, 'Stim Setup Issue');
    end

else
    % stop stimulus 
    handles.StimActive = false;
    set(hObject, 'String', 'Stim Off');
end

guidata(hObject, handles)


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


% --- Executes on selection change in pop_elecgrid.
function pop_elecgrid_Callback(hObject, eventdata, handles)
% hObject    handle to pop_elecgrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_elecgrid contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_elecgrid
contents = cellstr(get(hObject,'String'));
sel = contents{get(hObject,'Value')};
handles.showElecGrid = ~strcmp(sel, 'None');
try
if handles.showElecGrid
    if strcmp(sel, 'PAC')
        % ***** TO DO: phase amplitude coupling: handles.elecGridFunc = ... 
    else
        % __ band power
        if strcmp(sel, 'Selected Band Power')
            if ~handles.FilterSetUp
                error('Filter must be set for this selection.')
            end
            fbnd = [handles.locutoff, handles.hicutoff];
        elseif strcmp(sel, 'Beta Power')
            fbnd = [13, 30]; % Hz
        elseif strcmp(sel, 'Gamma Power')
            fbnd = [50, 200]; % Hz 
        end
        handles.elecGridFunc = @(data, fs) bandpower(data, fs, fbnd);
    end
end
catch ME4
    getReport(ME4);
    errordlg(ME4.message, 'Electrode Grid Selection Issue');
    handles.showElecGrid = false;
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



function txt_ArtStart_Callback(hObject, eventdata, handles)
% hObject    handle to txt_ArtStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_ArtStart as text
%        str2double(get(hObject,'String')) returns contents of txt_ArtStart as a double
handles.ArtifactStartBefore = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txt_ArtStart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_ArtStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_ArtDur_Callback(hObject, eventdata, handles)
% hObject    handle to txt_ArtDur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_ArtDur as text
%        str2double(get(hObject,'String')) returns contents of txt_ArtDur as a double
handles.ArtifactDuration = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txt_ArtDur_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_ArtDur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_StimulatorLagTime_Callback(hObject, eventdata, handles)
% hObject    handle to txt_StimulatorLagTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_StimulatorLagTime as text
%        str2double(get(hObject,'String')) returns contents of txt_StimulatorLagTime as a double
handles.StimulatorLagTime = str2double(get(hObject,'String'));
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function txt_StimulatorLagTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_StimulatorLagTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_PowerThresh_Callback(hObject, eventdata, handles)
% hObject    handle to txt_PowerThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_PowerThresh as text
%        str2double(get(hObject,'String')) returns contents of txt_PowerThresh as a double
handles.bpthresh = str2double(get(hObject,'String'));
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function txt_PowerThresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_PowerThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_ARlearnrate_Callback(hObject, eventdata, handles)
% hObject    handle to txt_ARlearnrate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_ARlearnrate as text
%        str2double(get(hObject,'String')) returns contents of txt_ARlearnrate as a double
handles.ARlearnrate = str2double(get(hObject,'String'));
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function txt_ARlearnrate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_ARlearnrate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
