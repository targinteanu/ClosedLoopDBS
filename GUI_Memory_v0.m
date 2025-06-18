function varargout = GUI_Memory_v0(varargin)
% GUI_MEMORY_V0 MATLAB code for GUI_Memory_v0.fig
%      GUI_MEMORY_V0, by itself, creates a new GUI_MEMORY_V0 or raises the existing
%      singleton*.
%
%      H = GUI_MEMORY_V0 returns the handle to a new GUI_MEMORY_V0 or the handle to
%      the existing singleton*.
%
%      GUI_MEMORY_V0('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_MEMORY_V0.M with the given input arguments.
%
%      GUI_MEMORY_V0('Property','Value',...) creates a new GUI_MEMORY_V0 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_Memory_v0_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_Memory_v0_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_Memory_v0

% Last Modified by GUIDE v2.5 14-Oct-2024 14:59:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_Memory_v0_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_Memory_v0_OutputFcn, ...
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


% --- Executes just before GUI_Memory_v0 is made visible.
function GUI_Memory_v0_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_Memory_v0 (see VARARGIN)

% Choose default command line output for GUI_Memory_v0
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


handles.cbmexStatus = false;

% start receiver serial communication from paradigm computer
handles.textSrl.String = 'attempting to start serial com here ...';
thisportname = FindMySerialPort();
noSerialSetup = isempty(thisportname);
if ~noSerialSetup
    receiverSerial = serialport(thisportname, 9600);
end
ud = struct('ReceivedData', '', ...
            'TrialNumber', -1, ...
            'StimOn', false, ...
            ...'ParadigmPhase', 'WAIT', ...
            'ParadigmPhase', 'HOLD', ...
            'ImageVisible', false);
receiverSerial.UserData = ud;
if ~noSerialSetup
configureCallback(receiverSerial,"terminator",...
    @(hsrl,evt)CharSerialCallbackReceiver_Memory_v0(hsrl,evt, ...
                    handles.textSrl, handles.ParadigmInfoTable));
end
handles.srl = receiverSerial; 

% initiate other vars ...
handles.StimActive = false;
handles.RunMainLoop = false;
handles.FilterSetUp = false;
handles.MdlSetUp    = false;
handles.showElecGrid = false;
handles.srlLastMsg  = ud.ReceivedData;
handles.stimLastTime = -inf;
handles.stimNewTime = -inf;
handles.artReplaceRemaining = [];
handles.ArtifactDuration = .04; % set artifact duration (seconds) 
handles.ArtifactStartBefore = .01; % artifact start uncertainty (seconds)
handles.stimind = -1;

emptyStorage = nan(100000,1);
handles.pkStorage1 = emptyStorage; handles.pkP1 = 1;
handles.pkStorage2 = emptyStorage; 
handles.trStorage1 = emptyStorage; handles.trP1 = 1;
handles.trStorage2 = emptyStorage; 
handles.stStorage1 = emptyStorage; handles.stP1 = 1;
ud.TimeStamp = nan;
handles.srlStorage1 = repmat(ud,[1000,1]);
handles.srlP1 = 1; 

svloc = ['Saved Data Memory',filesep,'Saved Data ',...
    datestr(datetime, 'yyyy-mm-dd HH.MM.SS')];
pause(1)
mkdir(svloc); 
%handles.SaveFileName = [svloc,filesep,'SaveFile'];
handles.SaveFileLoc = svloc;
handles.SaveFileN = 1;

handles.StimulatorLagTime = 0.03; 
handles.QueuedStim = timer(...
    'StartDelay', 10, ...
    'TimerFcn',   {@myPULSE, hObject}, ...
    'StopFcn',    {@finishPULSE, hObject}, ...
    'StartFcn',   {@schedulePULSE, hObject}, ...
    'UserData',   Inf);

% Update handles structure
guidata(hObject, handles);

clc
%clear all

% UIWAIT makes GUI_Memory_v0 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_Memory_v0_OutputFcn(hObject, eventdata, handles) 
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
% set channel popup meno to hold channels
set(handles.pop_channels,'String',handles.channelList);
set(handles.pop_channel1,'String',handles.channelList);
set(handles.pop_channel2,'String',handles.channelList);
% electrode grid 
axes(handles.ax_elecgrid)
ncol = 3;
nrow = 21;
img = nan(nrow, ncol);
[X,Y] = meshgrid(1:ncol, 1:nrow);
handles.elecGridImg = imagesc(img); colormap('parula'); colorbar; hold on;
for ch = 1:min(length(handles.channelList),63)
    chname = handles.channelList(ch); chname = num2str(chname);
    text(X(ch),Y(ch), chname, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'Color','#A2142F');
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

    % artifact removal (1)
    if handles.FilterSetUp
    if handles.MdlSetUp
    if handles.check_artifact.Value
        if numel(handles.artReplaceRemaining)
            try
            % continue replacing from last loop iter
            artLen = min(length(handles.artReplaceRemaining), length(newContinuousData));
            artReplaceRemaining = handles.artReplaceRemaining(1:artLen); 
            handles.artReplaceRemaining = handles.artReplaceRemaining((artLen+1):end);
            newContinuousData(1:artLen) = artReplaceRemaining;
            catch ME3
                getReport(ME3)
                % keyboard
                errordlg(ME3.message, 'Artifact Removal Issue');
                handles.check_artifact.Value = false;
                pause(.01);
            end
        end
    end
    end
    end

    handles.rawDataBuffer = cycleBuffer(handles.rawDataBuffer, newContinuousData);
    N = length(newContinuousData);
    t0 = handles.lastSampleProcTime; 
    handles.lastSampleProcTime = ...
        time + (length(newContinuousData)-1)/handles.fSample;
    T = handles.lastSampleProcTime - t0;
    if T < 0
        %warning('Reported sample time is negative.')
    end
    diffSampleProcTime = nan(size(newContinuousData)); diffSampleProcTime(end) = T;
    handles.diffSampleProcTime = cycleBuffer(handles.diffSampleProcTime, diffSampleProcTime);

    % artifact removal (2) 
    if handles.FilterSetUp
    if handles.MdlSetUp
    if handles.check_artifact.Value
        stimind = handles.stimind - N; % N samples have passed
        if stimind > 0
            try
            artInd = stimind;
            artStart = -ceil(handles.ArtifactStartBefore*handles.fSample);
            artEnd = ceil(handles.fSample*handles.ArtifactDuration) - artStart -1;
            artInd = (artStart:artEnd) + artInd ...
                + round(handles.StimulatorLagTime*handles.fSample);
            artInd = artInd(artInd > 0); % can't change the past
            artLen = length(artInd);
            artInd = artInd(artInd <= handles.bufferSize);
            artPastData = zeros(handles.PDSwin1,1);
            if ~isempty(artInd)
                artPastStart = artInd(1) - handles.PDSwin1;
                artPastStart = max(1, artPastStart);
                rawOffset = mean(handles.rawDataBuffer);
                artPastData1 = handles.rawDataBuffer(artPastStart:(artInd(1)-1),:);
                artPastData1 = artPastData1 - rawOffset;
                artPastN = size(artPastData1,1);
                artPastData((end-artPastN+1):end,:) = artPastData1;
                artReplace = myFastForecastAR(handles.Mdl, artPastData, artLen);
                artReplace = artReplace + rawOffset;
                handles.artReplaceRemaining = artReplace((length(artInd)+1):end); % continue replacing on next loop iter
                artReplace = artReplace(1:length(artInd));
                handles.rawDataBuffer(artInd) = artReplace;
            end
            catch ME3
                getReport(ME3)
                % keyboard
                errordlg(ME3.message, 'Artifact Removal Issue');
                handles.check_artifact.Value = false;
                pause(.01);
            end
        end
    end
    end
    end

    guidata(hObject,handles)

    % update serial log 
    % TO DO: there should be a better way to do this; serial callback
    % should trigger an event or listener that logs the info 
    ReceivedData = handles.srl.UserData.ReceivedData; 
    if ~strcmp(ReceivedData, handles.srlLastMsg)
        ud = handles.srl.UserData; 
        ud.TimeStamp = handles.lastSampleProcTime;
        if handles.srlP1 <= length(handles.srlStorage1)
            handles.srlStorage1(handles.srlP1) = ud;
            handles.srlP1 = handles.srlP1+1;
        else
            ud = struct('ReceivedData', '', ...
                'TrialNumber', -1, ...
                'StimOn', false, ...
                'ParadigmPhase', 'WAIT', ...
                'ImageVisible', false, ...
                'TimeStamp', nan);
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
            getReport(ME4)
            errordlg(ME4.message, 'Electrode Grid Issue');
            handles.showElecGrid = false;
            pause(.01);
        end
    end

    % update filtered data 
    if handles.FilterSetUp
        try
        [handles.filtDataBuffer, handles.filtCond] = FilterAndCycle(...
            handles.filtDataBuffer, newContinuousData, ...
            handles.BPF, handles.filtCond);
        set(handles.h_filtDataTrace,'YData',handles.filtDataBuffer)

        %{
        % ===== FILTER DEBUGGING =====
        % (This imposes unnecessary computational burden and should be
        % removed or commented-out when debugging is complete.)
        x1 =   filter(handles.BPF,1, handles.rawDataBuffer);
        x2 = filtfilt(handles.BPF,1, handles.rawDataBuffer);
        x1 = x1((handles.IndShiftFIR + 1):end);
        set(handles.h_fullfiltTrace, 'YData', x1);
        set(handles.h_filtfiltTrace, 'YData', x2);
        % ============================
        %}

        % update model-forecasted data 
        if handles.MdlSetUp
            try

            % use AR model to get some future data 
            StimulatorLagInd = round(handles.fSample*handles.StimulatorLagTime);
            dataPast = handles.filtDataBuffer; 
            dataPast = dataPast((end-handles.PDSwin1+1):end);
            dataFutu = myFastForecastAR(handles.Mdl, dataPast, handles.PDSwin1);
            dataFutu2 = dataFutu(1:handles.PDSwin2);
            if handles.check_polar.Value
            handles.predDataBuffer = OverwriteAndCycle(...
                handles.predDataBuffer, dataFutu, N);
            set(handles.h_predTrace,'YData',handles.predDataBuffer);
            end

            % find the time to next peak, trough and plot 
            bp = norm(dataPast,2)^2/numel(dataPast); % band power surrogate 
            M = handles.PDSwin1 - handles.IndShiftFIR;
            dataPk = false(M,1); % +1? 
            dataTr = dataPk; % dataSt = dataPk;  
            [t2,i2,phi_inst,f_inst] = ...
                blockPDS(dataPast,dataFutu2, handles.fSample, [0,pi], ...
                handles.TimeShiftFIR + handles.StimulatorLagTime, ...
                handles.locutoff, handles.hicutoff);
            t2 = t2 - handles.TimeShiftFIR - handles.StimulatorLagTime; 
            i2 = i2 - handles.IndShiftFIR;% - StimulatorLagInd;
            %t2 = max(t2,0); i2 = max(i2,1);
            t2peak = t2(1); t2trou = t2(2);
            i2peak = i2(1); i2trou = i2(2);
            if i2peak > 0
                dataPk(i2peak) = true;
            end
            if i2trou > 0
                dataTr(i2trou) = true;
            end
            [handles.peakDataBuffer, oldPeak] = CombineAndCycle(...
                handles.peakDataBuffer, dataPk, N, M-StimulatorLagInd); 
            [handles.trouDataBuffer, oldTrou] = CombineAndCycle(...
                handles.trouDataBuffer, dataTr, N, M-StimulatorLagInd);
            set(handles.h_peakTrace,'YData',0*plotLogical(handles.peakDataBuffer));
            set(handles.h_trouTrace,'YData',0*plotLogical(handles.trouDataBuffer));

            % time of stimulus 
            if handles.stimNewTime > 0
            stimtimerel = handles.stimNewTime - handles.lastSampleProcTime; 
                % rel to 0 on screen
                % lastSampleProcTime should be time 0 on the screen
            stimind = round(stimtimerel*handles.fSample); % ind rel to END of buffer
            stimind = stimind + handles.bufferSize; % ind rel to START of buffer 
            handles.stimind = stimind - N;
            handles.stimNewTime = -inf;
            if stimind > 0
                handles.stimDataBuffer(stimind) = true;
            end
            else
                handles.stimind = -1;
            end

            [handles.stimDataBuffer, oldStim] = CombineAndCycle(...
                handles.stimDataBuffer, [], N, 0);
            set(handles.h_stimTrace,'YData',0*plotLogical(handles.stimDataBuffer));

            % queue stimulus pulse, if applicable 
            % ***** TO DO: can this be moved elsewhere to avoid delays?  
            Stim2Q = false;
            if bp > 1000 % min band power cutoff; orig at 1000
            if handles.StimActive
                ParadigmPhase = handles.srl.UserData.ParadigmPhase;
                if ~strcmpi(ParadigmPhase,'WAIT')
                    try
                        StimMode = getfield(handles.StimMode, ParadigmPhase);
                    catch
                        warning(['ParadigmPhase ',ParadigmPhase,' unrecognized.'])
                        StimMode = 'None';
                    end
                    if strcmpi(StimMode,'Peak')
                        t2Q = t2peak;
                        Stim2Q = true;
                    end
                    if strcmpi(StimMode,'Trough')
                        t2Q = t2trou;
                        Stim2Q = true;
                    end
                end
            end
            end
            if Stim2Q && (t2Q >= 0)
                % possibly overwrite existing timer with new one
                t2Q = .001*floor(1000*t2Q); % round to nearest 1ms 
                if t2Q < (100/handles.locutoff + handles.TimeShiftFIR)
                    t2Qabs = t2Q + handles.lastSampleProcTime; % in NSP time "absolute"
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

            % plot sine wave 
            if handles.check_polar.Value
            tSin = (1:handles.PDSwin1)/handles.fSample; tSin = tSin';
            sinMag = handles.filtDataBuffer(end) / cos(phi_inst);
            sinData = sinMag*cos(2*pi*f_inst*tSin + phi_inst);
            handles.sineDataBuffer = OverwriteAndCycle(...
                handles.sineDataBuffer, sinData, N);
            set(handles.h_sineTrace,'YData',handles.sineDataBuffer);
            end

            % evaluate accuracy of above --> polar histogram
            if handles.check_polar.Value
            phiPk = handles.peakDataBuffer; 
            phiTr = handles.trouDataBuffer;
            phiPk = phiPk(1:length(handles.filtDataBuffer));
            phiTr = phiTr(1:length(handles.filtDataBuffer));
            phi = instPhaseFreq(handles.filtDataBuffer, handles.fSample);
            phiPk = phi(phiPk); phiTr = phi(phiTr);
            set(handles.h_peakPhase,'Data',phiPk);
            set(handles.h_trouPhase,'Data',phiTr);
            end

            % store peaks and troughs that have been buffered out
                % N samples of new continuous data have come in 
                % most recent time stamp is lastSampleProcTime 
            tOldBuffer = ((1-N):0)/handles.fSample + handles.lastSampleProcTime - ...
                handles.bufferSize/handles.fSample;
            tOldPeak = tOldBuffer(oldPeak); tOldTrou = tOldBuffer(oldTrou);
            [handles.pkStorage1, handles.pkP1,  handles.pkStorage2, p2] = ...
                cycleStorage(handles.pkStorage1, handles.pkP1, ...
                             handles.pkStorage2, tOldPeak);
            if ~handles.pkP1
                % storage full; save 
                PeakTime = handles.pkStorage1;
                svfn = [handles.SaveFileLoc,filesep,'SaveFile',num2str(handles.SaveFileN),'.mat'];
                disp(['Saving Peaks to ',svfn])
                save(svfn,'PeakTime');
                handles.SaveFileN = handles.SaveFileN + 1;
                handles.pkStorage1 = handles.pkStorage2; 
                handles.pkP1 = p2; 
                handles.pkStorage2 = nan(size(handles.pkStorage2));
            end
            [handles.trStorage1, handles.trP1,  handles.trStorage2, p2] = ...
                cycleStorage(handles.trStorage1, handles.trP1, ...
                             handles.trStorage2, tOldTrou);
            if ~handles.trP1
                % storage full; save 
                TroughTime = handles.trStorage1;
                svfn = [handles.SaveFileLoc,filesep,'SaveFile',num2str(handles.SaveFileN),'.mat'];
                disp(['Saving Troughs to ',svfn])
                save(svfn,'TroughTime');
                handles.SaveFileN = handles.SaveFileN + 1;
                handles.trStorage1 = handles.trStorage2; 
                handles.trP1 = p2; 
                handles.trStorage2 = nan(size(handles.trStorage2));
            end

            catch ME2
                getReport(ME2)
                errordlg(ME2.message, 'Model Prediction Issue');
                handles.MdlSetUp = false;
                pause(.01);
            end
        end

        catch ME1
            getReport(ME1)
            errordlg(ME1.message, 'Filtering Issue');
            handles.FilterSetUp = false;
            pause(.01);
        end
    end

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

    % Check which channel is selected and get some data to plot
    handles.channelIndex = get(handles.pop_channels,'Value');

    % get data from Central
    [events, time, continuousData] = cbmex('trialdata',1);
    
    % Check to make sure continuous sampling is enabled on at least one
    % channel
    if isempty(continuousData)
        % wait some time and try again before throwing error
        pause(.5)
        [events, time, continuousData] = cbmex('trialdata',1);
            if isempty(continuousData)
            errordlg(['Continuous acquisition not enabled.', ...
                'Select a sampling rate in Hardware Configuration'], ...
                'No Continuous Data')
            return
        end
    end
    newContinuousData = continuousData{handles.channelIndex,3};
    handles.fSample = continuousData{handles.channelIndex,2};

    % check units 
    config = cbmex('config', handles.channelIndex);
    unitname_is = lower(config{11,1});
    if contains(unitname_is, 'unit')
        unitname = config{11,2};
    else
        unitname_is = contains(lower(config(:,1)), 'unit');
        unitname_is = find(unitname_is);
        unitname_is = unitname_is(1); 
        unitname = config{unitname_is,2};
    end
    
    % Now that we know the sampling rate of the selected channel,
    % Create raw data buffer of zeros of the correct length
    handles.bufferSize = str2double(get(handles.txt_display,'String')) * handles.fSample;
    handles.rawDataBuffer = zeros(handles.bufferSize,1);

    if handles.FilterSetUp
        if handles.bufferSize > handles.IndShiftFIR
            handles.filtDataBuffer = ...
                zeros(handles.bufferSize - handles.IndShiftFIR, 1);
        else
            errordlg({'Data buffer is not long enough to be filtered.',...
                      'Try increasing display window or altering filter.'},...
                     'Filtering Issue');
            handles.FilterSetUp = false;
        end

        if handles.MdlSetUp
            sz = size(handles.filtDataBuffer); 
            sz(1) = sz(1) + handles.PDSwin1;
            handles.predDataBuffer = nan(sz);
            handles.peakDataBuffer = false(sz);
            handles.trouDataBuffer = false(sz);
            handles.sineDataBuffer = nan(sz);
            handles.stimDataBuffer = false(sz);
        end
    end

    % keep track of the proc time of the most recent data point.  This will
    % help if you want to match spike times with points in the buffer.
    % 'time' is the time at the first data point of the new chunk of
    % continuous data in seconds.
    handles.lastSampleProcTime = ...
        time + (length(newContinuousData)-1)/handles.fSample;
    handles.diffSampleProcTime = nan(handles.bufferSize,1);

    handles.rawDataBuffer = cycleBuffer(handles.rawDataBuffer, newContinuousData);
    xValues = linspace(-handles.bufferSize/handles.fSample,0,handles.bufferSize);
    axes(handles.ax_raw);
    handles.h_rawDataTrace = plot(xValues,handles.rawDataBuffer);
    grid on; title('Raw Channel Data'); xlabel('time (s)'); ylabel(unitname);
    common_xlim = xlim;

    axes(handles.ax_timing); 
    handles.h_timingTrace = stem(xValues,handles.diffSampleProcTime);
    grid on; title('Update Duration'); xlabel('time (s)'); ylabel('s');

    if handles.FilterSetUp
        try
        filtInitCond = zeros(handles.FilterOrder,1);
        [handles.filtDataBuffer, handles.filtCond] = FilterAndCycle(...
            handles.filtDataBuffer, newContinuousData, ...
            handles.BPF, filtInitCond);
        xValues2 = xValues(1:(end-handles.IndShiftFIR));
        common_xdiff = diff(common_xlim); 
        ext_xdiff = common_xdiff * handles.ax_filt.InnerPosition(3) / ...
            handles.ax_raw.InnerPosition(3); 
        ext_xlim = [0, ext_xdiff] + common_xlim(1); % align left 
        axes(handles.ax_filt); hold off; 
        handles.h_filtDataTrace = plot(xValues2, handles.filtDataBuffer); 
        grid on; hold on; 
        title('Filtered & Predicted Data'); xlabel('time (s)'); ylabel(unitname);
        xlim(ext_xlim);

        %{
        % ===== FILTER DEBUGGING =====
        % (This imposes unnecessary computational burden and should be
        % removed or commented-out when debugging is complete.)
        x1 =   filter(handles.BPF,1, handles.rawDataBuffer);
        x2 = filtfilt(handles.BPF,1, handles.rawDataBuffer);
        x1 = x1((handles.IndShiftFIR + 1):end);
        handles.h_fullfiltTrace = plot(xValues2, x1, ':', 'LineWidth',1);
        handles.h_filtfiltTrace = plot(xValues,  x2, ':', 'LineWidth',1);
        % ============================
        %}

        if handles.MdlSetUp
            xValues3 = (-handles.bufferSize) : ...
                       (handles.PDSwin1 - handles.IndShiftFIR - 1);
            xValues3 = xValues3/handles.fSample;
            handles.h_peakTrace = plot(xValues3,0*plotLogical(handles.peakDataBuffer), ...
                '^', 'Color',"#EDB120"); 
            handles.h_trouTrace = plot(xValues3,0*plotLogical(handles.trouDataBuffer), ...
                'v', 'Color',"#EDB120"); 
            handles.h_stimTrace = plot(xValues3,0*plotLogical(handles.stimDataBuffer), ...
                '*', 'Color','r'); 
            handles.h_predTrace = plot(xValues3,handles.predDataBuffer, ':');
            handles.h_sineTrace = plot(xValues3,handles.sineDataBuffer,'--');

            bedge = (-1:2:35)*pi/18; 
            axes(handles.ax_polar); hold off;
            handles.h_peakPhase = polarhistogram(nan, bedge);
            hold on; 
            handles.h_trouPhase = polarhistogram(nan, bedge);
            title('Actual Phase')
        end

        catch ME1
            getReport(ME1)
            errordlg(ME1.message, 'Filtering Issue');
            handles.FilterSetUp = false;
            pause(.01);
        end
    end

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
stimstatus = stimulator.getSequenceStatus(); 
if stimstatus == 2
    warning('Stimulator is already playing.')
end
stimtime1 = cbmex('time');
stimulator.play(1);
stimtime2 = cbmex('time');
dstimtime = stimtime2 - stimtime1; 
% disp time of pulse using eventdata
eventTime = datestr(eventdata.Data.time);
stimtime = .5*(stimtime1 + stimtime2);
stimschedtime = hTimer.UserData; 
disp(['Stimulus pulsed at ',eventTime,' within ',num2str(dstimtime),'s, '...
      num2str(stimtime - stimschedtime),'s late'])
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

% details from data 
srate = handles.fSample;
nyq            = srate*0.5;  % Nyquist frequency

% filtering bound rules 
minfac         = 2;    % this many (lo)cutoff-freq cycles in filter
min_filtorder  = 15;   % minimum filter length

% access desired parameters
hicutoff = str2double(handles.txt_hico.String); 
locutoff = str2double(handles.txt_loco.String); 

% check that specs are proper
if locutoff>0 & hicutoff > 0 & locutoff > hicutoff
    errordlg('locutoff > hicutoff ???', 'filter spec');
    return
end
if locutoff < 0 | hicutoff < 0
    errordlg('locutoff | hicutoff < 0 ???', 'filter spec');
    return
end
if locutoff>nyq
    errordlg('Low cutoff frequency cannot be > srate/2', 'filter spec');
    return
end
if hicutoff>nyq
    errordlg('High cutoff frequency cannot be > srate/2', 'filter spec');
    return
end
handles.hicutoff = hicutoff; handles.locutoff = locutoff; 

% filter order 
if locutoff>0
    filtorder = minfac*fix(srate/locutoff);
elseif hicutoff>0
    filtorder = minfac*fix(srate/hicutoff);
end
if filtorder < min_filtorder
    filtorder = min_filtorder;
end
handles.FilterOrder = filtorder;
handles.TimeShiftFIR = filtorder/(2*srate); % seconds 
handles.IndShiftFIR = ceil(filtorder/2); % samples ???

% build filter 
% usage i.e.: 
% >> filteredSignal = filter(filtwts, 1, unfilteredSignal) 
% -- OR --
% >> filteredSignal = filtfilt(filtwts, 1, unfilteredSignal)
filtwts = fir1(filtorder, [locutoff, hicutoff]./(srate/2));
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
handles.PDSwin2 = ceil(.1*PDSwin); 

try

    % catch mistakes 
    if PDSwin > handles.bufferSize
        error('Phase estimation window cannot be larger than display window.')
    end
    if handles.PDSwin1 <= handles.IndShiftFIR
        error('Forecast length does not overcome filter delay.')
    end

    y = handles.filtDataBuffer; 
    L = min(length(y), 3*PDSwin) - 1;
    y = y((end-L):end);
    y = iddata(y,[],1/handles.fSample);
    ARmdl = ar(y,n,'yw');
    
    handles.Mdl = ARmdl; 
    handles.MdlSetUp = true;
    
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

    chan1ind = handles.pop_channel1.Value; 
    chan2ind = handles.pop_channel2.Value; 
    channel1 = str2double(handles.pop_channel1.String(chan1ind,:)) 
    channel2 = str2double(handles.pop_channel2.String(chan2ind,:))
    %channel1 = chan1ind; 
    %channel2 = chan2ind;

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
        % ***** TO DO: ...
    else
        % ___ band power 
        if strcmp(sel, 'Selected Band Power')
            if ~handles.FilterSetUp
                error('Filter must be set for this selection.')
            end
            fbnd = [handles.locutoff, handles.hicutoff];
        elseif strcmp(sel, 'Theta Power')
            fbnd = [4, 9]; % Hz
        elseif strcmp(sel, 'Gamma Power')
            fbnd = [50, 200]; % Hz 
        end
        handles.elecGridFunc = @(data, fs) bandpower(data, fs, fbnd);
    end
end
catch ME4 
    getReport(ME4)
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
