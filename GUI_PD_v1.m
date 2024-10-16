function varargout = GUI_PD_v1(varargin)
% GUI_PD_V1 MATLAB code for GUI_PD_v1.fig
%      GUI_PD_V1, by itself, creates a new GUI_PD_V1 or raises the existing
%      singleton*.
%
%      H = GUI_PD_V1 returns the handle to a new GUI_PD_V1 or the handle to
%      the existing singleton*.
%
%      GUI_PD_V1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_PD_V1.M with the given input arguments.
%
%      GUI_PD_V1('Property','Value',...) creates a new GUI_PD_V1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_PD_v1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_PD_v1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_PD_v1

% Last Modified by GUIDE v2.5 16-Oct-2024 01:07:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_PD_v1_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_PD_v1_OutputFcn, ...
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


% --- Executes just before GUI_PD_v1 is made visible.
function GUI_PD_v1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_PD_v1 (see VARARGIN)

% Choose default command line output for GUI_PD_v1
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


handles.DAQstatus = false;

% start receiver serial communication from paradigm computer
handles.textSrl.String = 'attempting to start serial com here ...';
thisportname = FindMySerialPort();
noSerialSetup = isempty(thisportname);
if ~noSerialSetup
    receiverSerial = serialport(thisportname, 9600);
end
ud = struct('ReceivedData', '', ...
            'ParadigmPhase', 'Stopped');
receiverSerial.UserData = ud;
if ~noSerialSetup
configureCallback(receiverSerial,"terminator",...
    @(hsrl,evt)CharSerialCallbackReceiver_PD_v0(hsrl,evt, ...
                    handles.textSrl, handles.txt_Status)); 
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
handles.stimNewTime  = -inf;
handles.bufferSize     = 10; 
handles.bufferSizeGrid = 10;
handles.channelIndex = [];
handles.PhaseOfInterest = [0, pi];

ud.TimeStamp = nan;
handles.udBlank = ud;
handles.srlStorage1 = repmat(ud,[1000,1]);
handles.srlP1 = 1; 

svloc = ['Saved Data PD',filesep,'Saved Data ',...
    datestr(datetime, 'yyyy-mm-dd HH.MM.SS')];
pause(1)
mkdir(svloc); 
handles.SaveFileLoc = svloc;
handles.SaveFileName = [svloc,filesep,'SaveFile'];
handles.SaveFileN = 1;

% start parallel pool(s) 
handles.pool = gcp('nocreate');
if isempty(handles.pool)
    handles.pool = parpool;
end

% Set up a data queue(s)
handles.userQueue = parallel.pool.PollableDataQueue;
handles.dataQueue = parallel.pool.PollableDataQueue;
handles.stimQueue = parallel.pool.PollableDataQueue;
%handles.modlQueue = parallel.pool.PollableDataQueue;

% remove these fields when sending handles over the userQueue
handles.rmfieldList = {...
    'userQueue', 'dataQueue', 'stimQueue', 'modlQueue', ...
    'f_PhaseDetect'};

% Update handles structure
guidata(hObject, handles);

clc
%clear all

% UIWAIT makes GUI_PD_v1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_PD_v1_OutputFcn(hObject, eventdata, handles) 
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

connect_cbmex();
handles.time0 = datetime - seconds(cbmex('time'));
disconnect_cbmex();

handles.f_PhaseDetect = parfeval(handles.pool, @bg_PhaseDetect, 1, ...
    handles.userQueue, handles.dataQueue, handles.stimQueue, ...
    @InitializeRecording_cbmex, @disconnect_cbmex, @getNewRawData_cbmex, []);

handles.DAQstatus = true;

% Acquire some data to get channel information. Determine which channels
% are enabled
send(handles.userQueue, rmfield(handles, handles.rmfieldList));
[dataRecd, handles.SaveFileN, ~,~, ~,~,~,~, rawInfo] = ...
    pollDataQueue_PhaseDetect_v1(handles.dataQueue, ...
    handles.SaveFileName, handles.SaveFileN, handles.time0, 10);
if ~dataRecd
    handles.DAQstatus = false;
    send(handles.userQueue, rmfield(handles, handles.rmfieldList));
    warning('Data aquisition timed out.')
else

% set channel popup menu to hold channels
handles.fSamples = cellfun(@(ch) ch.SampleRate, rawInfo);
handles.channelIDlist = cellfun(@(ch) ch.IDnumber, rawInfo);
handles.channelList = cellfun(@(ch) [num2str(ch.IDnumber),': ',ch.Name], ...
    rawInfo, 'UniformOutput',false);
chL = [handles.channelList, 'None'];
set(handles.pop_channels, 'String', handles.channelList);
for pop_ = [handles.pop_channel1, ...
            handles.pop_channel2, ...
            handles.pop_channel3, ...
            handles.pop_channel4, ...
            handles.pop_channel5]
    set(pop_, 'String', chL);
end

% electrode grid 
gridmaxval = eval(handles.txt_gridmax.String);
gridminval = eval(handles.txt_gridmin.String);
axes(handles.ax_elecgrid)
ncol = 3; % # of columns
nrow = 21;
img = nan(nrow, ncol);
[X,Y] = meshgrid(1:ncol, 1:nrow);
handles.elecGridImg = imagesc(img, [gridminval, gridmaxval]); 
colormap('parula'); colorbar; hold on;
text(X(:),Y(:), handles.channelList, ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', ...
    'FontWeight', 'bold', ...
    'Color',[.8 0 0]);

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
handles.DAQstatus = false;
send(handles.userQueue, rmfield(handles, handles.rmfieldList));
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
    if ~handles.DAQstatus
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
% stop stim
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
cbmex('close')
handles.DAQstatus = false;
send(handles.userQueue, rmfield(handles, handles.rmfieldList));
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

try
% save stored data 
SerialLog = handles.srlStorage1;
svfn = [handles.SaveFileName,num2str(handles.SaveFileN),'.mat'];
disp(['Saving Serial to ',svfn])
save(svfn,'SerialLog');
disp('Saving all data...')
ConsolidateSavedData_v1(handles.SaveFileLoc)
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
function updateDisplay(hObject, eventdata)

handles = guidata(hObject);
while handles.RunMainLoop

try

    % timing 
    pause(.1); % s between displays 
    tNow = datetime;
    timeDisp2 = handles.timeDisp0 + toc(handles.timeDisp1);
    handles.timeDispBuff = bufferData(handles.timeDispBuff, timeDisp2);
    
    % ensure connection with hardware 
    handles = guidata(hObject);
    if ~handles.DAQstatus
        %stop(handles.timer)
        StopMainLoop(hObject,eventdata,handles)
    end
    
    % get data from Central
    [dataRecd, handles.SaveFileN, timeBuff, forBuff, ...
    tPltRng, rawPlt, fltPlt, forPlt, ...
    rawD1, rawD4, fltD1, fltD4, forD1, forD4] = ...
    pollDataQueue_PhaseDetect_v1(handles.dataQueue, ...
        handles.SaveFileName, handles.SaveFileN, handles.time0, 10);
    if ~dataRecd
        error('Data aquisition timed out.')
    end
    lastSampleProcTime = timeBuff(end);

    % update serial log 
    % TO DO: there should be a better way to do this; serial callback
    % should trigger an event or listener that logs the info 
    % Should this be in background instead?
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
            svfn = [handles.SaveFileName,num2str(handles.SaveFileN),'.mat'];
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
            for ch = 1:63
                x = rawD4{ch}(:,2); L = handles.bufferSizeGrid(ch);
                if height(x) > L
                    x = x((end-L+1):end, :);
                end
                if height(x) < L
                    warning(['Channel ',num2str(ch),' Electrode Grid buffer is not full length!'])
                end
                fSample_ch = handles.fSamples(ch);
                elecimg(ch) = handles.elecGridFunc(x, fSample_ch);
            end
            handles.elecGridImg.CData = elecimg;
        catch ME4 
            getReport(ME4);
            errordlg(ME4.message, 'Electrode Grid Issue');
            handles.showElecGrid = false;
            pause(.01);
        end
    end

    % update filtered data plot
    if handles.FilterSetUp
        try
        set(handles.h_filtDataTrace,'YData',fltPlt.Variables);
        set(handles.h_filtDataTrace,'XData',fltPlt.Time - tNow);

        % update model-forecasted data plot
        if handles.MdlSetUp
            try
            if handles.check_polar.Value
                set(handles.h_predTrace,'YData',forPlt.Variables);
                set(handles.h_predTrace,'XData',forPlt.Time - tNow);
            end
            tPk = forBuff(:,1); tTr = forBuff(:,2);
            set(handles.h_peakTrace,'YData',zeros(size(tPk)));
            set(handles.h_peakTrace,'XData',handles.time0 + seconds(tPk) - tNow);
            set(handles.h_trouTrace,'YData',zeros(size(tTr)));
            set(handles.h_trouTrace,'XData',handles.time0 + seconds(tTr) - tNow);

            % time of stimulus 
            if handles.stimNewTime > 0
            stimtimerel = handles.stimNewTime - handles.lastSampleProcTime; 
                % rel to 0 on screen
                % lastSampleProcTime should be time 0 on the screen
            stimind = round(stimtimerel*handles.fSample); % ind rel to END of buffer
            stimind = stimind + handles.bufferSize; % ind rel to START of buffer 
            handles.stimNewTime = -inf;
            if stimind > 0
                handles.stimDataBuffer(stimind) = true;

                % artifact removal 
                if handles.check_artifact.Value
                    try
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

            set(handles.h_stimTrace,'YData',0*plotLogical(handles.stimDataBuffer));

            % queue stimulus pulse, if applicable 
            % ***** TO DO: can this be moved elsewhere to avoid delays?  
            if handles.StimActive
                ParadigmPhase = handles.srl.UserData.ParadigmPhase;
                if ~strcmpi(ParadigmPhase,'WAIT')
                    if strcmpi(ParadigmPhase, 'Started') || strcmp(ParadigmPhase, 'gray')
                        % Started, gray, and red should all be the same.
                        ParadigmPhase = 'red';
                    end
                    try
                        StimMode = getfield(handles.StimMode, ParadigmPhase);
                    catch
                        warning(['ParadigmPhase ',ParadigmPhase,' unrecognized.'])
                        StimMode = 'None';
                    end
                end
            end

            % plot sine wave 
            if handles.check_polar.Value
                % set(handles.h_sineTrace,'YData', ...
                % set(handles.h_sineTrace,'XData', ...
            end

            % evaluate accuracy of above --> polar histogram
            if handles.check_polar.Value

                % row indexes of peak events 
                rowPk = nan(size(tPk)); 
                for r = 1:height(tPk)
                    % find time of current event relative to time now
                    tPk_r = handles.time0 + seconds(tPk(r)) - tNow ;
                    if (tPk_r >= fltPlt.Time(1)) && (tPk_r <= fltPlt.Time(end))
                        % current event is in time range shown on screen,
                        % so let row index be the nearest 
                        [~,rowPk(r)] = min(abs( tPk_r - fltPlt.Time ));
                    end
                end
                rowPk = rowPk(~isnan(rowPk));

                % row indexes of trough events 
                rowTr = nan(size(tTr)); 
                for r = 1:height(tTr)
                    % find time of current event relative to time now
                    tTr_r = handles.time0 + seconds(tTr(r)) - tNow ;
                    if (tTr_r >= fltPlt.Time(1)) && (tTr_r <= fltPlt.Time(end))
                        % current event is in time range shown on screen,
                        % so let row index be the nearest 
                        [~,rowTr(r)] = min(abs( tTr_r - fltPlt.Time ));
                    end
                end
                rowTr = rowTr(~isnan(rowTr));

                % calc phase and histogram
                phi = instPhaseFreq(fltPlt.Variables, handles.fSample);
                phiPk = phi(rowPk); phiTr = phi(rowTr);
                set(handles.h_peakPhase,'Data',phiPk);
                set(handles.h_trouPhase,'Data',phiTr);
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

    % update raw data and timing plots
    set(handles.h_rawDataTrace,'YData',rawPlt.Variables);
    set(handles.h_rawDataTrace,'XData',rawPlt.Time - tNow);
    set(handles.h_timingTrace,'YData',[nan; diff(timeBuff)]);
    set(handles.h_timingTrace,'XData', ...
        handles.time0 + seconds(timeBuff) - tNow );
    set(handles.h_timeDispTrace,'YData',[nan; diff(handles.timeDispBuff)]);
    set(handles.h_timeDispTrace,'XData', ...
        handles.time0 + seconds(handles.timeDispBuff) - tNow );

    guidata(hObject,handles)

catch ME 
    getReport(ME)
    StopMainLoop(hObject,eventdata,handles)
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

    tNow = datetime;

    % Check which channel is selected and get some data to plot
    handles.channelIndex = get(handles.pop_channels,'Value'); 
    % Now we know the sampling rate of the selected channel
    handles.fSample = handles.fSamples(handles.channelIndex);
    handles.bufferSize = str2double(get(handles.txt_display,'String')) * handles.fSample;
    handles.bufferSizeGrid = str2double(get(handles.txt_griddur,'String')) * handles.fSamples;

    % get data from Central
    send(handles.userQueue, rmfield(handles, handles.rmfieldList));
    [dataRecd, handles.SaveFileN, timeBuff, forBuff, ...
    tPltRng, rawPlt, fltPlt, forPlt, ...
    rawD1, rawD4, fltD1, fltD4, forD1, forD4] = ...
    pollDataQueue_PhaseDetect_v1(handles.dataQueue, ...
        handles.SaveFileName, handles.SaveFileN, handles.time0, 10);
    if ~dataRecd
        error('Data aquisition timed out.')
    end

    % keep track of the display time 
    handles.timeDisp1 = tic; handles.timeDisp0 = timeBuff(end, :);
    handles.timeDispBuff = nan(size(timeBuff));

    % initiate raw data plot
    axes(handles.ax_raw);
    handles.h_rawDataTrace = plot(rawPlt.Time - tNow, rawPlt.Variables);
    grid on; title('Raw Channel Data'); 
    xlabel('time'); ylabel(rawPlt.Properties.VariableNames{1});
    if sum(~isnat(tPltRng))
        common_xlim = tPltRng; 
        xlim(common_xlim);
    else
        common_xlim = xlim();
    end

    % initiate timing stem plot
    axes(handles.ax_timing); 
    handles.h_timingTrace = ...
        stem(handles.time0 + seconds(timeBuff) - tNow, ...
        [nan; diff(timeBuff)], 's');
    handles.h_timeDispTrace = ...
        stem(handles.time0 + seconds(handles.timeDispBuff) - tNow, ...
        [nan; diff(handles.timeDispBuff)], 's');
    grid on; title('Update Duration'); xlabel('time (s)'); ylabel('s');
    xlim(common_xlim);

    % initiate filtered data plot
    if handles.FilterSetUp
        try        
        common_xdiff = diff(common_xlim); 
        ext_xdiff = common_xdiff * handles.ax_filt.InnerPosition(3) / ...
            handles.ax_raw.InnerPosition(3); 
        ext_xlim = [0, ext_xdiff] + common_xlim(1); % align left 
        axes(handles.ax_filt); hold off; 
        handles.h_filtDataTrace = plot(fltPlt.Time - tNow, fltPlt.Variables); 
        grid on; hold on; 
        title('Filtered & Predicted Data'); xlabel('time (s)'); ylabel(unitname);
        xlim(ext_xlim);

        % initiate prediction & peak/trough indicators overlayed on
        % filtered plot
        if handles.MdlSetUp
            tPk = forBuff(:,1); tTr = forBuff(:,2);
            handles.h_peakTrace = plot(handles.time0 + seconds(tPk) - tNow, zeros(size(tPk)), ...
                '^', 'Color',"#EDB120"); 
            handles.h_trouTrace = plot(handles.time0 + seconds(tTr) - tNow, zeros(size(tTr)), ...
                'v', 'Color',"#EDB120"); 
            handles.h_stimTrace = plot(xValues3,0*plotLogical(handles.stimDataBuffer), ...
                '*', 'Color','r'); 
            handles.h_predTrace = plot(forPlt.Time - tNow, forPlt.Variables, ':');
            %handles.h_sineTrace = plot(xValues3,handles.sineDataBuffer,'--');

            % initiate polar histogram 
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
            send(handles.userQueue, rmfield(handles, handles.rmfieldList));
            pause(.01);
        end
    end

    handles.RunMainLoop = true; 
    send(handles.userQueue, rmfield(handles, handles.rmfieldList));
    guidata(hObject,handles)
    updateDisplay(hObject,eventdata)
    
catch ME
    getReport(ME)
    handles.RunMainLoop = false; 
    send(handles.userQueue, rmfield(handles, handles.rmfieldList));
    guidata(hObject, handles);
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
    svfn = [handles.SaveFileName,num2str(handles.SaveFileN),'.mat'];
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


function [newBuffer, lastBuffer] = CombineAndCycle(oldBuffer, newData, N)
% ?? does this still work when N is longer than length oldBuffer ??
M = length(newData); 
newBuffer = false(size(oldBuffer)); 
newBuffer(1:(end-N)) = oldBuffer((N+1):end);
lastBuffer = oldBuffer; 
if N < length(lastBuffer)
    lastBuffer = lastBuffer(1:N);
end
newBuffer((end-M+1):end) = newData;


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


function dataForecast = MdlForecast(MdlObj, dataPast, k, fs)
% This needs to be made faster. Since fs is same as model fit, can do
% by direct multiplication instead of built-in? 
% Forecast the next <k> datapoints <dataForecast> using the model
% system <MdlObj> and the previous data <dataPast> sampled at constant rate
% <fs>. All data is in columns. 
dataPast = iddata(dataPast,[],1/fs); 
dataForecast = forecast(MdlObj,dataPast,k);
dataForecast = dataForecast.OutputData;


function setTgl(hObject, eventdata, handles, hTgl, newValue)
% set a toggle button to a desired Value and activate its callback if it is
% not currently at that value. 
curValue = hTgl.Value; 
hTgl.Value = newValue; guidata(hObject, handles);
if ~(curValue == newValue)
    hTgl.Callback(hTgl, eventdata);
    guidata(hObject, handles);
end


function [t2phi, i2phi, phi_inst, f_inst] = ...
    blockPDS(pastData, futureData, fs, phi, tmin, fmin, fmax)
% Determine the time (s) and # samples to next desired phase phi from a
% block of data sampled at a constant rate fs (Hz). Also return the current
% inst. phase phi_inst (rad) and frequency f_inst (Hz).
% Block data should include some length of pastData and
% (forecasted/predicted) futureData to minimize edge effects at the present
% timepoint, which is the last element of pastData. Data should be input as
% columns. 
% phi [desired] is in radians, i.e. phi=0 for peak, phi=pi for trough
% frequency will be clipped within range [fmin, fmax] (Hz) 

N = size(pastData,1); M = size(futureData,1);
blockData = [pastData; futureData];

[phi_block, f_block] = instPhaseFreq(blockData, fs);
phi_inst = phi_block(N,:);
f_block = max(f_block, fmin); 
f_block = min(f_block, fmax);
fwinlen = floor(.03*N); fwinlen = min(fwinlen, M);
fwin = N + ((-fwinlen):fwinlen);
f_inst = mean(f_block(fwin,:));
T=1/f_inst;

% time to next [desired] phi 
t2phi = zeros(size(phi)); i2phi = t2phi;
for p = 1:length(phi)
    phi_ = phi(p);
    t = (mod(phi_+2*pi-phi_inst,2*pi)./f_inst)/(2*pi); 

    % account for minimum delay time tmin 
    nT = (tmin-t)/T; % how many periods needed to add 
    t = t + ceil(nT)*T; 

    t2phi(p) = t;
    i2phi(p) = floor(fs*t2phi(p));
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

% details from data 
srate = handles.fSample;
nyq            = srate*0.5;  % Nyquist frequency

% filtering bound rules 
minfac         = 1;    % this many (lo)cutoff-freq cycles in filter
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
handles.PDSwin2 = ceil(.02*PDSwin); 

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

    handles.QueuedStim = timer(...
        'StartDelay', 10, ...
        'TimerFcn',   {@myPULSE, hObject}, ...
        'StopFcn',    {@finishPULSE, hObject}, ...
        'StartFcn',   {@schedulePULSE, hObject}, ...
        'UserData',   -1);

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



function txt_griddur_Callback(hObject, eventdata, handles)
% hObject    handle to txt_griddur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_griddur as text
%        str2double(get(hObject,'String')) returns contents of txt_griddur as a double


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
