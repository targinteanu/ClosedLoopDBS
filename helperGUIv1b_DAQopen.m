function handles = helperGUIv1b_DAQopen(handles, nrow, ncol)
% TO DO: add pop_channel3-5 to memory GUI

if nargin < 2
    nrow = []; ncol = [];
end
if isempty(nrow)
    % default to 21-by-3 grid 
    nrow = 21;
end
if isempty(ncol)
    % default to 21-by-3 grid
    ncol = 3;
end

if handles.DAQstatus
    error('DSP already connected. Please disconnect first.')
end

handles.HardwareFuncs.SetupRecording();
pause(.1);
try
handles.time0 = datetime - seconds(handles.HardwareFuncs.GetTime(handles.initTic));
catch ME
    % restart
    pause(.1);
    handles.HardwareFuncs.ShutdownRecording();
    pause(.1);
    handles.HardwareFuncs.SetupRecording();
    pause(.1);
    handles.time0 = datetime - seconds(handles.HardwareFuncs.GetTime(handles.initTic));
end

handles.DAQstatus = true;

% Acquire some data to get channel information. Determine which channels
% are enabled
[~,~,~,allChannelInfo] = handles.HardwareFuncs.InitRawData(handles.allChannelIDs, handles.bufferSizeGrid);
handles.allChannelInfo = allChannelInfo;

%handles.HardwareFuncs.ShutdownRecording();

% set channel popup menu to hold channels
handles.fSamples = cellfun(@(ch) ch.SampleRate, allChannelInfo);
handles.channelIDlist = cellfun(@(ch) ch.IDnumber, allChannelInfo);
handles.allChannelIDs = handles.channelIDlist; % resets channel selection 
handles.channelList = cellfun(@(ch) [num2str(ch.IDnumber),': ',ch.Name], ...
    allChannelInfo, 'UniformOutput',false);
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
handles.elecGridImg = helperGUIv1_ElectrodeGridInit(handles.ax_elecgrid, ...
    handles.channelList, nrow, ncol, gridminval, gridmaxval);

% Set the Start/Stop toggle button to stopped state (String is 'Start' and
% Value is 1)
set(handles.tgl_StartStop,'String','Start', 'Value',0)

end