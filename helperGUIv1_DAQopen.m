function handles = helperGUIv1_DAQopen(handles)
% TO DO: add pop_channel3-5 to memory GUI

if handles.DAQstatus
    error('DSP already connected. Please disconnect first.')
end

handles.HardwareFuncs.SetupRecording();
handles.time0 = datetime - seconds(handles.HardwareFuncs.GetTime(handles.initTic));

handles.DAQstatus = true;

% Acquire some data to get channel information. Determine which channels
% are enabled
[~,~,~,allChannelInfo] = handles.HardwareFuncs.InitRawData(handles.allChannelIDs, handles.bufferSizeGrid);
handles.allChannelInfo = allChannelInfo;

handles.HardwareFuncs.ShutdownRecording();

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
ncol = 3; % # of columns
nrow = 21;
img = nan(nrow, ncol);
[X,Y] = meshgrid(1:ncol, 1:nrow);
handles.elecGridImg = imagesc(img, [gridminval, gridmaxval]); 
xticks([]); yticks([]);
colormap('parula'); colorbar; hold on;
chL_ = handles.channelList(1:min((ncol*nrow), length(handles.channelList)));
X = X(:); X = X(1:length(chL_));
Y = Y(:); Y = Y(1:length(chL_));
text(X,Y, chL_, ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', ...
    'FontWeight', 'bold', ...
    'Color',[.8 0 0]);

% Set the Start/Stop toggle button to stopped state (String is 'Start' and
% Value is 1)
set(handles.tgl_StartStop,'String','Start', 'Value',0)

end