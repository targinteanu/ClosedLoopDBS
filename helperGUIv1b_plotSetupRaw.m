function [handles, fltPlt, forPlt, forBuff, tSt, common_xlim, unitname] = ...
    helperGUIv1b_plotSetupRaw(handles, tNow)

    % Check which channel is selected and get some data to plot
    handles.channelIndex = get(handles.pop_channels,'Value'); 
    % Now we know the sampling rate of the selected channel
    handles.fSample = handles.fSamples(handles.channelIndex);
    handles.bufferSize = ceil(str2double(get(handles.txt_display,'String')) * handles.fSample);
    handles.bufferSizeGrid = ceil(str2double(get(handles.txt_griddur,'String')) * handles.fSamples);

    % assign channel indexes (NOT IDs!) to use for filtering and forecasting
    chInd = handles.channelIndex;
    selRaw2Flt = []; selFlt2For = []; selFor2Art = [];
    selRaw2Art = chInd;
    if handles.FilterSetUp
        selRaw2Flt = chInd;
        if handles.MdlSetUp
            selFlt2For = 1;
            selFor2Art = 1;
        end
    end
    selRaw2For = [];
    handles.selInds = struct(...
        'selRaw2Flt', selRaw2Flt, ...
        'selFlt2For', selFlt2For, ...
        'selFor2Art', selFor2Art, ...
        'selRaw2Art', selRaw2Art, ...
        'selRaw2For', selRaw2For);

    % (re-)init data structs
    buffSize = handles.bufferSizeGrid .* ones(size(handles.allChannelInfo)); % samples
    buffSize(chInd) = handles.bufferSize;
    if handles.FilterSetUp
        filtOrds = [handles.FilterOrder]; % array with chans as cols
    else
        filtOrds = [];
    end
    [rawD, artD, fltD, forD, timeBuffs, initTic] = ...
        InitializeRecording(handles.HardwareFuncs.InitRawData, ...
            buffSize, filtOrds, handles.PDSwin1, ...
            handles.allChannelIDs, selRaw2Art, selRaw2Flt, selRaw2For, selFlt2For);
    rawD1 = rawD(1,:); rawD4 = rawD(4,:);
    artD1 = artD(1,:); artD4 = artD(4,:);
    fltD1 = fltD(1,:); fltD4 = fltD(4,:);
    forD1 = forD(1,:); forD4 = forD(4,:);
    timeBuff = timeBuffs{chInd};
    buffSize2 = (handles.bufferSize / rawD1{chInd}.SampleRate) * .5 * handles.stimMaxFreq;
    buffSize2 = ceil(buffSize2);
    forBuff = nan(buffSize2, length(handles.PhaseOfInterest));
    tSt = nan(buffSize2,1);
    handles.allChannelInfo = rawD1;
    handles.channelIDlist = cellfun(@(ch) ch.IDnumber, handles.allChannelInfo);
    rawPlt = data2timetable(rawD4(chInd),rawD1(chInd),handles.time0); rawPlt = rawPlt{1};
    fltPlt = data2timetable(fltD4,fltD1,handles.time0); fltPlt = fltPlt{1};
    forPlt = data2timetable(forD4,forD1,handles.time0); forPlt = forPlt{1};
    artPlt = data2timetable(artD4,artD1,handles.time0); artPlt = artPlt{1};
    recDataStructs.forBuffs = {forBuff}; recDataStructs.stimBuff = tSt;
    for v = ["rawD", "artD", "fltD", "forD", "timeBuffs", "initTic"]
        eval("recDataStructs."+v+" = "+v+";");
    end
    handles.recDataStructs = recDataStructs;

    unitname = handles.allChannelInfo{handles.channelIndex}.Unit;

    % keep track of the display time 
    handles.timeDisp1 = tic; handles.timeDisp0 = timeBuff(end, :);
    handles.timeDispBuff = nan(size(timeBuff));

    % time axis 
    if isempty(rawPlt)
        error('Raw data was not initialized properly.')
    end
    tRaw = rawPlt.Time - tNow; 
    if ~handles.check_polar.Value
        % do not track time exactly 
        tRaw = seconds(((-length(tRaw)+1):0)/handles.fSample);
    end
    tPltRng = tRaw; 
    tPltRng = [min(tPltRng), max(tPltRng)];
    tPltRng = tPltRng + [-1,1]*.1*diff(tPltRng); % do we actually want this extra space? 

    % initiate raw data plot
    axes(handles.ax_raw);
    hold off; 
    handles.h_rawDataTrace = plot(tRaw, rawPlt.Variables);
    grid on; title('Raw Channel Data'); 
    xlabel('time'); ylabel(rawPlt.Properties.VariableNames{1});
    if sum(~isnan(tPltRng))
        common_xlim = tPltRng; 
        xlim(common_xlim);
    else
        common_xlim = xlim();
    end
    % add artifact removal if applicable 
    if handles.FilterSetUp && handles.MdlSetUp
        if handles.check_artifact.Value
            if isempty(artPlt)
                set(handles.check_artifact,'Value',false);
                handles.check_artifact_Value = false;
                error('Artifact removal was not actually set up. Something is wrong in the code.')
            end
            hold on;
            handles.h_artDataTrace = plot(artPlt.Time - tNow, artPlt.Variables, ':');
        end
    end

    % initiate timing stem plot
    %if handles.check_polar.Value
        tStem = handles.time0 + seconds(timeBuff) - tNow; % ?
    %else
    %    tStem = [seconds(nan), tRaw]; % TO DO: fix data2timetable eating one sample, then get rid of the nan
    %end
    axes(handles.ax_timing); hold off;
    handles.h_timingTrace = ...
        stem(tStem, [nan; diff(timeBuff)], 's');
    hold on;
    %{
    handles.h_timeDispTrace = ...
        stem(handles.time0 + seconds(handles.timeDispBuff) - tNow, ...
        [nan; diff(handles.timeDispBuff)], 'o');
    %}
    grid on; title('Update Duration'); xlabel('time (s)'); ylabel('s');
    xlim(common_xlim);

end