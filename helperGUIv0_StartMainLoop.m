function handles = helperGUIv0_StartMainLoop(handles)

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

    handles.rawDataBuffer = bufferData(handles.rawDataBuffer, newContinuousData);
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
        [filtData, handles.filtCond] = filter(...
            handles.BPF, 1, newContinuousData, filtInitCond);
        handles.filtDataBuffer = bufferData(handles.filtDataBuffer, filtData);
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

function plotData = plotLogical(logData)
% take in a logical array and output an array that will plot true as 1 
% and will not plot false
plotData = double(logData); 
plotData(~logData) = nan; 
end

end