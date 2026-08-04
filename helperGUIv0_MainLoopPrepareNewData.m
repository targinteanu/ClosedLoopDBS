function [handles, newContinuousData] = helperGUIv0_MainLoopPrepareNewData(handles, newContinuousData, time)

    N = length(newContinuousData);

    %{
    % artifact removal (1)
    if handles.FilterSetUp
    if handles.MdlSetUp
    if handles.check_artifact.Value
        if numel(handles.artReplaceRemaining)
            try
            % continue replacing from last loop iter
            artLen = min(length(handles.artReplaceRemaining), N);
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
    %}

    rawDataBuffer = bufferData(handles.rawDataBuffer, newContinuousData);
    % new data starts at rawDataBuffer(end-N+1)
    t0 = handles.lastSampleProcTime; 
    handles.lastSampleProcTime = ...
        time + (N-1)/handles.fSample;
    T = handles.lastSampleProcTime - t0;
    if T < 0
        %warning('Reported sample time is negative.')
    end
    diffSampleProcTime = nan(size(newContinuousData)); diffSampleProcTime(end) = T;
    handles.diffSampleProcTime = bufferData(handles.diffSampleProcTime, diffSampleProcTime);

    % artifact removal (2) 
    if handles.FilterSetUp
    if handles.MdlSetUp
    if handles.check_artifact.Value
        stimind = handles.stimind - N; % N samples have passed
        if stimind > 0 % rel to start of buffer
            try
            rawOffset = mean(rawDataBuffer);
            rawDataBuffer = rawDataBuffer - rawOffset;
            artStart = stimind - handles.ArtifactStartOffsetSamples;
                % rel to start of buffer
            artDur = handles.ArtifactDurationSamples;
            Mdl = handles.Mdl; 
            wAR = -Mdl(2:end)/Mdl(1); wAR = fliplr(wAR); 
            preL = max(artDur, length(wAR)); % extra length needed for art rem
            dataStartIdx = min(handles.bufferSize-N+1, artStart) -preL;
            dataPadded = dataStartIdx < 1;
            if dataPadded
                padL = 1-dataStartIdx;
                rawDataBuffer = [rawDataBuffer(1,:)*ones(padL,1); rawDataBuffer];
                dataStartIdx = 1; artStart = artStart+padL;
            end
            artIdx = dataStartIdx:min((artStart+artDur), handles.bufferSize); % to replace
                % change endpoint to buffer end to make smoother?
            artStart = artStart - dataStartIdx +1; % rel to replace start
            stepsize = 0.5; % make user set instead? 
            if numel(artIdx)
            [rawDataBuffer(artIdx), handles.kalP, handles.wLMS] = iterKalmanLMS(...
                rawDataBuffer(artIdx), artStart, wAR, ...
                handles.kalP, handles.MdlErrVar, handles.wLMS, stepsize, true);
            end
            if dataPadded
                rawDataBuffer = rawDataBuffer((padL+1):end,:);
            end
            rawDataBuffer = rawDataBuffer + rawOffset;
            %{
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
            %}
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

    handles.rawDataBuffer = rawDataBuffer;
end