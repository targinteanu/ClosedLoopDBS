function handles = helperGUIv1b_MainStim(handles)

Stim2Q = false;
forBuff = handles.recDataStructs.forBuffs{1}; 
forBuffNew = [max(forBuff(:,1)), max(forBuff(:,2))]; 
timeBuffs = handles.recDataStructs.timeBuffs;
forBuffNew = forBuffNew - timeBuffs{handles.channelIndex}(end,:); % [t2p, t2t]
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

end