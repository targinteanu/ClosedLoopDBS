function stimulator = stimPulse_cerestim(stimulator, ~)
% Send a single stimulation pulse to the cerestim that was predefined in
% the stimSetup function. Adapted from myPULSE in verson 0

%% warnings 
% Is this necessary? Will removing improve speed? 
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

%% Stimulation Pulse 
stimulator.play(1); 
% consider changing to groupStimulus or manualStim to save time

end