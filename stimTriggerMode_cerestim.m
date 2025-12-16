function stimulator = stimTriggerMode_cerestim(stimulator, setupargs)

% setup the CereStim with desired stim params 
if isempty(stimulator)
    stimulator = stimSetup_cerestim(setupargs);
end

% set the CereStim stimulator to trigger mode
stimulator.trigger(1); % rising edge 

end