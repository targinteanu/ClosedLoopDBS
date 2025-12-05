function StimArgs = stimSetup_AO(StimSetupArgs)

% StimSetupArgs comes from a selection of the user's args from the GUI 
% using getStimSetupArgs. It should contain fields channel1 & channel2
% (annode & cathode; mat be lists), amp1 & amp2 (first & second phase ampl 
% [uA]), width1 & width2 (first & second phase width [us]), interphase
% (second phase delay + first phase width [us]), frequency [hz], pulses
% [#]. 

% adjust amplitude based on # of channels 
N1 = length(StimSetupArgs.channel1); 
N2 = length(StimSetupArgs.channel2);
if (N1 > 1) || (N2 > 1)
    error('Multi-channel stimulation on AO is not yet supported.')
end

% StimSetupArgs was made for CereStim; convert these values to be more
% useful for AO
dur = StimSetupArgs.pulses/StimSetupArgs.frequency; % s
del = StimSetupArgs.interphase - StimSetupArgs.width1; % us
StimArgs = struct(...
    'StimChannel', StimSetupArgs.channel2, ...
    'ReturnChannel', StimSetupArgs.channel1, ...
    'FirstPhaseAmpl_mA',  StimSetupArgs.amp1/1000, ...
    'SecondPhaseAmpl_mA', -StimSetupArgs.amp2/1000, ...
    'FirstPhaseWidth_ms',  StimSetupArgs.width1/1000, ...
    'SecondPhaseWidth_ms', StimSetupArgs.width2/1000, ...
    'SecondPhaseDelay_ms', del/1000, ...
    'Freq_hz', StimSetupArgs.frequency, ...
    'Duration_sec', dur);

end
