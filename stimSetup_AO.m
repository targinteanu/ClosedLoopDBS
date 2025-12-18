function StimArgs = stimSetup_AO(StimSetupArgs)

% StimSetupArgs comes from a selection of the user's args from the GUI 
% using getStimSetupArgs. It should contain fields channel1 & channel2
% (annode & cathode; mat be lists), amp1 & amp2 (first & second phase ampl 
% [uA]), width1 & width2 (first & second phase width [us]), interphase
% (second phase delay + first phase width [us]), frequency [hz], pulses
% [#]. 

% adjust amplitude based on # of channels 
ch1 = StimSetupArgs.channel1; ch2 = StimSetupArgs.channel2;
ch1 = ch1(~isnan(ch1)); ch2 = ch2(~isnan(ch2));
N1 = length(ch1); 
N2 = length(ch2);
if (N1 > 1) || (N2 > 1)
    error('Multi-channel stimulation on AO is not yet supported.')
end

% StimSetupArgs was made for CereStim; convert these values to be more
% useful for AO
dur = StimSetupArgs.pulses/StimSetupArgs.frequency; % s
dur = .9*dur;
del = StimSetupArgs.interphase - StimSetupArgs.width1; % us
StimArgs = struct(...
    'StimChannel', ch2, ...
    'ReturnChannel', ch1, ...
    'FirstPhaseAmpl_mA',  StimSetupArgs.amp1/1000, ...
    'SecondPhaseAmpl_mA', -StimSetupArgs.amp2/1000, ...
    'FirstPhaseWidth_ms',  StimSetupArgs.width1/1000, ...
    'SecondPhaseWidth_ms', StimSetupArgs.width2/1000, ...
    'SecondPhaseDelay_ms', del/1000, ...
    'Freq_hz', StimSetupArgs.frequency, ...
    'Duration_sec', dur);

end
