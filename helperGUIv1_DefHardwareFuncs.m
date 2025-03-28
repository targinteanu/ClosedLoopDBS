function [HardwareFuncs, StimTriggerMode] = helperGUIv1_DefHardwareFuncs()

RecOpts = {'Blackrock NSP'; 'Alpha Omega AlphaRS'}; % TO DO: add Neuro Omega ?
StimOpts = {'CereStim API'; 'CereStim Trigger / Cedrus c-pod'; 'AlphaRS'; 'None'}; 
RecSel = listdlg("PromptString", "Recording Hardware Configuration:", ...
    "ListString",RecOpts, "SelectionMode","single");
StimSel = listdlg("PromptString", "Stimulator Hardware Configuration:", ...
    "ListString",StimOpts, "SelectionMode","single");
if RecSel == 1
    % Blackrock NSP 
    HardwareFuncs = struct(...
        'SetupRecording', @connect_cbmex, ...
        'ShutdownRecording', @disconnect_cbmex, ...
        'InitRawData', @initRawData_cbmex, ...
        'InitRecording', @InitializeRecording_cbmex, ...
        'GetNewRawData', @getNewRawData_cbmex, ...
        'GetTime', @getTime_cbmex); 
elseif RecSel == 2
    % AO AlphaRS; might also work with Neuro Omega (untested) 
    HardwareFuncs = struct(...
        'SetupRecording', @connect_AO, ...
        'ShutdownRecording', @disconnect_AO, ...
        'InitRawData', @initRawData_AO, ...
        'InitRecording', @InitializeRecording_AO, ...
        'GetNewRawData', @getNewRawData_AO, ...
        'GetTime', @getTime_AO);
else
    error('Recording option must be specified.')
end
if StimSel == 1
    % Blackrock CereStim API
    StimTriggerMode = false;
    HardwareFuncs.SetupStimulator = @stimSetup_cerestim;
    HardwareFuncs.ShutdownStimulator = @stimShutdown_cerestim;
    HardwareFuncs.CheckConnectionStimulator = @stimCheckConnection_cerestim;
    HardwareFuncs.CalibrateStimulator = @stimCalibrate_cerestim;
    HardwareFuncs.SetupStimTTL = @(~) []; % no TTL enabled
    HardwareFuncs.PulseStimulator = @stimPulse_cerestim;
    HardwareFuncs.SetStimTriggerMode = @stimTriggerMode_cerestim;  
elseif StimSel == 2
    % c-pod + CereStim in trig mode 
    StimTriggerMode = true;
    HardwareFuncs.SetupStimulator = @stimSetup_cerestim;
    HardwareFuncs.ShutdownStimulator = @stimShutdown_cerestim;
    HardwareFuncs.CheckConnectionStimulator = @stimCheckConnection_cerestim;
    HardwareFuncs.CalibrateStimulator = @stimCalibrate_cerestim;
    HardwareFuncs.SetupStimTTL = @srlSetup_cpod;
    HardwareFuncs.PulseStimulator = @stimPulse_cpod; 
    HardwareFuncs.SetStimTriggerMode = @stimTriggerMode_cerestim;       
elseif StimSel == 3
    % AO AlphaRS; might also work with Neuro Omega (untested) 
    StimTriggerMode = false; 
    HardwareFuncs.SetupStimTTL = @(~) []; % no TTL enabled
    HardwareFuncs.SetStimTriggerMode = @(s) s; 
else
    % dummy mode / no stimulation 
    StimTriggerMode = false;
    HardwareFuncs.SetupStimulator = @(~) 0;
    HardwareFuncs.ShutdownStimulator = @(~,~) 0;
    HardwareFuncs.CheckConnectionStimulator = @() 0;
    HardwareFuncs.CalibrateStimulator = @(~,~,~,~,~,~,~,~,~,~) 0;
    HardwareFuncs.SetupStimTTL = @(~) []; % no TTL enabled
    HardwareFuncs.PulseStimulator = @(~,~) 0;
    HardwareFuncs.SetStimTriggerMode = @(s) s;                          
end

end