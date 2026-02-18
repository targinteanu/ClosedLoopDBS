function [HardwareFuncs, StimTriggerMode, StimLagTime] = ...
    helperGUIv1_DefHardwareFuncs(dostim)

if nargin < 1
    dostim = true;
end

RecOpts = {'Blackrock NSP'; 'Alpha Omega AlphaRS'; 'Neuro Omega'; 'Neuralynx'}; 
StimOpts = {'CereStim API'; 'CereStim Trigger / Cedrus c-pod'; 'Alpha Omega'; 'None'}; 
RecSel = listdlg("PromptString", "Recording Hardware Configuration:", ...
    "ListString",RecOpts, "SelectionMode","single");
if dostim
StimSel = listdlg("PromptString", "Stimulator Hardware Configuration:", ...
    "ListString",StimOpts, "SelectionMode","single");
end

if RecSel == 1
    % Blackrock NSP 
    HardwareFuncs = struct(...
        'SetupRecording', @connect_cbmex, ...
        'ShutdownRecording', @disconnect_cbmex, ...
        'InitRawData', @initRawData_cbmex, ...
        'GetNewRawData', @getNewRawData_cbmex, ...
        'GetTime', @getTime_cbmex); 
elseif RecSel == ((RecSel==2) || (RecSel==3))
    % AO 
    if RecSel==2
        % AlphaRS
        devicename = 'aRS';
    else
        % Neuro Omega
        devicename = 'NO';
    end
    chantypes = {'LFP', 'SPK', 'RAW', 'AI', 'SEG'};
    chantypesel = listdlg("PromptString","Channel Type:", ...
        "ListString",chantypes, "SelectionMode","single");
    chantype = chantypes{chantypesel};
    HardwareFuncs = struct(...
        'SetupRecording', @() connect_AO(devicename), ...
        'ShutdownRecording', @disconnect_AO, ...
        'InitRawData', @(a,b) initRawData_AO(a,b,chantype), ...
        'GetNewRawData', @getNewRawData_AO, ...
        'GetTime', @getTime_AO);
elseif RecSel == 4
    % NeuraLynx
else
    error('Recording option must be specified.')
end

if dostim
if StimSel == 1
    % Blackrock CereStim API
    StimTriggerMode = false;
    HardwareFuncs.SetupStimulator = @stimSetup_cerestim;
    HardwareFuncs.ShutdownStimulator = @stimShutdown_cerestim;
    HardwareFuncs.CheckConnectionStimulator = @stimCheckConnection_cerestim;
    HardwareFuncs.CalibrateStimulator = @stimCalibrate_cerestim;
    HardwareFuncs.PulseStimulator = @stimPulse_cerestim;
    HardwareFuncs.SetStimTriggerMode = @stimTriggerMode_cerestim; 
    StimLagTime = 0.03; % s 
elseif StimSel == 2
    % c-pod + CereStim in trig mode 
    StimTriggerMode = true;
    HardwareFuncs.SetupStimulator = @stimSetup_cpod;
    HardwareFuncs.ShutdownStimulator = @stimShutdown_cerestim_cpod;
    HardwareFuncs.CheckConnectionStimulator = @stimCheckConnection_cerestim;
    HardwareFuncs.CalibrateStimulator = @stimCalibrate_cerestim;
    HardwareFuncs.PulseStimulator = @stimPulse_cpod; 
    HardwareFuncs.SetStimTriggerMode = @stimTriggerMode_cerestim; 
    StimLagTime = 0.03; % s - not tested 
elseif StimSel == 3
    % AO AlphaRS or Neuro Omega
    if ~((RecSel==2) || (RecSel==3))
        error('Cannot Connect: Alpha Omega stimulator and recording device must be the same.')
    end
    StimTriggerMode = false; 
    HardwareFuncs.SetupStimulator = @stimSetup_AO;
    HardwareFuncs.ShutdownStimulator = @stimShutdown_AO;
    HardwareFuncs.CheckConnectionStimulator = @stimCheckConnection_AO;
    HardwareFuncs.CalibrateStimulator = @stimCalibrate_cerestim; % should also work for AO
    HardwareFuncs.PulseStimulator = @stimPulse_AO;
    HardwareFuncs.SetStimTriggerMode = @(~) error('Trigger mode not implemented AO.'); 
    StimLagTime = 0; % ?
else
    % dummy mode / no stimulation 
    StimTriggerMode = false;
    HardwareFuncs.SetupStimulator = @(~) 0;
    HardwareFuncs.ShutdownStimulator = @(~,~) 0;
    HardwareFuncs.CheckConnectionStimulator = @() 0;
    HardwareFuncs.CalibrateStimulator = @(~,~,~,~,~,~,~,~,~,~) 0;
    HardwareFuncs.PulseStimulator = @(~) 0;
    HardwareFuncs.SetStimTriggerMode = @(s) s;    
    StimLagTime = 0; 
end
end

end