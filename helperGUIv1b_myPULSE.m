function handles = helperGUIv1b_myPULSE(hTimer, eventdata, handles)

% perform the pulse and measure the time 
initTic = handles.initTic;
stimtime1 = handles.HardwareFuncs.GetTime(initTic);
handles.stimulator = handles.HardwareFuncs.PulseStimulator(handles.stimulator);
stimtime2 = handles.HardwareFuncs.GetTime(initTic);

% disp time of pulse using eventdata
eventTime = datestr(eventdata.Data.time);
stimtime = .5*(stimtime1 + stimtime2);
dstimtime = stimtime2 - stimtime1; 
stimschedtime = hTimer.UserData; 
disp(['Stimulus pulsed at ',eventTime,' within ',num2str(dstimtime),'s, ',...
      num2str(stimtime - stimschedtime),' s late'])

% record stim timing in the GUI 
handles.stimLastTime = stimtime; handles.stimNewTime = stimtime; 
if handles.stP1 <= length(handles.stStorage1)
    handles.stStorage1(handles.stP1) = stimtime; 
    handles.stP1 = handles.stP1 + 1;
else
    % storage full; save
    StimTime = handles.stStorage1;
    svfn = [handles.SaveFileLoc,filesep,'SaveFile',num2str(handles.SaveFileN),'.mat'];
    disp(['Saving Stimulus to ',svfn])
    save(svfn,'StimTime');
    handles.SaveFileN = handles.SaveFileN + 1;
    handles.stP1 = 1;
    handles.stStorage1 = nan(size(handles.stStorage1));
end

% note stim timing for artifact removal
handles.artRemArgs.StimTimes = {stimtime + handles.StimulatorLagTime};

end