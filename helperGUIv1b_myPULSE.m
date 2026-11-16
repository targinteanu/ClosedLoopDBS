function handles = helperGUIv1b_myPULSE(hTimer, eventdata, handles, hFigure)

% perform the pulse and measure the time 
initTic = handles.initTic;
stimtime1 = handles.HardwareFuncs.GetTime(initTic);
handles.stimulator = handles.HardwareFuncs.PulseStimulator(handles.stimulator);
stimtime2 = handles.HardwareFuncs.GetTime(initTic);

% disp time of pulse using eventdata
stimtime = .5*(stimtime1 + stimtime2);
%{
eventTime = datestr(eventdata.Data.time);
dstimtime = stimtime2 - stimtime1; 
stimschedtime = hTimer.UserData; 
disp(['Stimulus pulsed at ',eventTime,' within ',num2str(dstimtime),'s, ',...
      num2str(stimtime - stimschedtime),' s late'])
%}

% record stim timing in the GUI 
hFigure.UserData.stimLastTime = stimtime; hFigure.UserData.stimNewTime = stimtime; 
if hFigure.UserData.stP1 <= length(hFigure.UserData.stStorage1)
    hFigure.UserData.stStorage1(hFigure.UserData.stP1) = stimtime; 
    hFigure.UserData.stP1 = hFigure.UserData.stP1 + 1;
else
    % storage full; save
    StimTime = hFigure.UserData.stStorage1;
    svfn = [handles.SaveFileLoc,filesep,'SaveFile',num2str(handles.SaveFileN),'.mat'];
    disp(['Saving Stimulus to ',svfn])
    save(svfn,'StimTime');
    handles.SaveFileN = handles.SaveFileN + 1;
    hFigure.UserData.stP1 = 1;
    hFigure.UserData.stStorage1 = nan(size(hFigure.UserData.stStorage1));
end

% note stim timing for artifact removal
handles.artRemArgs.StimTimes = {stimtime + handles.foreArgs.StimulatorLagTime};

end