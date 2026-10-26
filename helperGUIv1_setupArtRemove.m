function handles = helperGUIv1_setupArtRemove(handles)

handles.artRemArgs.stepsize = 0.1; % make user set instead?

if ~handles.MdlSetUp
    error('Model is needed for artifact removal.')
end
handles.artRemArgs.ARmdls = handles.foreArgs.ARmdls;
ARords = cellfun(@length, handles.foreArgs.ARmdls) - 1;

handles.artRemArgs.ArtifactDurationSamples = ...
    ceil(handles.fSample*handles.artRemArgs.StimDur);

handles.artRemArgs.ArtifactStartOffsetSamples = ...
    ceil(handles.fSample*(...
        handles.artRemArgs.ArtifactStartBefore - handles.foreArgs.StimulatorLagTime));

handles.artRemArgs.kalP = cell(size(ARords)); 
handles.artRemArgs.wLMS = cell(size(ARords));
for ch = 1:length(ARords)

    % init LMS wts to 0
    handles.artRemArgs.wLMS{ch} = zeros(handles.artRemArgs.ArtifactDurationSamples, 1); 

    % init Kalman P (process cov) to Q (process noise cov)
    handles.artRemArgs.kalP{ch} = zeros(ARords(ch)); 
    handles.artRemArgs.kalP{ch}(end,end) = handles.artRemArgs.MdlErrVar{ch};
end

end