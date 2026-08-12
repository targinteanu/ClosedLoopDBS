function handles = helperGUIv0_setupArtRemove(handles)

% bandaid to make compatible with v1 and v2
if isfield(handles, 'foreArgs')
    StimLagTime = handles.foreArgs.StimulatorLagTime;
else
    StimLagTime = handles.StimulatorLagTime;
end

if ~handles.MdlSetUp
    error('Model is needed for artifact removal.')
end
ARord = length(handles.Mdl)-1;

handles.ArtifactDurationSamples = ...
    ceil(handles.fSample*handles.ArtifactDuration);

handles.ArtifactStartOffsetSamples = ...
    ceil(handles.fSample*(...
        handles.ArtifactStartBefore - StimLagTime));

% init LMS wts to 0
handles.wLMS = zeros(handles.ArtifactDurationSamples, 1); 

% init Kalman P (process cov) to Q (process noise cov)
handles.kalP = zeros(ARord); handles.kalP(end,end) = handles.MdlErrVar;

end