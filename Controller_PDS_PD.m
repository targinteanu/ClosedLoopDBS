function ind2Q = Controller_PDS_PD(srl, handles, bpData, bpThresh)
% Based on latest serial comm, GUI-set preferences, and a chunk of recent
% data, decide which index of the forecast buffer times is the time to next
% stim, or return 0 if no stim. For now, this only supports peak (1) or
% trough (2), but in the future it should be ammended to point to any index
% in PhasesOfInterest field. 

if nargin < 4
    bpThresh = 1000; % min band power cutoff; orig at 1000
end

if numel(bpData) == 1
    bp = bpData;
else
    bp = norm(bpData,2)^2/numel(bpData); % band power surrogate
end

ind2Q = 0; 

if bp > bpThresh
    if handles.StimActive

        ParadigmPhase = srl.UserData.ParadigmPhase;
        if ~strcmpi(ParadigmPhase,'WAIT')
            if strcmpi(ParadigmPhase, 'Started') || strcmp(ParadigmPhase, 'gray')
                % Started, gray, and red should all be the same.
                ParadigmPhase = 'red';
            end
            try
                StimMode = getfield(handles.StimMode, ParadigmPhase);
            catch
                warning(['ParadigmPhase ',ParadigmPhase,' unrecognized.'])
                StimMode = 0;
            end
            ind2Q = StimMode;
        end
    end
end

end