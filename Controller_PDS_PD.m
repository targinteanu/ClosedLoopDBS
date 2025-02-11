function ind2Q = Controller_PDS_PD(srl, handles, dataPast)
% Based on latest serial comm, GUI-set preferences, and a chunk of recent
% data, decide which index of the forecast buffer times is the time to next
% stim, or return 0 if no stim. For now, this only supports peak (1) or
% trough (2), but in the future it should be ammended to point to any index
% in PhasesOfInterest field. 

% change this - should output index of [PeakTroughTimes] to select, or zero
% if no stim2q 

bp = norm(dataPast,2)^2/numel(dataPast); % band power surrogate

ind2Q = 0; 

if bp > 200 % min band power cutoff; orig at 1000
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
                StimMode = 'None';
            end
            if strcmpi(StimMode,'Peak')
                ind2Q = 1;
            end
            if strcmpi(StimMode,'Trough')
                ind2Q = 2;
            end
        end
    end
end

end