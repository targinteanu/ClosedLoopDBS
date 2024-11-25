function [t2Q, Stim2Q] = Controller_PDS_PD(handles, dataPast, PeakTroughTimes)

bp = norm(dataPast,2)^2/numel(dataPast); % band power surrogate
t2peak = PeakTroughTimes(:,1); 
t2trou = PeakTroughTimes(:,2);

Stim2Q = false;
t2Q = inf; 

if bp > 10 % min band power cutoff; orig at 1000
    if handles.StimActive
        ParadigmPhase = handles.srl.UserData.ParadigmPhase;
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
                t2Q = t2peak;
                Stim2Q = true;
            end
            if strcmpi(StimMode,'Trough')
                t2Q = t2trou;
                Stim2Q = true;
            end
        end
    end
end

end