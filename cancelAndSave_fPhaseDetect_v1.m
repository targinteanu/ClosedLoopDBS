function [svN, fut] = cancelAndSave_fPhaseDetect_v1(fut, svname, svN)
% Cancel a future fut running bg_PhaseDetect, but save the buffer data. 

try
    bgArgOut = fetchOutputs(fut);
catch ME
    if strcmpi(ME.identifier, 'MATLAB:parallel:future:ExecutionErrors') && ...
            strcmpi(ME.message, 'One or more futures resulted in an error.')
        bgArgOut = [];
    else
        % are there other possibilities? 
        keyboard
    end
end
cancel(fut);

if ~isempty(bgArgOut)
    ForStimSv = bgArgOut; 
    svN = savebg_PhaseDetect_v1(ForStimSv, [], [], svname, svN);
end

end