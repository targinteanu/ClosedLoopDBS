function connect_cbmex()

try
    cbmex('open');
catch ME
    if strcmp(ME.identifier,'MATLAB:unassignedOutputs')
        % Dont need to do anything because cbmex is already open and it
        % already sends a message stating that
    elseif strcmp(ME.message,'Unknown error')
        error(['Unknown cbmex error; ' ...
            'ensure the correct cerebus SDK is in MATLAB path.'])
    else
        rethrow(ME)
    end
end

cbmex('trialconfig',1,'absolute','double')

end