function connect_AO(device)

if nargin < 1
    error('Device name must be specified for AO.')
end

if strcmp(device, 'aRS')
    MAC_addr = 'C8:DF:84:F8:9F:E2'; % MAC address of Functional Neuro Lab's AlphaRS 
elseif strcmp(device, 'NO')
    MAC_addr = 'bc:6a:29:d3:4a:63'; % MAC address of Zayed 3 Neuro Omega 
else
    disp(device)
    error('Device name above not recognized.')
end

Result = AO_DefaultStartConnection(MAC_addr);

if Result == 10
    % Dont need to do anything because connection is already open and it
    % already sends a message stating that
elseif Result
    error(['Connection to AO failed with code ',num2str(Result)])
end

if ~AO_IsConnected()
    error('Did not connect.')
end

end