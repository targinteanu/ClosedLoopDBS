function connect_Nlx()

%Load NetCom into MATLAB, and connect to the NetCom server
%If you are running MATLAB on a PC other than the one with the NetCom
%server, you will need to change the server name to the name of the server
%PC.
serverName = 'localhost';
fprintf('Connecting to %s...', serverName);
succeeded = NlxConnectToServer(serverName);
if succeeded ~= 1
    error('FAILED to connect. Exiting script.');
else
    fprintf('Connect successful.\n');
end

serverIP = NlxGetServerIPAddress();
fprintf('Connected to IP address: %s\n', serverIP);

serverPCName = NlxGetServerPCName();
fprintf('Connected to PC named: %s\n', serverPCName);

serverApplicationName = NlxGetServerApplicationName();
fprintf('Connected to the NetCom server application: %s\n', serverApplicationName);

%Identify this program to the server we're connected to.
succeeded = NlxSetApplicationName('ClosedLoopDBS Matlab Script');
if succeeded ~= 1
    error('FAILED to set the application name');
end

end