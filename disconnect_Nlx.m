function disconnect_Nlx()

[succeeded, ~] = NlxSendCommand('-StopAcquisition');
if succeeded == 0
    warning('Failed to stop acquisition');
end

%get a list of all objects in the DAS, along with their types.
[succeeded, dasObjects, dasTypes] = NlxGetDASObjectsAndTypes;
if succeeded == 0
    warning('FAILED get DAS objects and types; cannot close streams');
else
  
%close all open streams before disconnecting
for index = 1:length(dasObjects)
     if strcmp(char(dasTypes(index)), 'AcqSource') ~= 1 
        succeeded = NlxCloseStream(dasObjects(index));
        if succeeded == 0
            warning(['FAILED to close stream for ', char(dasObjects(index))]);
        end
     end
end;

%Disconnects from the server and shuts down NetCom
succeeded = NlxDisconnectFromServer();
if succeeded ~= 1
    error('FAILED disconnect from server');
else
    disp('Disconnected Nlx from server.');
end

end