function CharSerialCallbackReceiver_PD_v0(hSrlPrt, evt, hTxtBox, hStatus)
% by Toren Arginteanu, 7/5/2024
% Callback function to be assigned to a Serial Port object as follows: 
% >> hSrlPrt = serialPort( "COM??" , 9600 ); 
% >> configureCallback( hSrlPrt , "terminator" , @CharSerialCallback );
%   -- OR --
% >> configureCallback( hSrlPrt , "terminator" , ... 
%             @(hSrlPrt,evt) CharSerialCallback(hSrl,evt,hTxtBox) );
% hTxtBox: handle to text box that has field 'String'
% hTbl: handle to display table that has field 'Data'
% Received lines of data will be stored in hSrlPrt object's UserData field
% as char arrays. The following keywords will also trigger actions: 
%   'TEST' -> respond with 'COPY'
%   'COPY' -> prints elapsed time from stopwatch [i.e., toc]
%   'FLUSH' -> flushes hSrlPrt buffer
%   'DELETE' -> deletes and clears hSrlPrt
%   'EVT' -> displays info about event evt

%% handle inputs 
output_TxtBox = ~isempty(hTxtBox);
output_Status = ~isempty(hStatus);

%% read serial 
receivedData = readline(hSrlPrt); 
receivedData = char(receivedData); 
dspMsg = ['Serial Port ',char(hSrlPrt.Port),' Received Message: ',receivedData];
disp(dspMsg);
if output_TxtBox
    hTxtBox.String = dspMsg;
end

%% set user data 
ud = hSrlPrt.UserData;
flds = {'ReceivedData', ...
        'ParadigmPhase'};
ud.ReceivedData = receivedData;

if contains(receivedData, 'STATUS')
    % phase of experiment 
    ph = sscanf(receivedData, 'STATUS %s');
    ud.ParadigmPhase = ph;
end

%{
% check fields; handle empty 
for fld = flds
    fld_ = fld{1};
    if ~isfield(ud,fld_)
        ud.(fld_) = [];
    end
end
%}

hSrlPrt.UserData = ud;
if output_Status
    hStatus.String = ud.ParadigmPhase;
end

%% special behaviors 
if strcmpi(receivedData, 'TEST')
    writeline(hSrlPrt, 'COPY');
end
if strcmpi(receivedData, 'COPY')
    toc
end
if strcmpi(receivedData, 'FLUSH')
    flush(hSrlPrt);
end
%{
if strcmpi(receivedData, 'DELETE')
    delete(hSrlPrt);
end
%}
if strcmpi(receivedData, 'EVT')
    evt
end

end