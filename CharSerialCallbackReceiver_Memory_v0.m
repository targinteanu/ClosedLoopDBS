function CharSerialCallbackReceiver_Memory_v0(hSrlPrt, evt, hTxtBox, hTbl)
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
output_TxtBox = (nargin >= 3) && ~isempty(hTxtBox);
output_Tbl    = (nargin >= 4) && ~isempty(hTbl);

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
        'TrialNumber', ...
        'StimOn', ...
        'ParadigmPhase', ...
        'ImageVisible'};
ud.ReceivedData = receivedData;
 
if contains(receivedData, 'Trial') && contains(receivedData, 'stimon')
    % trial and stimon
    trl_stim = sscanf(receivedData, 'Trial %f of %f; stimon %f');
    ud.TrialNumber = trl_stim(1); 
    ud.StimOn = logical(trl_stim(3));

elseif contains(receivedData, 'Photo')
    % photo visible? 
    if contains(receivedData, 'Shown')
        ud.ImageVisible = true;
    elseif contains(receivedData, 'Gone')
        ud.ImageVisible = false;
    end
    
else
    % phase of experiment 
    expPhase = {'ENCODE', 'DECODE', 'HOLD'};
    noPhase = 'WAIT';
    for ep = expPhase
        ep_ = ep{1};
        if contains(receivedData, ep_)
            if contains(receivedData, 'Start')
                ud.ParadigmPhase = ep_;
            elseif contains(receivedData, 'End')
                ud.ParadigmPhase = noPhase;
            end
            break
        end
    end
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
if output_Tbl
    hTbl.Data = {ud.TrialNumber, ud.StimOn, ud.ParadigmPhase, ud.ImageVisible};
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