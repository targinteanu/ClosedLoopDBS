function SerialLogOut = fixSerialLog_Memory_v0(SerialLog)

ud0 = struct('ReceivedData', '', ...
        'TrialNumber', 0, ...
        'StimOn', false, ...
        'ParadigmPhase', '', ...
        'ImageVisible', false, ...
        'TimeStamp', nan);
SerialLogOut = repmat(ud0, size(SerialLog));
ud = ud0;

for idx = 1:length(SerialLog)
    receivedData = SerialLog(idx).ReceivedData; 
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

ud.TimeStamp = SerialLog(idx).TimeStamp; 
SerialLogOut(idx) = ud;

end

end