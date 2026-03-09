function t = getTime_Nlx(~)

[succeeded, t] = NlxSendCommand('-GetTimestamp');

if ~succeeded
    error('Failed to get Nlx timestamp.')
end

t = str2double(t{1})/1e6; % us -> s

end