function t = getTime_Nlx(~)

[succeeded, t] = NlxSendCommand('-GetTimestamp');

if ~succeeded
    error('Failed to get Nlx timestamp.')
end

t = t/1e6; % uS -> s

end