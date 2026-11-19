function t = getTime_AO(~)

[Result, t] = AO_GetLatestTimeStamp(); 
t = double(t);
t = t/(3e4); % is this correct time in seconds?
% should be in 30kHz system clocks 

if Result
    error(['Failed to get latest time stamp: code ',num2str(Result)])
end

end