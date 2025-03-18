function t = getTime_AO(~)

[Result, t] = AO_GetLatestTimeStamp(); 
t = double(t)/(3.02e4); % TO DO: make sure this is correct time in seconds!

if Result
    error(['Failed to get latest time stamp: code ',num2str(Result)])
end

end