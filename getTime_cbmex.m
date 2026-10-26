function t = getTime_cbmex(~)

% TO DO: the time units may be different (microseconds?) on the gemini!

t = cbmex('time');

end