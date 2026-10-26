ch = floor(width(Tbl)/2);
chname = Tbl.Properties.VariableNames(ch)
x = Tbl{:,ch};

x = x-min(x); x = x/max(x); x = x*4;

t = Tbl.Time; 
t = seconds(t-t(1));

figure; plot(t, x); grid on; xlabel('s'); ylabel('V');

Fs = Tbl.Properties.SampleRate;
if isnan(Fs)
    Fs = 1/median(diff(t));
end

dq = daq('ni');
addoutput(dq, "Dev1", "ao0", "Voltage");
%%
tic
for k = 1:length(x)
    tcurr = toc; 
    icurr = round(tcurr*Fs);
    if icurr > 0 & icurr < length(x)
        x(icurr)
        write(dq, x(icurr));
    end
    pause(1/Fs);
end