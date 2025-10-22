Tbl = sortrows(Tbl, 'Time');
Tbl.Properties.Events = sortrows(Tbl.Properties.Events, 'Time');
t = Tbl.Time; ev = Tbl.Properties.Events;

Ts = seconds(Tbl.Properties.TimeStep);
if isnan(Ts)
    Ts = median(seconds(diff(t)));
end

d = nan(size(t));
for ei = 1:height(ev)
    [dt,di] = min(abs(seconds(t-ev.Time(ei))));
    if dt < 5*Ts
        lblvals = sscanf(ev.EventLabels(ei), '%c%c%f%c%f%c%f%s');
        d(di) = lblvals(5);
    end
end

dnan = isnan(d); dnanind = find(~dnan);
for di = 1:height(d)
    if dnan(di)
        d(di) = d(dnanind(1));
    else
        dnanind = dnanind(dnanind > di);
    end
end
figure; plot(t,d);

%% time selection
selind = (t < seconds(9250)) & (t > seconds(8000));
Tblsel = Tbl(selind,:);

%% matrix of PACs

x = Tblsel.CLFP_01; 

phif = 10:1:35; ampf = 30:5:200;
phibw = 2; ampbw = 5;
phiflt = arrayfun(@(f) buildFIRBPF(1/Ts, f-phibw, f+phibw), ...
    phif, 'UniformOutput',false);
ampflt = arrayfun(@(f) buildFIRBPF(1/Ts, f-ampbw, f+ampbw), ...
    ampf, 'UniformOutput',false);

Xa = repmat(x,1,length(ampf)); Xp = repmat(x,1,length(phif));
for a = 1:width(Xa)
    Xa(:,a) = filtfilt(ampflt{a},1,x);
end
for p = 1:width(Xp)
    Xp(:,p) = filtfilt(phiflt{p},1,x);
end

PAC = nan(length(ampf), length(phif));
for p = 1:width(PAC)
    for a = 1:height(PAC)
        PAC(a,p) = calcPAC(Xp(:,p),Xa(:,a));
    end
end

figure; imagesc(PAC); colorbar;
xticks(1:length(phif)); xticklabels(string(phif)); xlabel('Phase center freq.');
yticks(1:length(ampf)); yticklabels(string(ampf)); ylabel('Amp center freq.');
title('PAC modulation index');