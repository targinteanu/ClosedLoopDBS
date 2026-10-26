function [phiTbl, fTbl, ATbl] = instPhaseFreqTbl(tbl)
% 
% get instantaneous phase, frequency, & amp at each timepoint in a timetable 
% 

fs = tbl.Properties.SampleRate; 
if isnan(fs)
    fs = 1/mean(seconds(diff(tbl.Time)));
end
[phi_block, f_block, A_block] = instPhaseFreq(tbl.Variables, fs);

phiTbl = tbl; phiTbl.Variables = phi_block; 
for c = 1:width(phiTbl)
    phiTbl.Properties.VariableUnits{c} = 'Radians';
    phiTbl.Properties.VariableNames{c} = ...
        [tbl.Properties.VariableNames{c},' phase'];
end
phiTbl.Properties.Description = [tbl.Properties.Description,' phase'];

fTbl = tbl; fTbl.Variables = f_block; 
for c = 1:width(fTbl)
    fTbl.Properties.VariableUnits{c} = 'Hz';
    fTbl.Properties.VariableNames{c} = ...
        [tbl.Properties.VariableNames{c},' frequency'];
end
fTbl.Properties.Description = [tbl.Properties.Description,' frequency'];

ATbl = tbl; ATbl.Variables = A_block; 
varunits = tbl.Properties.VariableUnits;
for c = 1:width(ATbl)
    if length(varunits) >= c
        ATbl.Properties.VariableUnits{c} = varunits{c};
    end
    ATbl.Properties.VariableNames{c} = ...
        [tbl.Properties.VariableNames{c},' amplitude'];
end
fTbl.Properties.Description = [tbl.Properties.Description,' frequency'];

end