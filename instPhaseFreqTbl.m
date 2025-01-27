function [phiTbl, fTbl] = instPhaseFreqTbl(tbl)
% 
% get instantaneous phase and frequency at each timepoint in a timetable 
% 

fs = tbl.Properties.SampleRate; 
[phi_block, f_block] = instPhaseFreq(tbl.Variables, fs);

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

end