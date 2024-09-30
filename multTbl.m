function tbl = multTbl(scalar, tbl)
% works like .*, but on a timetable
% might not be necessary in newer versions of matlab
tbl.Variables = scalar .* tbl.Variables;
end