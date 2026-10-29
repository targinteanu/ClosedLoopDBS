function T2 = myRetime(T1, fs, constval, tol)
% Retime the timetable T1 to a constant sampling rate fs. Any point within
% tol (default 2) samples of a value in the original table will take on
% that value; otherwise, it will be substituted with a constant value
% (default 0).

if nargin < 3
    constval = int16(0);
end
if nargin < 4
    tol = 2;
end

[~,sortind] = sort(T1.Time);
T1 = T1(sortind,:);

T2 = [];
dt = seconds(diff(T1.Time));
iinsert = dt > tol/fs; iinsert = find(iinsert);
i1 = 1;
for i2 = iinsert'
    T2 = [T2; T1(i1:i2, :)]; % pre-insert
    i1 = i2+1;
    ti = T1.Time(i2):seconds(1/fs):T1.Time(i1);
    ti = ti(2:(end-1));
    Di = repmat(constval, length(ti), width(T1));
    Ti = array2timetable(Di,"RowTimes",ti,"VariableNames",T1.Properties.VariableNames);
    T2 = [T2; Ti]; % insert 
end
T2 = [T2; T1(i1:end, :)];

end