fs = 1000; % Hz
[fn,fp] = uigetfile('*.csv');
tbl = readtable(fullfile(fp,fn));
t = tbl.dataTimestamp/fs; % s
TT = table2timetable(tbl, 'RowTimes', seconds(t));
y = TT(:,2); 
figure; plot(y, 'data')
ARmdl = ar(y, 10, 'yw'); 
ar_coeffs = -1*(ARmdl.A(2:end));
fractionalScaling = 13; 
formatted_ar_coeffs = int32(round(ar_coeffs*2^fractionalScaling))'