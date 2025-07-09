function ConsolidateSavedData(filepath)
% Collapse multiple "SaveFile" .mat files output from OnlineDisplay into a
% single "OnlineDisplaySavedData.mat" file 

if nargin < 1
    filepath = uigetdir;
end

varnames = {'SerialLog','PeakTime','TroughTime','StimTime',...
    'Phase3Time','Phase4Time','Phase5Time'};
% init temp vars 
for V = varnames
    v = V{:};
    eval([v,'_temp = [];'])
end

mydir = dir([filepath,filesep,'*SaveFile*.mat']); 
myfiles = mydir(~[mydir.isdir]);
cd0 = cd;
cd(filepath)

for myfile = myfiles'
    for V = varnames
        v = V{:};
        clearvars(v)
    end
    load(myfile.name)
    for V = 1:length(varnames)
        v = varnames{V};
        if exist(v) == 1
            % omit nan
            if V == 1 % SerialLog
                eval([v,' = ',v,'(~isnan([',v,'.TimeStamp]));'])
            else
                eval([v,' = ',v,'(~isnan(',v,'));'])
            end
            % store to temp
            eval([v,'_temp = [',v,'_temp; ',v,'];'])
        end
    end
end

% retrieve temp 
for V = varnames
    v = V{:};
    clearvars(v)
    eval([v,' = ',v,'_temp;'])
    clearvars([v,'_temp'])
end
clearvars v V myfile

% save to new filename 
fn = ['OnlineDisplaySavedData_',...
    datestr(datetime, 'yyyy-mm-dd_HH-MM-SS'),'.mat'];
pause(1)
save(fn)

% delete files that have been consolidated 
for myfile = myfiles
    delete(myfile.name)
end
cd(cd0)

end