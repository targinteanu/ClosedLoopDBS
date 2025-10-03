%% setup

% find files and determine data fields 
fp = uigetdir;
filelist = dir(fp);
file1 = filelist(~[filelist.isdir]);
file1 = file1(arrayfun(@(f) f.bytes > 0, file1));
file1 = file1(arrayfun(@(f) f.name(1)~='.', file1));
file1 = file1(1); load(fullfile(file1.folder, file1.name));
channelNames = {Channel_ID_Name_Map.Name};
fileDataFields = {...
    'SF_DRIVE_CONF', ...
    ... 'SF_DRIVE_DOWN', ...
    'SF_DRIVE_SET_POS', ...
    ... 'SF_DRIVE_STOP', ...
    ... 'SF_DRIVE_UP', ...
    'SF_LEVEL' ...
    }; 
channelDataFields = {'BitResolution', 'Gain'};

% channel selection 
% For now, only select ANALOG_IN and LFP channels; ignore spike, RAW, and
% SEG. In the future, spike sorting should be performed later in this
% script instead of throwing out spikes. Not sure what to do with RAW/SEG
chincl = contains(channelNames, 'LFP') | ...
         contains(channelNames, 'ANALOG_IN');
channelNames = channelNames(chincl); 

% group channels by sampling rate 
Fs = cellfun(@getsamplerate, channelNames);
FsGrouped = unique(Fs);
channelNamesGrouped = cell(size(FsGrouped));
for grp = 1:length(FsGrouped)
    grpInds = Fs == FsGrouped(grp);
    channelNamesGrouped{grp} = channelNames(grpInds);
end
clear grp grpInds
Tbls0 = cell(size(FsGrouped));
Tbls = Tbls0;

% set file saving location
svloc = [fp,filesep,'Saved To Table',filesep,'Table Data ',...
    datestr(datetime, 'yyyy-mm-dd HH.MM.SS'),' ',...
    getFileVersion(mfilename("fullpath"))];
svN = 1;
sizethresh = 2e10; % size (bytes) at which to save and clear
pause(1)
mkdir(svloc); 

%% run1 - populate "horizontal" 

for f = filelist'
    clearvars -except ...
        Tbls Tbls0 channelNames channelDataFields fileDataFields ...
        f filelist svloc svN sizethresh ...
        Fs FsGrouped channelNamesGrouped
    if (~f.isdir) && (f.bytes > 0)
        fnfull = fullfile(f.folder, f.name); 
        [fp,fn,fe] = fileparts(fnfull);
        fni = find(fn == '.');
        fn1 = fn(1:(fni-1)); fn2 = fn((fni+1):end);
        fn1 = string(fn1); fn2 = string(fn2);
        if strcmpi(fe, '.mat')
            load(fnfull)

            %{
            % anticipate if the memory will become full; 
            % if so, save and clear 
            sz = whos('Tbls'); sz = sz.bytes;
            if (sz + f.bytes > sizethresh) && (f.bytes < sizethresh)
                save([svloc,filesep,'SavedTables',num2str(svN),'.mat'], "Tbls","channelNames");
                svN = svN+1;
                clear Tbls
                Tbls = Tbls0;
            end
            %}

            % look at channels 
            chNames = {Channel_ID_Name_Map.Name};
            chincl = contains(chNames, 'LFP') | ...
                     contains(chNames, 'ANALOG_IN');
            chNames = chNames(chincl);
            if ~strcmpwrapper(channelNames, chNames)
                warning(['On ',fn,': channel names do not match'])
            end
            FileData = varnames2struct(fileDataFields, '');
            for FSGRP = 1:length(channelNamesGrouped)
                TT = []; T1 = inf; T2 = -inf;
                chNamesGrp = channelNamesGrouped{FSGRP};
                for CHNAME = chNamesGrp
                chName = CHNAME{:};
                try
                t1 = eval([chName,'_TimeBegin']); % s
                t2 = eval([chName,'_TimeEnd']); % s
                fs = eval([chName,'_KHz'])*1000; % Hz
                T1 = min(T1, t1); T2 = max(T2, t2);
                t = t1:(1/fs):t2; 
                if exist([chName,'_BitResolution'],'var')
                    DataRes = eval([chName,'_BitResolution']);
                else
                    DataRes = 1;
                end
                if exist([chName,'_Gain'],'var')
                    DataGain = eval([chName,'_Gain']);
                else
                    DataGain = 1;
                end
                Data = eval(chName); % int
                Data = double(Data)*DataRes/DataGain;
                if length(t) ~= length(Data)
                    err = abs(length(Data) - length(t))/length(Data);
                    warning([chName,' mismatch in time and data length by ', num2str(100*err),'%'])
                    % treat times as more accurate than sample rate 
                    t = linspace(t1, t2, length(Data));
                end

                t = seconds(t);
                T = array2timetable(Data', 'RowTimes', t'); T.Properties.VariableNames{1} = chName;
                if isempty(TT)
                    TT = T;
                else
                    TT = synchronize(TT, T);
                end
                
                catch ME
                    warning(['On ',fn,' - ',chName,': ',ME.message])
                end
                clear chName t1 t2 t Data err FileName T
                end

                % mark table with file details 
                TT.Properties.Events = eventtable(seconds([T1;T2]), ...
                    "EventLabels", string(fn)+[" Start";" End"]);

                Tbls{FSGRP} = tblvertcat(Tbls{FSGRP}, TT);
                % if the memory is getting full, save and clear 
                sz = whos('Tbls'); sz = sz.bytes;
                if sz > sizethresh
                    save([svloc,filesep,'SavedTables',num2str(svN),'.mat'], "Tbls");
                    svN = svN+1;
                    clear Tbls
                    Tbls = Tbls0;
                end
            end
            clear FileData
        end
    end
end

% save if haven't yet
if sum(cellfun(@numel, Tbls))
    save([svloc,filesep,'SavedTables',num2str(svN),'.mat'], "Tbls");
    svN = svN+1;
end
clear Tbls
Tbls = Tbls0;

%% run2 - consolidate saved files "vertically" 
% if svN > 2
% there are multiple "Horizontal" files to consolidate

clearvars -except channelNames svloc svN sizethresh ...
    FsGrouped channelNamesGrouped
channelNames0 = channelNames; clear channelNames
filelist = dir(svloc);
filelist(~[filelist.isdir]);
filelist = filelist(arrayfun(@(f) f.bytes > 0, filelist));
filelist = filelist(contains({filelist.name}, 'SavedTables'));
svN = 1;

for FSGRP = 1:length(channelNamesGrouped)
    fs = FsGrouped(FSGRP);
    Tbl = [];
    for f = filelist'
        fnfull = fullfile(f.folder, f.name); 
        [fp,fn,fe] = fileparts(fnfull);
        if strcmpi(fe, '.mat')
            load(fnfull);
            T = Tbls{FSGRP}; clear Tbls;
            if isempty(Tbl)
                Tbl = T;
            else
                Tbl = tblvertcat(Tbl, T);
            end

            % if the memory is getting full, save and clear
            sz = whos('Tbl'); sz = sz.bytes;
            if sz > sizethresh
                save([svloc,filesep,'SavedTable',num2str(fs),'Hz',num2str(svN),'.mat'], "Tbl");
                svN = svN+1;
                clear Tbl
                Tbl = [];
            end
        end
    end

    % save if haven't yet 
    if numel(Tbl)
        save([svloc,filesep,'SavedTable',num2str(fs),'Hz',num2str(svN),'.mat'], "Tbl");
        svN = svN+1;
    end
end

% end

%% helper 

function fs = getsamplerate(chName)
fs = 0;
    try 
        fs = evalinwrapper([chName,'_KHz'])*1000; % Hz
    catch ME
        warning(['On ',chName,': ',ME.message])
    end
end

function S = varnames2struct(varnames, header)
if nargin < 2
    header = '';
end
vars = cellfun(@(vn) evalinwrapper([header,vn]), varnames, 'UniformOutput', false);
S = cell2struct(vars, varnames, 2);
end

function yn = strcmpwrapper(strs1, strs2)
yn = length(strs1) == length(strs2); 
if yn 
    yn = prod(strcmp(strs1, strs2));
end
end

function var = evalinwrapper(varname)
try 
    var = evalin('base',varname);
catch ME 
    if contains(ME.message, 'Unrecognized function or variable')
        var = nan;
    else
        rethrow(ME)
    end
end
end

function T = tblvertcat(T1, T2)
try
    T = [T1; T2];
catch ME
    if contains(ME.identifier, 'SizeMismatch')
        T1small = width(T1) < width(T2);
        if T1small
            Tsmall = T1; Tlarge = T2;
        else
            Tsmall = T2; Tlarge = T1;
        end
        smallvars = Tsmall.Properties.VariableNames;
        for VAR = Tlarge.Properties.VariableNames
            var = VAR{:};
            if ~sum(strcmp(var, smallvars))
                Tsmall.(var) = nan(height(Tsmall),1);
            end
        end
        if T1small
            T = [Tsmall; Tlarge];
        else
            T = [Tlarge; Tsmall];
        end
    else
        rethrow(ME)
    end
end
end