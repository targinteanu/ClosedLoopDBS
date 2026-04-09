function events = parseNEVSerialCodes(nevFilePath, codesToKeep)
% parseNEVSerialCodes: Extract SerialDigitalIO events from a Blackrock .nev file
% Inputs:
%   nevFilePath  - path to .nev file
%   codesToKeep  - vector of serial codes to extract (e.g. [61 10 11 1:3 5 6 70 71 7 8 9])
% Outputs:
%   events       - table with columns:
%                   Time      : datetime
%                   Time_s    : seconds relative to first NEV timestamp
%                   EndTime   : optional end-time (datetime)
%                   EndTime_s : seconds relative to first NEV timestamp
%                   Code      : numeric serial code
%                   EventName : human-readable label
% Notes:
%   - Expects nev2table to return a table with columns Time, EventLabels, EventEnds
%   - Only extracts single numeric tokens from EventLabels
%   - Uses default code-to-label mapping for EventName

if nargin < 2 || isempty(codesToKeep)
    error('You must supply a vector of codes to keep, e.g. [61 10 11 1:3 5 6 70 71 7 8 9].');
end

% Default code-to-label mapping
codeLabelMap = containers.Map('KeyType','double','ValueType','char');
codeLabelMap(61) = 'start_experiment';
codeLabelMap(10) = 'break_in_experiment';
codeLabelMap(11) = 'fixation';
codeLabelMap(1)  = 'image_1';
codeLabelMap(2)  = 'image_2';
codeLabelMap(3)  = 'image_3';
codeLabelMap(5)  = 'encoding_mid';
codeLabelMap(6)  = 'end_encoding';
codeLabelMap(70) = 'stim_start';
codeLabelMap(71) = 'stim_end';
codeLabelMap(7)  = 'decoding_start';
codeLabelMap(8)  = 'decoding_end';
codeLabelMap(9)  = 'decoding_end';

% Load events table from nev2table
T = nev2table(nevFilePath);
% only select rows with "SerialDigitalIO" events
T = T(contains(T.EventLabels, 'SerialDigitalIO'), :);
T = sortrows(T, 'Time');
if isa(T, 'eventtable')
    timeVals = T.Time;
    rawLabels = T.EventLabels;
    endVals = T.EventEnds;
elseif istable(T) && all(ismember({'Time','EventLabels','EventEnds'}, T.Properties.VariableNames))
    timeVals = T.Time;
    rawLabels = T.EventLabels;
    endVals = T.EventEnds;
else
    disp('--- nev2table output diagnostic ---');
    disp(['Class: ' class(T)]);
    if istable(T) || istimetable(T)
        disp('Variable names:');
        disp(T.Properties.VariableNames);
    elseif isstruct(T)
        disp('Struct fields:');
        disp(fieldnames(T));
    else
        disp('Object properties:');
        try
            disp(properties(T));
        catch
            disp('Could not list properties.');
        end
        disp('Object fields:');
        try
            disp(fieldnames(T));
        catch
            disp('Could not list fields.');
        end
    end
    error('nev2table output must be a table or eventtable with Time, EventLabels, EventEnds.');
end

% Normalize EventLabels to cellstr
if isstring(rawLabels)
    rawLabels = cellstr(rawLabels);
elseif ischar(rawLabels)
    rawLabels = cellstr(rawLabels);
elseif iscategorical(rawLabels)
    rawLabels = cellstr(cellstr(rawLabels));
elseif ~iscell(rawLabels)
    rawLabels = cellfun(@num2str, rawLabels, 'uni', false);
end

N = numel(rawLabels);
codes = nan(N,1);
for k = 1:N
    s = rawLabels{k};
    m = regexp(s, '([0-9]+)', 'match');
    if ~isempty(m)
        codes(k) = str2double(m{end});
    end
end

% Filter requested codes
keepIdx = ismember(codes, codesToKeep);
if ~any(keepIdx)
    warning('No matching codes found in NEV table for the requested codes. Returning empty table.');
    events = table([], [], [], [], [], [], 'VariableNames', {'Time','Time_s','EndTime','EndTime_s','Code','EventName'});
    return;
end

keptTimes = timeVals(keepIdx);
keptEnds = endVals(keepIdx);
keptCodes = codes(keepIdx);

% Convert times to seconds relative to first NEV timestamp
if isdatetime(timeVals)
    t0 = timeVals(1);
    time_s = seconds(keptTimes - t0);
    end_s = nan(numel(keptEnds),1);
    for k = 1:numel(keptEnds)
        evEnd = keptEnds(k);
        if isempty(evEnd) || (isnumeric(evEnd) && all(isnan(evEnd)))
            end_s(k) = NaN;
        elseif isdatetime(evEnd)
            end_s(k) = seconds(evEnd - t0);
        elseif ischar(evEnd) || isstring(evEnd)
            try
                evdt = datetime(evEnd, 'TimeZone', 'UTC');
                end_s(k) = seconds(evdt - t0);
            catch
                end_s(k) = NaN;
            end
        else
            end_s(k) = NaN;
        end
    end
else
    error('Time column must be datetime.');
end

% event names
labelsOut = arrayfun(@(c) sprintf('code_%d', c), keptCodes, 'uni', false);
eventNames = cell(size(keptCodes));
for k = 1:numel(keptCodes)
    c = double(keptCodes(k));
    if isKey(codeLabelMap, c)
        eventNames{k} = codeLabelMap(c);
    else
        eventNames{k} = labelsOut{k};
    end
end

events = table(keptTimes(:), time_s(:), keptEnds(:), end_s(:), keptCodes(:), eventNames(:), ...
    'VariableNames', {'Time','Time_s','EndTime','EndTime_s','Code','EventName'});
end
