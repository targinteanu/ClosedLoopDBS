function [rawD, fltD, forD, timeBuffs] = ...
    InitializeRecording_cbmex(buffSize, filtorder, forecastwin, ...
    selRaw, selRaw2Flt, selRaw2For, selFlt2For)
% Initialize the multichannel data structures used in ClosedLoopDBS using
% BlackRock hardware with the cbmex function.
% Return structures for raw, filtered, and forecast data, as well as timing
% buffers.

IndShiftFIR = ceil(filtorder/2); % samples
selFor = sort(unique([selRaw2For, selFlt2For]));

connect_cbmex(); 
pause(1);

%% defining the data structures
% much of this is generic; can it be replaced with a non-cbmex-exclusive
% function that gets called?
try
    [rawH, rawT, rawB, rawN] = initRawData_cbmex(selRaw, buffSize);
catch ME
    warning(['Error on first attempt: ',ME.message]);
    pause(1);
    [rawH, rawT, rawB, rawN] = initRawData_cbmex(selRaw, buffSize);
end
fltH = initFilteredData(rawH, IndShiftFIR); 
[forH, forT, forB] = initForecastData(fltH, forecastwin);

rawD = [rawN; rawH; rawT; rawB]; 
forD = [rawN; forH; forT; forB];

fltT = cellfun(@(D) [1,0].*D, rawT, 'UniformOutput',false);
fltB = cellfun(@(D) [1,0].*D, rawB, 'UniformOutput',false);
fltD = [rawN; fltH; fltT; fltB];

%% timing buffers
% this is generic; can it be replaced with a non-cbmex-exclusive
% function that gets called?
timeBuffs = cell(size(rawH));
for ch = 1:length(rawH)
    t_ch = rawH{ch}(:,1);
    dt_ch = 1/rawN{ch}.SampleRate;
    for it = 2:length(t_ch)
        if isnan(t_ch(it))
            t_ch(it) = t_ch(it-1) + dt_ch;
        end
    end
    timeBuffs{ch} = t_ch;% + t0;
end

%% resizing buffers 
% should this be placed within the respective init scripts instead? 
if ~isempty(selRaw2Flt)
    fltD = fltD(:,selRaw2Flt); 
end
if ~isempty(selFor)
    forD = forD(:,selFor);
end

end