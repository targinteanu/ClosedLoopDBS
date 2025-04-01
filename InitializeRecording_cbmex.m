function [rawD, artD, fltD, forD, timeBuffs, startTic] = ...
    InitializeRecording_cbmex(buffSize, filtorder, forecastwin, ...
    selRawID, selRaw2Art, selRaw2Flt, selRaw2For, selFlt2For)
% Initialize the multichannel data structures used in ClosedLoopDBS using
% BlackRock hardware with the cbmex function.
% Return structures for raw, filtered, and forecast data, as well as timing
% buffers.
% 
% Buffer Size(s) [buffSize], Filter Order(s) [filtorder], and Forecast
% Window(s) [forecastwin] can be a single input that applies to all 
% channels, an array of inputs for each channel, or empty to bypass the
% respective operation where applicable. 
% selxxx2xxx is a horizontal array of selected column indexes, but selRawID
% is ID number(s) and can be passed as empty to select all.

IndShiftFIR = ceil(filtorder/2); % samples
selFor = [selRaw2For, selFlt2For];

%% defining the data structures
% much of this is generic; can it be replaced with a non-cbmex-exclusive
% function that gets called?

emptyOut = {[]; []; []; []};

try
    [rawH, rawT, rawB, rawN, startTic] = initRawData_cbmex(selRawID, buffSize);
catch ME
    warning(['Error on first attempt: ',ME.message]);
    pause(1);
    [rawH, rawT, rawB, rawN, startTic] = initRawData_cbmex(selRawID, buffSize);
end
rawD = [rawN; rawH; rawT; rawB]; 

if numel(selRaw2Art)
    artD = rawD(:,selRaw2Art);
else
    artD = emptyOut;
end

if numel(IndShiftFIR) && numel(selRaw2Flt)
    fltH = initFilteredData(rawH(selRaw2Flt), IndShiftFIR); 
    fltT = cellfun(@(D) [1,0].*D, rawT(selRaw2Flt), 'UniformOutput',false);
    fltB = cellfun(@(D) [1,0].*D(1:(end-IndShiftFIR),:), rawB(selRaw2Flt), 'UniformOutput',false);
    fltN = rawN(selRaw2Flt);
    fltD = [fltN; fltH; fltT; fltB];
else
    fltD = emptyOut;
end

if numel(forecastwin) && numel(selFor)
    [forH, forT, forB] = initForecastData(...
        [rawH(selRaw2For), fltH(selFlt2For)], forecastwin);
    forD = [[rawN(selRaw2For), fltN(selFlt2For)]; forH; forT; forB];
else
    forD = emptyOut;
end

%% timing buffers
% this is generic; can it be replaced with a non-cbmex-exclusive
% function that gets called?
timeBuffs = cell(size(rawH));
for ch = 1:width(rawH)
    t_ch = rawH{ch}(:,1);
    dt_ch = 1/rawN{ch}.SampleRate;
    for it = 2:length(t_ch)
        if isnan(t_ch(it))
            t_ch(it) = t_ch(it-1) + dt_ch;
        end
    end
    timeBuffs{ch} = t_ch;% + t0;
end

end