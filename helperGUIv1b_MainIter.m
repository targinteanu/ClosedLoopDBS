function handles = helperGUIv1b_MainIter(handles)
% wrap and run iterReadBrain on the main CPU process for GUI version 1b

selRaw2Flt = handles.selInds.selRaw2Flt;
selFlt2For = handles.selInds.selFlt2For;
selFor2Art = handles.selInds.selFor2Art;
selRaw2Art = handles.selInds.selRaw2Art;
selRaw2For = handles.selInds.selRaw2For;

rawData = handles.recDataStructs.rawD;
artRemData = handles.recDataStructs.artD;
fltData = handles.recDataStructs.fltD;
forData = handles.recDataStructs.forD;
timeBuffs = handles.recDataStructs.timeBuffs;
initTic = handles.recDataStructs.initTic;
forBuffs = handles.recDataStructs.forBuffs;

% enable defined filtering/forecasting funcs only if ready
if handles.FilterSetUp
    filtfun = @filtFun;
else
    filtfun = [];
end
if handles.MdlSetUp
    forefun = @foreFun;
else
    forefun = [];
end
if handles.check_artifact.Value && handles.MdlSetUp
    artremfun = @artRemFun;
else
    artremfun = [];
end

[...
timeBuffs, rawData, ...
artRemData, handles.artRemArgs, ...
fltData, handles.filtArgs, ...
forBuffs, forData, handles.foreArgs] = ...
iterReadBrain(...
    timeBuffs, rawData, @() handles.HardwareFuncs.GetNewRawData(handles.allChannelIDs), ...
    selRaw2Art, selFor2Art, selRaw2Flt, selRaw2For, selFlt2For, ...
    artRemData, artremfun, handles.artRemArgs, ...
    fltData, filtfun, handles.filtArgs, ...
    forBuffs, forData, forefun, handles.foreArgs);

handles.recDataStructs.rawD = rawData;
handles.recDataStructs.artD = artRemData;
handles.recDataStructs.fltD = fltData;
handles.recDataStructs.forD = forData;
handles.recDataStructs.timeBuffs = timeBuffs;
handles.recDataStructs.forBuffs = forBuffs;


% --- PhaseDetect helpers ---

function [artRemTails, artRemArgs] = ...
        artRemFun(artRemArgs, rawTails, forTails)
% forTails must be the forecast data starting at the same time as rawTails
artRemTails = cell(size(rawTails));
FsArt = artRemArgs.SampleRates;
StimDur = artRemArgs.StimDur; % seconds to remove
StimLen = ceil(StimDur.*FsArt); % #samples to remove 
StimTimesTail = artRemArgs.StimTimes; 
for ch_art = 1:size(rawTails,2)
    tXfor = forTails{ch_art};
    tX = rawTails{ch_art}; % [time, data]
    stimtimes = StimTimesTail{ch_art}; % time to stim (sec)
    stiminds = round(stimtimes * FsArt(ch_art));
    stiminds = stiminds(stiminds > 0);
    for i1 = stiminds'
        i2 = i1 + StimLen(ch_art);
        if i2 > height(tX)
            % artifact will carry over into the next packet 
            nO = i2 - height(tX);
            artRemArgs.nOverlap(ch_art) = nO;
            tX = [tX; nan(nO,2)]; 
        else
            % artifact limited to this packet 
            artRemArgs.nOverlap(ch_art) = 0;
        end
        i2 = min(i2, height(tX)); % why is this necessary??
        tX(i1:i2,2) = tXfor(i1:i2,2);
    end
    artRemTails{ch_art} = tX;
end
end

function [foreTails, foreBuffsAdd, foreArgs] = foreFun(foreArgs, inData)
% inData should just be the (filtered) channel(s) of interest
% foreTails will be the forecast tails starting at the current time minus
% any filtering delay
if isempty(foreArgs)
    error('Attempted forecast function with empty arguments.')
end
K = foreArgs.K; k = foreArgs.k;
ARmdls = foreArgs.ARmdls;
fs = foreArgs.SampleRates;
Ts = foreArgs.TimeShift;
Fco = foreArgs.FreqRange;
phis = foreArgs.PhaseOfInterest;
TstimLag = foreArgs.StimulatorLagTime;
foreTails = cell(size(inData)); foreBuffsAdd = foreTails; 
for ch_fore = 1:size(inData,2)
    armdl = ARmdls{ch_fore};
    FT = myFastForecastAR(armdl, inData{ch_fore}(:,2), K);
    foreTails{ch_fore} = [nan(height(FT),1), FT];
    foreTails{ch_fore}(1,1) = foreArgs.TimeStart(ch_fore);

    FT = FT(1:k,:); % use limited duration for hilbert padding
    %[t2,i2,phi_inst,f_inst] = blockPDS(...
    t2 = blockPDS(...
        inData{ch_fore}(:,2), FT, fs(ch_fore), phis, ...
        Ts(ch_fore)+TstimLag, Fco(1), Fco(2));
    t2 = t2-Ts(ch_fore); % [t2peak, t2trough]
    % t2 = max(t2,0); % Should this be necessary? Will this cause problems?
    t2(t2 < 0) = nan;
    foreBuffsAdd{ch_fore} = t2; 
end
end

function [fltTails, fltArgs] = filtFun(fltArgs, rawTails)
filtObj = fltArgs.fltObj;
filtInit = fltArgs.fltInit;
filtFin = cell(size(filtInit));
fltTails = cell(size(rawTails));
for ch_flt = 1:size(rawTails,2)
    a = filtObj{1,ch_flt}; b = filtObj{2,ch_flt};
    [FT,filtFin{ch_flt}] = filter(b,a,rawTails{ch_flt}(:,2:end),filtInit{ch_flt});
    fltTails{ch_flt} = [rawTails{ch_flt}(:,1) - fltArgs.TimeShift(ch_flt), FT];
end
fltArgs.fltInit = filtFin;
end


end