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
if handles.StimActive
    noRecentStim = ...
        handles.stimLastTime + handles.artRemArgs.StimDur - handles.artRemArgs.ArtifactStartBefore ...
        < handles.lastSampleProcTime;
else
    noRecentStim = true;
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
    forBuffs, forData, forefun, handles.foreArgs, ...
    @(args,data) ARmdlUpdate(args,data,noRecentStim));

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
artRemArgs.nOverlap = zeros(1, size(rawTails,2));
for ch_art = 1:size(rawTails,2)
    tXfor = forTails{ch_art};
    tX = rawTails{ch_art}; % [time, data]
    [tN, tProc] = get_tProc(tX);
    [tForN, tForProc] = get_tProc(tXfor);
    tDiff = tForProc - tProc; iDiff = round(tDiff * FsArt(ch_art)); 
    iDiff = iDiff + tForN - tN;
    stimtimes = StimTimesTail{ch_art}; % time to stim FROM STARTUP (sec)
    stimtimes = stimtimes - tProc; % from last proc
    stimtimes = stimtimes - artRemArgs.ArtifactStartBefore;
    stiminds = floor(stimtimes * FsArt(ch_art));
    stiminds = stiminds + tN; % from tail start
    stiminds = stiminds(stiminds > 0);
    for i1 = stiminds'
        if ~isnan(i1)
        i2 = i1 + StimLen(ch_art) - 1;
        if i2 > height(tX)
            % artifact will carry over into the next packet 
            nO = i2 - height(tX);
            tX = [tX; nan(nO,2)]; 
        else
            % artifact limited to this packet 
            nO = 0;
        end
        artRemArgs.nOverlap(ch_art) = max(artRemArgs.nOverlap(ch_art), nO);
        if ~isnan(iDiff)
            i3 = i1 + iDiff; i4 = i2 + iDiff; 
            i3 = max(i3, 1); i3 = min(i3, height(tXfor));
            i4 = max(i4, 1); i4 = min(i4, height(tXfor));
            tXreplace = tXfor(i3:i4,2); 
            i2 = min(i2, i1+height(tXreplace)-1);
            tX(i1:i2,2) = tXreplace;
        end
        end
    end
    artRemTails{ch_art} = tX;
end
end

function [tN, tProc] = get_tProc(tX)
    % find the last proc time in a time-data buffer (in absolute seconds
    % from hardware startup) and the sample-offset from buffer start. 
    % Start with the assumption that the first sample has a valid t: 
    tN = 1; tProc = tX(tN,1); % this should not be NaN by def
    if isnan(tProc)
        % If above fails, find the last valid proc time explicitly:
        tN = find(~isnan(tX(:,1)));
        if ~isempty(tN)
            tN = tN(end); tProc = tX(tN,1);
        else
            tN = 1;
        end
    end
    tN = tN - 1; % tProc sample offset from tail start
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

function mdlArgs = ARmdlUpdate(mdlArgs, inData, noRecentStim)
    % inData should just be the (filtered) channel(s) of interest
    if noRecentStim
    for ch_upd = 1:size(inData, 2)
        ARupdated = false;
        if mdlArgs.ARlearnrate(ch_upd) > 0
            Mdl = mdlArgs.ARmdls{ch_upd};
            w = -Mdl(2:end)/Mdl(1);
            w = fliplr(w); 
            dataPast = inData{ch_upd}((end-length(w)):end, 2);
            E = dataPast(end,:); x = dataPast(1:(end-1),:);
            ypred = w*x; E = E - ypred;
            del = x*E;
            %del = del./(x'*x + eps); % normalize
            w = w + mdlArgs.ARlearnrate(ch_upd) * del';
            r = roots([1, -fliplr(w)]);
            if max(abs(r)) < 1 % ensure stability
                Mdl = [1, -fliplr(w)];
                ARupdated = true;
            end
            % simulate forward to determine if model blowing up
            if ARupdated
                dataFutu = myFastForecastAR(Mdl, dataPast, height(dataPast));
                if norm(dataFutu) <= 10*norm(dataPast)
                    mdlArgs.ARmdls{ch_upd} = Mdl;
                end
            end
        end
    end
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