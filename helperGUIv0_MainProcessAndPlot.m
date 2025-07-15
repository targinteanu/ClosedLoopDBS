function [handles, phaseTraceHandles, phaseBuffers, phaseStorage] = ...
    helperGUIv0_MainProcessAndPlot(handles, newContinuousData, ...
    controllerFun, phaseTraceHandles, phaseBuffers, phaseStorage)

    N = length(newContinuousData);

    % update filtered data 
    if handles.FilterSetUp
        try
        [handles.filtDataBuffer, handles.filtCond] = FilterAndCycle(...
            handles.filtDataBuffer, newContinuousData, ...
            handles.BPF, handles.filtCond);
        set(handles.h_filtDataTrace,'YData',handles.filtDataBuffer)

        % update model-forecasted data 
        if handles.MdlSetUp
            try
            Mdl = handles.Mdl;

            % update AR coeffs using gradient descent 
            ARupdated = false;
            if handles.stimLastTime + handles.ArtifactDuration - handles.ArtifactStartBefore ...
                    < handles.lastSampleProcTime
            if handles.ARlearnrate > 0
                w = -Mdl(2:end)/Mdl(1);
                w = fliplr(w); 
                x = handles.filtDataBuffer((end-length(w)):(end-1));
                ypred = w*x; 
                E = handles.filtDataBuffer(end) - ypred; del = x*E;
                del = del./(x'*x + eps); % normalize 
                w = w + handles.ARlearnrate * del';
                r = roots([1, -fliplr(w)]);
                if max(abs(r)) < 1 % ensure stability
                    Mdl = [1, -fliplr(w)];
                    ARupdated = true;
                end
            end
            end

            % use AR model to get some future data 
            StimulatorLagInd = round(handles.fSample*handles.StimulatorLagTime);
            dataPast = handles.filtDataBuffer; 
            dataPast = dataPast((end-handles.PDSwin1+1):end);
            dataFutu = myFastForecastAR(Mdl, dataPast, handles.PDSwin1);
            if ARupdated
            if norm(dataFutu) > 10*norm(dataPast)
                % revert to prevent blowing up
                dataFutu = myFastForecastAR(handles.Mdl, dataPast, handles.PDSwin1);
            else
                handles.Mdl = Mdl;
            end
            end
            dataFutu2 = dataFutu(1:handles.PDSwin2);
            if handles.check_polar.Value
            handles.predDataBuffer = OverwriteAndCycle(...
                handles.predDataBuffer, dataFutu, N);
            set(handles.h_predTrace,'YData',handles.predDataBuffer);
            end

            % find the time to next peak, trough and plot 
            bp = norm(dataPast,2)^2/numel(dataPast); % band power surrogate 
            M = handles.PDSwin1 - handles.IndShiftFIR;
            dataPh = false(M,1); oldPh = cell(size(phaseBuffers));
            [t2,i2,phi_inst,f_inst] = ...
                blockPDS(dataPast,dataFutu2, handles.fSample, handles.PhaseOfInterest, ...
                handles.TimeShiftFIR + handles.StimulatorLagTime, ...
                handles.locutoff, handles.hicutoff);
            t2 = t2 - handles.TimeShiftFIR - handles.StimulatorLagTime; 
            i2 = i2 - handles.IndShiftFIR;% - StimulatorLagInd;
            %t2 = max(t2,0); i2 = max(i2,1);
            %t2peak = t2(1); t2trou = t2(2);
            %i2peak = i2(1); i2trou = i2(2);
            for iph = 1:min(length(i2), length(phaseBuffers))
                dataPh_ = dataPh; i2_ = i2(iph);
                if (~isnan(i2_)) && ( (i2_ > 0) && (i2_ <= M) )
                    dataPh_(i2_) = true;
                end
                [phaseBuffers{iph}, oldPh{iph}] = CombineAndCycle(...
                    phaseBuffers{iph}, dataPh_, N, M-StimulatorLagInd-1);
                set(phaseTraceHandles{iph},'YData',0*plotLogical(phaseBuffers{iph}));
            end

            % time of stimulus 
            if handles.stimNewTime > 0
            stimtimerel = handles.stimNewTime - handles.lastSampleProcTime; 
                % rel to 0 on screen
                % lastSampleProcTime should be time 0 on the screen
            stimind = round(stimtimerel*handles.fSample); % ind rel to END of buffer
            stimind = stimind + handles.bufferSize; % ind rel to START of buffer 
            handles.stimind = stimind - N;
            handles.stimNewTime = -inf;
            if stimind > 0
                handles.stimDataBuffer(stimind) = true;
            end
            else
                handles.stimind = -1;
            end

            [handles.stimDataBuffer, oldStim] = CombineAndCycle(...
                handles.stimDataBuffer, [], N, 0);
            set(handles.h_stimTrace,'YData',0*plotLogical(handles.stimDataBuffer));

            % queue stimulus pulse, if applicable 
            % ***** TO DO: can this be moved elsewhere to avoid delays?  
            Stim2Q = false;
            i2Q = controllerFun(handles.srl, handles, bp, handles.bpthresh);
            if i2Q
                Stim2Q = true;
                t2Q = t2(i2Q);
            end
            if Stim2Q && (t2Q >= 0)
                % possibly overwrite existing timer with new one
                t2Q = .001*floor(1000*t2Q); % round to nearest 1ms 
                if t2Q < (100/handles.locutoff + handles.TimeShiftFIR)
                    t2Qabs = t2Q + handles.lastSampleProcTime; % in NSP time "absolute"
                    Dt2Q = t2Qabs - handles.stimLastTime; 
                    if 1/Dt2Q <= handles.stimMaxFreq
                        stim2Q_proceed = true;
                        if strcmp(handles.QueuedStim.Running, 'on')
                            % last queued stim has not yet fired 
                            if t2Qabs > handles.QueuedStim.UserData
                                % new requested point is later than current timer
                                stim2Q_proceed = t2Q > handles.StimulatorLagTime; 
                                    % is there enough time to make a change
                                stim2Q_proceed = stim2Q_proceed && ...
                                    (t2Qabs - handles.QueuedStim.UserData) > handles.StimulatorLagTime; 
                                    % is the change outside margin of error
                                stim2Q_proceed = stim2Q_proceed && ...
                                    (t2Qabs - handles.QueuedStim.UserData) < 1/handles.hicutoff; 
                                    % is it trying to target the next cycle
                            end
                        end
                        if stim2Q_proceed
                            if strcmp(handles.QueuedStim.Running, 'on')
                                stop(handles.QueuedStim);
                            end
                            handles.QueuedStim.StartDelay = t2Q;
                            handles.QueuedStim.UserData = t2Qabs;
                            start(handles.QueuedStim);
                        end
                    end
                end
            else
                if strcmp(handles.QueuedStim.Running, 'on')
                    stop(handles.QueuedStim);
                end
            end

            % plot sine wave 
            if handles.check_polar.Value
            tSin = (1:handles.PDSwin1)/handles.fSample; tSin = tSin';
            sinMag = handles.filtDataBuffer(end) / cos(phi_inst);
            sinData = sinMag*cos(2*pi*f_inst*tSin + phi_inst);
            handles.sineDataBuffer = OverwriteAndCycle(...
                handles.sineDataBuffer, sinData, N);
            set(handles.h_sineTrace,'YData',handles.sineDataBuffer);
            end

            % evaluate accuracy of above --> polar histogram
            if handles.check_polar.Value
            phiPk = handles.peakDataBuffer; 
            phiTr = handles.trouDataBuffer;
            phiPk = phiPk(1:length(handles.filtDataBuffer));
            phiTr = phiTr(1:length(handles.filtDataBuffer));
            phi = instPhaseFreq(handles.filtDataBuffer, handles.fSample);
            phiPk = phi(phiPk); phiTr = phi(phiTr);
            set(handles.h_peakPhase,'Data',phiPk);
            set(handles.h_trouPhase,'Data',phiTr);
            end

            % store peaks and troughs that have been buffered out
                % N samples of new continuous data have come in 
                % most recent time stamp is lastSampleProcTime 
            tOldBuffer = ((1-N):0)/handles.fSample + handles.lastSampleProcTime - ...
                handles.bufferSize/handles.fSample;
            for iph = 1:length(oldPh)
                tOld = tOldBuffer(oldPh{iph});
                if numel(tOld)
                [phaseStorage{1,iph},phaseStorage{2,iph},phaseStorage{3,iph},p2] = ...
                    cycleStorage(phaseStorage{1,iph},phaseStorage{2,iph},phaseStorage{3,iph},...
                    tOld);
                if ~phaseStorage{2,iph}
                    % storage full; save
                    phname = handles.PhaseOfInterestName(iph);
                    phname = phname+"Time";
                    %assignin("base",phname,phaseStorage{1,iph});
                    %eval(phname+" = phaseStorage{1,iph};");
                    dataNeedsName = phaseStorage{1,iph}; dataName = phname;
                    svfn = [handles.SaveFileLoc,filesep,'SaveFile',num2str(handles.SaveFileN),'.mat'];
                    disp("Saving "+phname+" to "+svfn)
                    save(svfn,'dataNeedsName','dataName');
                    handles.SaveFileN = handles.SaveFileN + 1;
                    phaseStorage{1,iph} = phaseStorage{3,iph};
                    phaseStorage{2,iph} = p2;
                    phaseStorage{3,iph} = nan(size(phaseStorage{3,iph}));
                end
                end
            end

            catch ME2
                getReport(ME2)
                errordlg(ME2.message, 'Model Prediction Issue');
                handles.MdlSetUp = false;
                pause(.01);
            end
        end

        catch ME1
            getReport(ME1)
            errordlg(ME1.message, 'Filtering Issue');
            handles.FilterSetUp = false;
            pause(.01);
        end
    end

%% helpers 


function newBuffer = cycleBuffer(oldBuffer, newData)
L = length(newData);
if L >= length(oldBuffer)
    newBuffer = newData(end-length(oldBuffer)+1:end);
else
    newBuffer = [oldBuffer(L+1:end); newData];
end
end


function plotData = plotLogical(logData)
% take in a logical array and output an array that will plot true as 1 
% and will not plot false
plotData = double(logData); 
plotData(~logData) = nan; 
end


function [storage1, p1, storage2, p2] = ...
    cycleStorage(storage1, p1, storage2, newData)
L = length(newData); 
if p1+L-1 > length(storage1)
    % storage 1 is now full 
    if L > length(storage2)
        warning('Data overloaded save buffer; some data may not be saved.')
        L = length(storage2);
        newData = newData(1:L);
    end
    p1 = 0; 
    storage2(1:L) = newData; 
    p2 = L+1;
else
    storage1(p1:(p1+L-1)) = newData;
    p1 = p1+L;
    p2 = [];
end
end


function [newBuffer, lastBuffer] = CombineAndCycle(oldBuffer, newData, N, M)
% ?? does this still work when N is longer than length oldBuffer ??
% M = length(newData); 
newBuffer = false(size(oldBuffer)); 
newBuffer(1:(end-N)) = oldBuffer((N+1):end);
lastBuffer = oldBuffer; 
if N < length(lastBuffer)
    lastBuffer = lastBuffer(1:N);
end
newBuffer((end-M+1):end) = newData((end-M+1):end);
end


function newBuffer = OverwriteAndCycle(oldBuffer, newData, N)
% N = # of points of data that is actually new 
if N <= length(newData)
    if N >= length(oldBuffer)
        newBuffer = newData(end-length(oldBuffer)+1:end);
    else
        % buffer AND overwrite 
        L = length(newData)-N; % length to overwrite
        newBuffer = [oldBuffer((N+1):(end-L)); newData];
    end
else
    % there is no data to overwrite; in fact, there is not enough new data
    % nan-pad newData to length N and cycle buffer 
    newData = [newData; nan(N-length(newData),1)];
    newBuffer = cycleBuffer(oldBuffer, newData);
end
end


function [newFiltBuffer, filterFinalCond] = FilterAndCycle(...
    oldFiltBuffer, newUnfilt, filtobj, filterInitCond)
[newFilt,filterFinalCond] = filter(filtobj,1,newUnfilt,filterInitCond);
newFiltBuffer = cycleBuffer(oldFiltBuffer, newFilt);
end

end