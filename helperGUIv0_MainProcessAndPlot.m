function handles = helperGUIv0_MainProcessAndPlot(handles, newContinuousData, controllerFun)

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

            % use AR model to get some future data 
            StimulatorLagInd = round(handles.fSample*handles.StimulatorLagTime);
            dataPast = handles.filtDataBuffer; 
            dataPast = dataPast((end-handles.PDSwin1+1):end);
            dataFutu = myFastForecastAR(handles.Mdl, dataPast, handles.PDSwin1);
            dataFutu2 = dataFutu(1:handles.PDSwin2);
            if handles.check_polar.Value
            handles.predDataBuffer = OverwriteAndCycle(...
                handles.predDataBuffer, dataFutu, N);
            set(handles.h_predTrace,'YData',handles.predDataBuffer);
            end

            % find the time to next peak, trough and plot 
            bp = norm(dataPast,2)^2/numel(dataPast); % band power surrogate 
            M = handles.PDSwin1 - handles.IndShiftFIR;
            dataPk = false(M,1); 
            dataTr = dataPk; % dataSt = dataPk;  
            [t2,i2,phi_inst,f_inst] = ...
                blockPDS(dataPast,dataFutu2, handles.fSample, [0,pi], ...
                handles.TimeShiftFIR + handles.StimulatorLagTime, ...
                handles.locutoff, handles.hicutoff);
            t2 = t2 - handles.TimeShiftFIR - handles.StimulatorLagTime; 
            i2 = i2 - handles.IndShiftFIR;% - StimulatorLagInd;
            %t2 = max(t2,0); i2 = max(i2,1);
            %t2peak = t2(1); t2trou = t2(2);
            i2peak = i2(1); i2trou = i2(2);
            if i2peak > 0
                dataPk(i2peak) = true;
            end
            if i2trou > 0
                dataTr(i2trou) = true;
            end
            [handles.peakDataBuffer, oldPeak] = CombineAndCycle(...
                handles.peakDataBuffer, dataPk, N, M-StimulatorLagInd); 
            [handles.trouDataBuffer, oldTrou] = CombineAndCycle(...
                handles.trouDataBuffer, dataTr, N, M-StimulatorLagInd);
            set(handles.h_peakTrace,'YData',0*plotLogical(handles.peakDataBuffer));
            set(handles.h_trouTrace,'YData',0*plotLogical(handles.trouDataBuffer));

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
            tOldPeak = tOldBuffer(oldPeak); tOldTrou = tOldBuffer(oldTrou);
            [handles.pkStorage1, handles.pkP1,  handles.pkStorage2, p2] = ...
                cycleStorage(handles.pkStorage1, handles.pkP1, ...
                             handles.pkStorage2, tOldPeak);
            if ~handles.pkP1
                % storage full; save 
                PeakTime = handles.pkStorage1;
                svfn = [handles.SaveFileLoc,filesep,'SaveFile',num2str(handles.SaveFileN),'.mat'];
                disp(['Saving Peaks to ',svfn])
                save(svfn,'PeakTime');
                handles.SaveFileN = handles.SaveFileN + 1;
                handles.pkStorage1 = handles.pkStorage2; 
                handles.pkP1 = p2; 
                handles.pkStorage2 = nan(size(handles.pkStorage2));
            end
            [handles.trStorage1, handles.trP1,  handles.trStorage2, p2] = ...
                cycleStorage(handles.trStorage1, handles.trP1, ...
                             handles.trStorage2, tOldTrou);
            if ~handles.trP1
                % storage full; save 
                TroughTime = handles.trStorage1;
                svfn = [handles.SaveFileLoc,filesep,'SaveFile',num2str(handles.SaveFileN),'.mat'];
                disp(['Saving Troughs to ',svfn])
                save(svfn,'TroughTime');
                handles.SaveFileN = handles.SaveFileN + 1;
                handles.trStorage1 = handles.trStorage2; 
                handles.trP1 = p2; 
                handles.trStorage2 = nan(size(handles.trStorage2));
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
N = length(newData);
if N >= length(oldBuffer)
    newBuffer = newData(end-length(oldBuffer)+1:end);
else
    newBuffer = [oldBuffer(N+1:end); newData];
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
N = length(newData); 
if p1+N-1 > length(storage1)
    % storage 1 is now full 
    if N > length(storage2)
        warning('Data overloaded save buffer; some data may not be saved.')
        N = length(storage2);
        newData = newData(1:N);
    end
    p1 = 0; 
    storage2(1:N) = newData; 
    p2 = N+1;
else
    storage1(p1:(p1+N-1)) = newData;
    p1 = p1+N;
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