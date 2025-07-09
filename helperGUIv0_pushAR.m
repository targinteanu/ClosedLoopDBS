function handles = helperGUIv0_pushAR(handles, PDSwin, n)

    handles.PDSwin2 = ceil(.8*PDSwin); % set model forecast duration

    % catch mistakes 
    if PDSwin > handles.bufferSize
        error('Phase estimation window cannot be larger than display window.')
    end
    if handles.PDSwin1 <= handles.IndShiftFIR
        error('Forecast length does not overcome filter delay.')
    end

    y = handles.filtDataBuffer; 
    L = min(length(y), 3*PDSwin) - 1;
    y = y((end-L):end);
    y = iddata(y,[],1/handles.fSample);
    ARmdl = ar(y,n,'yw');
    
    handles.Mdl = ARmdl; 
    handles.MdlSetUp = true;

end