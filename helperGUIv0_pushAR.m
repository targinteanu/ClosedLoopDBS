function handles = helperGUIv0_pushAR(handles, PDSwin, n, y)

if nargin < 4
    y = handles.filtDataBuffer; % only works for v0
end

    handles.PDSwin2 = ceil(.8*PDSwin); % set model forecast duration
    %handles.PDSwin2 = PDSwin; % minimize edge effects by putting point in middle

    % catch mistakes 
    if PDSwin > handles.bufferSize
        error('Phase estimation window cannot be larger than display window.')
    end
    if handles.PDSwin1 <= handles.IndShiftFIR
        error('Forecast length does not overcome filter delay.')
    end
    if handles.PDSwin2 < 2*handles.fSample/handles.locutoff
        error('Phase estimation window must allow at least two full cycles.')
    end

    L = min(length(y), 3*PDSwin) - 1;
    y = y((end-L):end);
    y = iddata(y,[],1/handles.fSample);
    ARmdl = ar(y,n,'yw');
    
    handles.Mdl = ARmdl.A; 
    handles.MdlSetUp = true;

end