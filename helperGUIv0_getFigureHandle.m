function currentHandle = helperGUIv0_getFigureHandle(currentHandle)
% For when the handle to the main GUI figure is needed, but some other
% object (e.g. a uicontrol that may be nested in a pannel) calls the
% function

while ~strcmp(class(currentHandle), 'matlab.ui.Figure')

    if isempty(currentHandle) || strcmp(class(currentHandle), 'matlab.ui.Root')
        error('Could not find figure handle.')
    end

    currentHandle = currentHandle.Parent;

end

end