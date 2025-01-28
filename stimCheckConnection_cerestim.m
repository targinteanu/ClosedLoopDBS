function [stimulator, dl] = stimCheckConnection_cerestim()
% check cerestim connection before proceeding
    stimulator = cerestim96();
    dl = stimulator.scanForDevices();
    if isempty(dl)
        warnmsg = ['Difficulty connecting to CereStim. ' ...
            'It is recommended to turn the device off and on ' ...
            'again before proceeding.'];
        resp = questdlg(warnmsg, 'CereStim connection issue', ...
            'Do Not Proceed', 'Proceed Anyway', 'Do Not Proceed');
        if ~strcmp(resp, 'Proceed Anyway')
            error('Try again after turning the CereStim off and on.')
        end
    end
end