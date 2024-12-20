thisportname = FindMySerialPort();
receiverSerial = serialport(thisportname, 9600);
ud = struct('ReceivedData', '', ...
            'ParadigmPhase', 'Stopped');
receiverSerial.UserData = ud;
configureCallback(receiverSerial,"terminator",...
    @(hsrl,evt)CharSerialCallbackReceiver_PD_v0(hsrl,evt, ...
                    handles.textSrl, handles.txt_Status)); 
pause(1);

fut = parfeval(@bg_Serial, 1, receiverSerial);
pause(1);
fetchOutputs(fut)

delete(receiverSerial)

function out = bg_Serial(srl)
out = getfield(srl, 'UserData');
end