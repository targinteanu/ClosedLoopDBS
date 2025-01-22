function srlLastMsg = test_bgSrl(DQ)

thisportname = FindMySerialPort();
receiverSerial = serialport(thisportname, 9600);
configureCallback(receiverSerial,"terminator", ...
    @(hSrl,evt)CharSerialCallbackReceiver(hSrl,evt,DQ) );

%{
while ~strcmp(receiverSerial, 'stop')
    pause(.5);
end
%}
pause(20);

%srlLastMsg = receiverSerial.UserData;
srlLastMsg = 'goodbye';

function CharSerialCallbackReceiver(hSrlPrt, evt, DQ)
receivedData = readline(hSrlPrt); 
receivedData = char(receivedData); 
dspMsg = ['Serial Port ',char(hSrlPrt.Port),' Received Message: ',receivedData];
disp(dspMsg);
send(DQ, dspMsg);
hSrlPrt.UserData = dspMsg;
end

end