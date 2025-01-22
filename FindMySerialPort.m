function thisport = FindMySerialPort()

thiscomputer = getenv("COMPUTERNAME");
thiscomputertype = computer; 

if contains(thiscomputertype, 'MAC') || contains(thiscomputertype, 'mac')
    % Toren Arginteanu - Macbook 
    thisport = "/dev/tty.usbserial-A10K10JC";
else
    
if strcmp(thiscomputer, 'WSA_DELL_LAPTOP')
    % Anderson Lab - Dell Latitude Laptop 
    thisport = "COM4";
elseif strcmp(thiscomputer, 'LAPTOP-Q2F2EVC4')
    % Anderson Lab - Lenovo Thinkpad paradigm laptop 
    thisport = "COM6";
elseif strcmp(thiscomputer, 'LAPTOP-T2VLTA0O')
    % Kelly Mills - movement disorders Lenovo laptop
    thisport = "COM3";
elseif strcmp(thiscomputer, 'DESKTOP-SDRL657')
    % Yousef Salimpour - Surface Book 
    thisport = "COM3";
elseif strcmp(thiscomputer, 'DESKTOP-PF4UFTT')
    % Toren Arginteanu - Surface Book
    %thisport = '';    
    thisport = 'COM9';
else
    thisport = '';
end

end

end