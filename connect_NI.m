function DAQobj = connect_NI(device)

if nargin < 1
    device = 'USB-6212 (BNC)';
end

% consider adding daqreset

DL = daqlist("ni");

if isempty(DL)
    daqreset;
    error('No NI DAQs found. Check connection, or try resetting DAQ toolbox (daqreset) and try again.')
end

% find the desired device 
devicefound = find(strcmp(DL.Model, device));
if isempty(devicefound)
    warning(['No device ',device,' found. Will attempt to connect to a different ' ...
        'device. If this is not expected, check connection, or try resetting DAQ ' ...
        'toolbox (daqreset) and try again.']);
    devicefound = 1;
elseif length(devicefound) > 1
    warning(['Multiple NI DAQs found matching search criteria ''',device,'''.' ...
        'Will attempt to use the first available.']);
    devicefound = devicefound(1);
end

% find device ID/name for connecting 
deviceName = DL.DeviceID(devicefound);
deviceInfo = DL.DeviceInfo(devicefound);
device = char(DL.Model(devicefound));
disp(['Attempting to connect to NI DAQ ''',device,'''...']);
DAQobj = daq("ni");

% find analog in chan names and details 
subsystemtypes = arrayfun(@class, deviceInfo.Subsystems, 'UniformOutput',false);
subsys = find(contains(lower(subsystemtypes), 'analoginputinfo'));
if isempty(subsys)
    error('No analog input channels found.')
end
subsys = deviceInfo.Subsystems(subsys);
for sys = subsys
    for ch = string(sys.ChannelNames)'
        addinput(DAQobj, deviceName, ch, sys.DefaultMeasurementType);
    end
end

% start continuous data stream/buffer 
start(DAQobj, "continuous");

% success: 
disp('...Connected to:');
disp(DAQobj);

end