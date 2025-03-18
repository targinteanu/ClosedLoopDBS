% copied from Get_channel_Data_test.m

MAC_addr = 'C8:DF:84:F8:9F:E2';

sampling_rate  = 44000;  % in Hz, 44000 for SPK
BufferSizemSec = 20000; % size of recording buffer in milliseconds (not sure what it does, does not change size of data packet)
matlab_pause   = 0.00001;   % how long to wait between data requests (in seconds) \\up to how much data you want to see (should take proccessing time of your code in calculation, if proccessing time is bigger than the pause  )
time_to_display= 1;     % seconds
mVolt_to_display = 500;
First_channel_ID=10128;
Channels_Count=63; % was 128
Channels=[];

for i=0:(Channels_Count-1)
    Channels(i+1) = First_channel_ID+i;
end
Channel_Gain = 20;
Bit_resolution = 2500000/(2^16*Channel_Gain);

%END of USER DEFINED VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AO_DefaultStartConnection(MAC_addr)
%%
Result=1;
for i=1:Channels_Count
    while Result~=0
        Result=AO_AddBufferingChannel(Channels(i),BufferSizemSec); %ChannelID,BufferSizemSec
    end
    Result=1;
end
Result=1;
AO_ClearChannelData();
pause(1);
pData=[];
for i=1:Channels_Count
    %while Result~=0 || isnan(DC)
        [Result1(i),pData(i,:),DC(i)]=AO_GetChannelData(Channels(i)); %ChannelID
    %end
    %Result=1;
end