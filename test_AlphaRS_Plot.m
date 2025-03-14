% copied from plotting_channels.m

%AO_DefaultStartConnection('bc:6a:29:5f:da:15');
sampling_rate  = 22000;  % in Hz, 44000 for SPK
BufferSizemSec = 20000; % size of recording buffer in milliseconds (not sure what it does, does not change size of data packet)
matlab_pause   = 0.00001;   % how long to wait between data requests (in seconds) \\up to how much data you want to see (should take proccessing time of your code in calculation, if proccessing time is bigger than the pause  )
time_to_display= 1.5;     % seconds
mVolt_to_display = 500;
First_channel_ID=10256;
Channels_Count=94;
Channels_Names=[];
Channels=[];

for i=0:Channels_Count-1
    Channels(i+1) = First_channel_ID+i;
%     if i<=8
%         Channels_Names=[Channels_Names;strcat('ECOG_HF0',num2str(i+1))];
%     else
%         Channels_Names=[Channels_Names;strcat('ECOG_HF',num2str(i+1))];
%     end
end
%Channels_Names = ['LFP 01';'LFP 02';'LFP 03';'LFP 04';'LFP 05'];
Channel_Gain = 20;
Bit_resolution = 2500000/(2^16*Channel_Gain);

%END of USER DEFINED VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:Channels_Count
    R(i)=AO_AddBufferingChannel(Channels(i),BufferSizemSec); %ChannelID,BufferSizemSec
end

AO_ClearChannelData(); %we clear the old data so we can have the new data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % prepare the plot
my_time    = 0:1/sampling_rate:time_to_display - 1/sampling_rate;
% 
y  = -inf*ones(length(my_time),Channels_Count);
h=figure(1);
for m = 1:64
    if Channels_Count<=64
        set(0, 'CurrentFigure', h);
        subplot(ceil(Channels_Count/8),8,m)
        xlim([0,time_to_display])
     %  xlabel('[s]')
        ylim([-mVolt_to_display,mVolt_to_display])
    %     ylabel('[mV]')
        Data_Line(m) = line(my_time,y(:,m),...
                        'markersize',5,...
                        'linestyle','-');
    else
        set(0, 'CurrentFigure', h);
        subplot(8,8,m)
        xlim([0,time_to_display])
     %  xlabel('[s]')
        ylim([-mVolt_to_display,mVolt_to_display])
    %     ylabel('[mV]')
        Data_Line(m) = line(my_time,y(:,m),...
                        'markersize',5,...
                        'linestyle','-');
    end
end
if Channels_Count>64
    h2=figure(2);
%     set(0, 'CurrentFigure', h2);
    for m = 65:Channels_Count
        set(0, 'CurrentFigure', h2);
        subplot(ceil((Channels_Count-64)/8),8,m-64)
        xlim([0,time_to_display])
        %  xlabel('[s]')
        ylim([-mVolt_to_display,mVolt_to_display])
        %        ylabel('[mV]')

        Data_Line(m) = line(my_time,y(:,m),...
            'markersize',5,...
            'linestyle','-');
        %    Vertical_Line(m) =  line([inf,inf],[-mVolt_to_display,mVolt_to_display]);

    end
end


shg; % Show most recent graph

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r = 0;
tic
Data=[];
Old_TS=0;
Old_DC=0;
while(toc < 5) %verifing proccessing time is lower than 1 sec.
tic
%     pause(matlab_pause);
    
    while size(Data,1)<=time_to_display*sampling_rate
        DataCapture=0;
        Result=0;
        pause(matlab_pause);
        
       [Result,pData,DataCapture,TS_FirstSample] = AO_GetAlignedData(Channels);        
        if DataCapture==0 || Result~=0
            continue;
        end
        if Old_TS~=0
            Diff=double(TS_FirstSample)-(Old_DC/Channels_Count);
            if Diff>Old_TS
                Miss=strcat('Missing data =',num2str(double(TS_FirstSample)-(Old_DC/Channels_Count)-double(Old_TS)))
            end
        end
        Old_TS=TS_FirstSample;
        Old_DC=DataCapture;
        mt_Channels_Data = reshape(pData(1:DataCapture),...
                           [DataCapture/Channels_Count,Channels_Count]);
        Data=[Data;mt_Channels_Data];
    end
    my_samples = (Bit_resolution)*Data(1:time_to_display*sampling_rate,:);
    Data=Data(time_to_display*sampling_rate+1:end,:);
%        my_samples = (Bit_resolution)*mt_Channels_Data(:,:);
   
%     for i = 1:size(my_samples,1)       
% 
%         r = r+1;
%         s = rem(r-1,size(y,1))+1;
%         y(s,:) = my_samples(i,:);
%     
%     end
    if ~ishandle(h) || (~ishandle(h2) && Channels_Count>65)
        return
    end
    figure(1);
    for m = 1:Channels_Count
        if Channels_Count<=64
            subplot(ceil(Channels_Count/8),8,m);
            set(Data_Line(m),...
              'xdata',my_time,...
              'ydata',my_samples(:,m));
        else
            if m<=64
                set(0, 'CurrentFigure', h);
                subplot(8,8,m);
                set(Data_Line(m),...
                  'xdata',my_time,...
                  'ydata',my_samples(:,m));
            elseif m==65
                set(0, 'CurrentFigure', h2);
                subplot(ceil((Channels_Count-64)/8),8,m-64);
                set(Data_Line(m),...
                  'xdata',my_time,...
                  'ydata',my_samples(:,m));
            else
                set(0, 'CurrentFigure', h2);
                subplot(ceil((Channels_Count-64)/8),8,m-64);
                set(Data_Line(m),...
                  'xdata',my_time,...
                  'ydata',my_samples(:,m));
            end
        end
%         subplot(ceil(Channels_Count/8),8,m);
%         set(Data_Line(m),...
%           'xdata',my_time,...
%           'ydata',my_samples(:,m));
%         set(Vertical_Line(m),...
%           'xdata',[my_time(s),my_time(s)]);
    
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%