function varargout = test_ParallelGUIDE(varargin)
% TEST_PARALLELGUIDE MATLAB code for test_ParallelGUIDE.fig
%      TEST_PARALLELGUIDE, by itself, creates a new TEST_PARALLELGUIDE or raises the existing
%      singleton*.
%
%      H = TEST_PARALLELGUIDE returns the handle to a new TEST_PARALLELGUIDE or the handle to
%      the existing singleton*.
%
%      TEST_PARALLELGUIDE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TEST_PARALLELGUIDE.M with the given input arguments.
%
%      TEST_PARALLELGUIDE('Property','Value',...) creates a new TEST_PARALLELGUIDE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before test_ParallelGUIDE_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to test_ParallelGUIDE_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help test_ParallelGUIDE

% Last Modified by GUIDE v2.5 01-Nov-2024 17:48:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @test_ParallelGUIDE_OpeningFcn, ...
                   'gui_OutputFcn',  @test_ParallelGUIDE_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before test_ParallelGUIDE is made visible.
function test_ParallelGUIDE_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to test_ParallelGUIDE (see VARARGIN)

% Choose default command line output for test_ParallelGUIDE
handles.output = hObject;

handles.dataBuffer = zeros(1000,1);
handles.QLength = nan(1000,1);
handles.dataMean = 0; handles.dataStd = 1;

axes(handles.axes1);
handles.h_plt = plot(handles.dataBuffer, '.');
grid on; hold on; 
handles.h_QL = plot(handles.QLength, '.');

handles.pool = gcp('nocreate');
if isempty(handles.pool)
    handles.pool = parpool;
end

handles.dataQueue = parallel.pool.PollableDataQueue;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes test_ParallelGUIDE wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = test_ParallelGUIDE_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in StartButton.
function StartButton_Callback(hObject, eventdata, handles)
% hObject    handle to StartButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject, handles);
handles.f = parfeval(@bg_Data, 1, ...
    handles.dataQueue, true, handles.dataMean, handles.dataStd);
guidata(hObject, handles);
while strcmp(handles.f.State, 'running') && handles.checkEnable.Value
    handles = guidata(hObject);
    pause(.5);
    while handles.dataQueue.QueueLength > 0
        x = poll(handles.dataQueue);
        l = nan(size(x)); 
        if ~isempty(x)
            l(end) = handles.dataQueue.QueueLength;
            handles.QLength = bufferData(handles.QLength, l);
            handles.dataBuffer = bufferData(handles.dataBuffer, x);
            handles.h_plt.YData = handles.dataBuffer;
            handles.h_QL.YData = handles.QLength;
            guidata(hObject, handles);
        end
    end
end
guidata(hObject, handles);


% --- Executes on button press in StopButton.
function StopButton_Callback(hObject, eventdata, handles)
% hObject    handle to StopButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cancel(handles.f)
guidata(hObject, handles);


function MeanBox_Callback(hObject, eventdata, handles)
% hObject    handle to MeanBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MeanBox as text
%        str2double(get(hObject,'String')) returns contents of MeanBox as a double
handles.dataMean = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function MeanBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MeanBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function StdBox_Callback(hObject, eventdata, handles)
% hObject    handle to StdBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of StdBox as text
%        str2double(get(hObject,'String')) returns contents of StdBox as a double
handles.dataStd = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function StdBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StdBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function x = bg_Data(DQ, RunMainLoop, xMean, xStd)
while RunMainLoop
    pause(.01);
    x = xStd*randn + xMean;
    send(DQ,x);
end


% --- Executes on button press in checkEnable.
function checkEnable_Callback(hObject, eventdata, handles)
% hObject    handle to checkEnable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkEnable
