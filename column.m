function varargout = column(varargin)
% COLUMN MATLAB code for column.fig
%      COLUMN, by itself, creates a new COLUMN or raises the existing
%      singleton*.
%
%      H = COLUMN returns the handle to a new COLUMN or the handle to
%      the existing singleton*.
%
%      COLUMN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COLUMN.M with the given input arguments.
%
%      COLUMN('Property','Value',...) creates a new COLUMN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before column_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to column_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help column

% Last Modified by GUIDE v2.5 15-Aug-2017 23:28:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @column_OpeningFcn, ...
                   'gui_OutputFcn',  @column_OutputFcn, ...
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


% --- Executes just before column is made visible.
function column_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to column (see VARARGIN)

% Choose default command line output for column
handles.output = hObject;
global clm

handles.n_tray.String=clm.n_tray;
handles.f_tray.String=clm.f_tray;
handles.rr.String=clm.rr;
handles.dp.String=clm.dp;
handles.d_rate.String=clm.distil_rate;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes column wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = column_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function n_tray_Callback(hObject, eventdata, handles)
% hObject    handle to n_tray (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of n_tray as text
%        str2double(get(hObject,'String')) returns contents of n_tray as a double


% --- Executes during object creation, after setting all properties.
function n_tray_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n_tray (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function f_tray_Callback(hObject, eventdata, handles)
% hObject    handle to f_tray (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of f_tray as text
%        str2double(get(hObject,'String')) returns contents of f_tray as a double


% --- Executes during object creation, after setting all properties.
function f_tray_CreateFcn(hObject, eventdata, handles)
% hObject    handle to f_tray (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rr_Callback(hObject, eventdata, handles)
% hObject    handle to rr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rr as text
%        str2double(get(hObject,'String')) returns contents of rr as a double


% --- Executes during object creation, after setting all properties.
function rr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dp_Callback(hObject, eventdata, handles)
% hObject    handle to dp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dp as text
%        str2double(get(hObject,'String')) returns contents of dp as a double


% --- Executes during object creation, after setting all properties.
function dp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
global clm


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global clm

clm.n_tray=str2double(handles.n_tray.String);
clm.f_tray=str2double(handles.f_tray.String);
clm.rr=str2double(handles.rr.String);
clm.dp=str2double(handles.dp.String);
clm.distil_rate=str2double(handles.d_rate.String);
% Hint: delete(hObject) closes the figure
delete(hObject);



function d_rate_Callback(hObject, eventdata, handles)
% hObject    handle to d_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of d_rate as text
%        str2double(get(hObject,'String')) returns contents of d_rate as a double


% --- Executes during object creation, after setting all properties.
function d_rate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to d_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
