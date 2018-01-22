function varargout = stream(varargin)
% STREAM MATLAB code for stream.fig
%      STREAM, by itself, creates a new STREAM or raises the existing
%      singleton*.
%
%      H = STREAM returns the handle to a new STREAM or the handle to
%      the existing singleton*.
%
%      STREAM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STREAM.M with the given input arguments.
%
%      STREAM('Property','Value',...) creates a new STREAM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before stream_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to stream_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help stream

% Last Modified by GUIDE v2.5 14-Aug-2017 02:02:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @stream_OpeningFcn, ...
    'gui_OutputFcn',  @stream_OutputFcn, ...
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


% --- Executes just before stream is made visible.
function stream_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to stream (see VARARGIN)

% Choose default command line output for stream
handles.output = hObject;



% UIWAIT makes stream wait for user response (see UIRESUME)
% uiwait(handles.stream);

% hossein added
global feed D_v D_l B_l sel_comp comp results

handles.flash.Visible='off';
% {'T (K)','';
%     'P (Kpa)','';
% 'Molar Flow (kgmol/s)','';
% 'Molar fraction',''};
data={'T (C)','';
    'P (kpa)','';
    'Molar Flow rate','';
    'Liquid molar fraction','';
    'Vapor molar fraction',''};
handles.uitable1.Data=data;
handles.uitable1.ColumnName = {'Conditions','Value'};
handles.uitable1.ColumnWidth={120,120};
handles.uitable1.ColumnEditable=logical([0 1]);



handles.uitable2.Data={'','','',''};
handles.uitable2.ColumnName = {'Molar fractions','Total','Liquid','Vapor'};
handles.uitable2.ColumnWidth={90,80,80,80};
if strcmp(varargin{1},'feed')
    handles.uitable2.ColumnEditable=logical([0 1 0 0]);
else
    handles.uitable2.ColumnEditable=logical([0 0 0 0]);
end
if isempty(sel_comp)
    warndlg('Please select component.','!! Warning !!','modal')
else
    
    if strcmp(varargin{1},'feed')
        handles.flash.Visible='on'
        handles.stream.Name='feed';
        handles.uitable2.Data=feed.comp;
        handles.uitable1.Data=feed.data;

    elseif strcmp(varargin{1},'distil_v')
        handles.stream.Name='distil_vapor';
        handles.uitable2.Data=D_v.comp;
    elseif strcmp(varargin{1},'distil_l')
        handles.stream.Name='distil_liquid';
        handles.uitable1.Data=data;
        handles.uitable1.Data(1,2)=cellstr(num2str(results.Tjk(1,end)-273));
        handles.uitable1.Data(2,2)=cellstr(num2str(results.Pj(1)/1e3));
        handles.uitable1.Data(3,2)=cellstr(num2str(results.Uj(1)));
        handles.uitable1.Data(4,2)=cellstr(num2str(1));
        handles.uitable1.Data(5,2)=cellstr(num2str(0));
        handles.uitable2.Data=(cellfun(@cellstr,comp.name))';
        handles.uitable2.Data=[(cellfun(@cellstr,comp.name))',cellstr(num2str(results.xji(1,:)')),cellstr(num2str(results.xji(1,:)'*1)),cellstr(num2str(results.xji(1,:)'*0))];

    elseif strcmp(varargin{1},'bottom')
        handles.stream.Name='bottom';
        handles.uitable1.Data=data;
        handles.uitable1.Data(1,2)=cellstr(num2str(results.Tjk(end,end)-273));
        handles.uitable1.Data(2,2)=cellstr(num2str(results.Pj(end)/1e3));
        handles.uitable1.Data(3,2)=cellstr(num2str(results.Lj(end)));
        handles.uitable1.Data(4,2)=cellstr(num2str(1));
        handles.uitable1.Data(5,2)=cellstr(num2str(0));
        handles.uitable2.Data=(cellfun(@cellstr,comp.name))';
        handles.uitable2.Data=[(cellfun(@cellstr,comp.name))',cellstr(num2str(results.xji(end,:)')),cellstr(num2str(results.xji(end,:)'*1)),cellstr(num2str(results.xji(end,:)'*0))];
        
    end
end

% Update handles structure
guidata(hObject, handles);
% --- Outputs from this function are returned to the command line.
function varargout = stream_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles;% hossein removed .output



% --- Executes when user attempts to close stream.
function stream_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to stream (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if strcmp(handles.stream.Name,'feed')
    global feed
    feed.data=handles.uitable1.Data;
    feed.comp=handles.uitable2.Data;
elseif strcmp(handles.stream.Name,'distil_vapor')
    2
elseif strcmp(handles.stream.Name,'distil_liquid')
    3
elseif strcmp(handles.stream.Name,'bottom')
    4
end
delete(hObject);


% --- Executes on key press with focus on uitable1 and none of its controls.
function uitable1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
% figure
% plot(1,1)


% --- Executes when entered data in editable cell(s) in uitable1.
function uitable1_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
% global feed
% if strcmp(handles.stream.Name,'feed')
%     feed.data=handles.uitable1.Data;
%     feed.comp=handles.uitable2.Data;
% end


% --- Executes on button press in flash.
function flash_Callback(~, eventdata, handles)
% hObject    handle to flash (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global feed comp
feed.data=handles.uitable1.Data;
feed.comp=handles.uitable2.Data;
T=str2double(cell2mat(feed.data(1,2)))+273; %K
P=str2double(cell2mat(feed.data(2,2)))*1000; %pa
z=str2double(cellfun(@cellstr,feed.comp(:,2)))';
z=abs(z)./sum(abs(z));
[B,x,y,zx,zy]=VL(T,P,z,comp);
feed.comp=[feed.comp(:,1),cellstr(num2str(z')),cellstr(num2str(x')),cellstr(num2str(y'))];
handles.uitable2.Data=feed.comp;
feed.data(4,2)=num2cell(B);
feed.data(5,2)=num2cell(1-B);
handles.uitable1.Data=feed.data;
