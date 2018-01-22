function varargout = components(varargin)
% COMPONENTS MATLAB code for components.fig
%      COMPONENTS, by itself, creates a new COMPONENTS or raises the existing
%      singleton*.
%
%      H = COMPONENTS returns the handle to a new COMPONENTS or the handle to
%      the existing singleton*.
%
%      COMPONENTS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COMPONENTS.M with the given input arguments.
%
%      COMPONENTS('Property','Value',...) creates a new COMPONENTS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before components_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to components_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help components

% Last Modified by GUIDE v2.5 13-Aug-2017 12:10:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @components_OpeningFcn, ...
    'gui_OutputFcn',  @components_OutputFcn, ...
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


% --- Executes just before components is made visible.
function components_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to components (see VARARGIN)

% Choose default command line output for components
handles.output = hObject;
global all_comp sel_comp
handles.all_comp.Value=size(all_comp,2);
if isempty(sel_comp)
    handles.sel_comp.Value=0;
else
    handles.sel_comp.Value=size(sel_comp,2);
end
handles.all_comp.String=all_comp;
handles.sel_comp.String=sel_comp;




% Update handles structure
guidata(hObject, handles);

% UIWAIT makes components wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = components_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in all_comp.
function all_comp_Callback(hObject, eventdata, handles)
% hObject    handle to all_comp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns all_comp contents as cell array
%        contents{get(hObject,'Value')} returns selected item from all_comp


% --- Executes during object creation, after setting all properties.
function all_comp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to all_comp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in sel_comp.
function sel_comp_Callback(hObject, eventdata, handles)
% hObject    handle to sel_comp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns sel_comp contents as cell array
%        contents{get(hObject,'Value')} returns selected item from sel_comp


% --- Executes during object creation, after setting all properties.
function sel_comp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sel_comp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over all_comp.
function all_comp_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to all_comp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global all_comp sel_comp
indx=find(strcmp(all_comp,hObject.String{hObject.Value}));
sel_comp=[sel_comp;all_comp{indx}];
all_comp(indx)=[];
handles.all_comp.Value=size(all_comp,2);
handles.sel_comp.Value=size(sel_comp,2);
all_comp(~cellfun('isempty',all_comp));
sel_comp(~cellfun('isempty',sel_comp));
all_comp=sort(all_comp);
sel_comp=sort(sel_comp);
handles.all_comp.String=all_comp;
handles.sel_comp.String=sel_comp;



% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over sel_comp.
function sel_comp_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to sel_comp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global all_comp sel_comp
indx=find(strcmp(sel_comp,hObject.String{hObject.Value}));
all_comp=[all_comp;sel_comp{indx}];
sel_comp(indx)=[];
handles.all_comp.Value=size(all_comp,2);
handles.sel_comp.Value=size(sel_comp,2);
all_comp(~cellfun('isempty',all_comp));
sel_comp(~cellfun('isempty',sel_comp));
all_comp=sort(all_comp);
sel_comp=sort(sel_comp);
handles.all_comp.String=all_comp;
handles.sel_comp.String=sel_comp;


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global feed D_v D_l B_l sel_comp data comp
com{size(sel_comp,1),2}=[];

feed.comp=[sel_comp,com];
D_v.comp=[sel_comp,com];
D_l.comp=[sel_comp,com];
B_l.comp=[sel_comp,com];
comp=[];
for i=1:size(sel_comp,1)
    idx=find(strcmp(data(:,1), sel_comp(i)));
    comp.name{i}=data(idx,1);
    comp.Tc(i)=cell2mat(data(idx,4))+273;%K
    comp.Pc(i)=cell2mat(data(idx,3))*1e3;%pa
    comp.w(i)=cell2mat(data(idx,2));
    comp.Mw(i)=cell2mat(data(idx,6));
    comp.hoffset(i)=cell2mat(data(idx,19));   
    comp.hv{i}=data(idx,7:12);
    comp.an{i}=data(idx,13:18);   
end

% Hint: delete(hObject) closes the figure
delete(hObject);
