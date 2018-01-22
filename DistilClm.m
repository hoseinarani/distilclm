function varargout = DistilClm(varargin)
% DISTILCLM MATLAB code for DistilClm.fig
%      DISTILCLM, by itself, creates a new DISTILCLM or raises the existing
%      singleton*.
%
%      H = DISTILCLM returns the handle to a new DISTILCLM or the handle to
%      the existing singleton*.
%
%      DISTILCLM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DISTILCLM.M with the given input arguments.
%
%      DISTILCLM('Property','Value',...) creates a new DISTILCLM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DistilClm_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DistilClm_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DistilClm

% Last Modified by GUIDE v2.5 16-Aug-2017 21:21:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @DistilClm_OpeningFcn, ...
    'gui_OutputFcn',  @DistilClm_OutputFcn, ...
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


% --- Executes just before DistilClm is made visible.
function DistilClm_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DistilClm (see VARARGIN)

% Choose default command line output for DistilClm
handles.output = hObject;
global feed D_v D_l B_l clm sel_comp all_comp prop data clm
clm.n_tray=10;
clm.f_tray=3;
clm.rr=1.4;
clm.dp=.7;
clm.distil_rate=10;

sel_comp={};
% which Activate
%system('taskkill /F /IM EXCEL.EXE');
% num = readtable('data.xlsx');
data=loaddata;
all_comp=data(:,1);
feed.data={'T (C)','';
    'P (kpa)','';
    'Molar flow rate','';
    'Liquid molar fraction','';
    'Vapor molar fraction',''};;
% feed.comp{1,4}=[];
% D_v.comp{1,4}=[];
% D_l.comp{1,4}=[];
% B_l.comp{1,4}=[];

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DistilClm wait for user response (see UIRESUME)
% uiwait(handles.figure1);
axes(handles.column)
matlabImage = imread('distillation-column.png');
image(matlabImage)
axis off
axis image

% --- Outputs from this function are returned to the command line.
function varargout = DistilClm_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in feed.
function feed_Callback(hObject, eventdata, handles)
% hObject    handle to feed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% uf = uifigure;
% t = uitable(uf);
% d = {'Male',52,'kcal/mol';'Male',40,'K';'Female',25,'Kg'};
% t.Data = d;
% t.Position = [20 20 258 78];
% t.ColumnName = {'Gender','Age','Authorized'};
% t.ColumnEditable = true;
stream('feed');



% --- Executes on button press in distil_v.
function distil_v_Callback(hObject, eventdata, handles)
% hObject    handle to distil_v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
stream('distil_v');

% --- Executes on mouse press over axes background.
function column_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to column (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

plot(1,1)


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cp=get (hObject, 'CurrentPoint');

if campare_position(handles.column.Position,cp)
    column();
else

end


% --- Executes on button press in distil_l.
function distil_l_Callback(hObject, eventdata, handles)
% hObject    handle to distil_l (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
stream('distil_l');


% --- Executes on button press in bottom_l.
function bottom_l_Callback(hObject, eventdata, handles)
% hObject    handle to bottom_l (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
stream('bottom');
%setAlwaysOnTop(f.stream,'true')
%uistack(f.stream,'top');
% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over feed.
function feed_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to feed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on mouse motion over figure - except title and menu.
function figure1_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cp=get (hObject, 'CurrentPoint');
% if handles.feed.Position(1)<cp(1) && handles.feed.Position(1)+handles.feed.Position(3)>cp(1) ...
%         && handles.feed.Position(2)<cp(2) && handles.feed.Position(2)+handles.feed.Position(4)>cp(2)
%     handles.feed.Position
if campare_position(handles.feed.Position,cp)
    handles.feed.BackgroundColor=[1 1 1];
else
    handles.feed.BackgroundColor=[0.94 0.94 0.94];
end
if campare_position(handles.distil_v.Position,cp)
    handles.distil_v.BackgroundColor=[1 1 1];
else
    handles.distil_v.BackgroundColor=[0.94 0.94 0.94];
end
if campare_position(handles.distil_l.Position,cp)
    handles.distil_l.BackgroundColor=[1 1 1];
else
    handles.distil_l.BackgroundColor=[0.94 0.94 0.94];
end
if campare_position(handles.bottom_l.Position,cp)
    handles.bottom_l.BackgroundColor=[1 1 1];
else
    handles.bottom_l.BackgroundColor=[0.94 0.94 0.94];
end

if campare_position(handles.energy_d.Position,cp)
    handles.energy_d.BackgroundColor=[1 1 1];
else
    handles.energy_d.BackgroundColor=[0.94 0.94 0.94];
end
if campare_position(handles.energy_b.Position,cp)
    handles.energy_b.BackgroundColor=[1 1 1];
else
    handles.energy_b.BackgroundColor=[0.94 0.94 0.94];
end

column_Position=[handles.column.Position(1) handles.column.Position(2) , ...
    handles.column.Position(3)/2, handles.column.Position(4)];
if campare_position(column_Position,cp)
    axes(handles.column)
    matlabImage = imread('distillation-column1.png');
    image(matlabImage)
    axis off
    axis image
else
    axes(handles.column)
    matlabImage = imread('distillation-column.png');
    image(matlabImage)
    axis off
    axis image
end

function answ=campare_position(position,cp)
if  position(1)<cp(1) && position(1)+position(3)>cp(1) ...
        && position(2)<cp(2) && position(2)+position(4)>cp(2)
    answ=1;
else
    answ=0;
end


% --- Executes on mouse press over figure background.
function figure1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function select_comp_Callback(hObject, ~, handles)
% hObject    handle to select_comp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
components();

% --------------------------------------------------------------------
function file_Callback(hObject, eventdata, handles)
% hObject    handle to file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);
clear
clear all


% --- Executes on button press in runclm.
function runclm_Callback(hObject, eventdata, handles)
% hObject    handle to runclm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global comp clm feed results
[Lj,Vj,xji,yji,Tjk,Qj,Uj,Wj,Fj,Pj,N,nc]=columncalc(clm,comp,feed);
results.Lj=Lj;
results.Vj=Vj;
results.xji=xji;
results.yji=yji;
results.Tjk=Tjk;
results.Qj=Qj;
results.Uj=Uj;
results.Wj=Wj;
results.Fj=Fj;
results.Pj=Pj;
results.N=N;
results.nc=nc;

% --- Executes on button press in plt_T.
function plt_T_Callback(hObject, eventdata, handles)
% hObject    handle to plt_T (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global results
fprintf('\n***************************************************\n')
fprintf('\n temperature K :\n ');
fprintf('===============');
for j=1:results.N
    fprintf('\n stage %d \t %f',j,results.Tjk(j,end));
end
fprintf('\n***************************************************\n')
figure;
plot(1:results.N,results.Tjk(:,end)-273.15);
grid 
title(['Stage 1 is condenser and stage ',num2str(results.N),' is reboiler'])
xlabel('Stage')
ylabel('Temperature (C)')
% --- Executes on button press in plt_flow.
function plt_flow_Callback(hObject, eventdata, handles)
% hObject    handle to plt_flow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global results
%fprintf('feed composition :\n C3=%f \n n-C4=%f \n n-C5=%f \n',zF(1),zF(2),zF(3))
fprintf('liquid distilate flow rate kmol/h = %f \n',results.Uj(1))
%fprintf('vapore distilate flow rate kmol/h = %f \n',V1)

figure
subplot(2,1,1);
plot(1:results.N,results.Lj(:))
grid
xlabel('Stage')
ylabel('Liquid flow rate')
title(['Stage 1 is condenser and stage ',num2str(results.N),' is reboiler'])

subplot(2,1,2);
plot(1:results.N,results.Vj(:))
grid
xlabel('Stage')
ylabel('Vapore flow rate')

%fprintf('\nvapore composition: \n \t\t\t C3 \t\t n-C4 \t\t n-C5 ');
fprintf('\nvapore composition: \n');

% --- Executes on button press in plt_vap_comp.
function plt_vap_comp_Callback(hObject, eventdata, handles)
% hObject    handle to plt_vap_comp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global results comp
fprintf('***************************************************\n')
disp('Output Data:')
fprintf('============\n');
fprintf('\n Vapor composition: \n ');
fprintf('***************************************************\n')
fprintf('\t \t \t');

for i=1:results.nc
    fprintf('%s \t', string(comp.name{i}));
end
fprintf('\n');
for j=1:results.N
    fprintf('\n stage %d \t',j);
    for i=1:results.nc
        fprintf('%f \t',results.yji(j,i));
    end
end
fprintf('\n***************************************************\n')
figure
h=plot(1:results.N,results.yji(:,:));
grid
xlabel('stage')
ylabel('Vapore composition')
for i=1:size(comp.name,2)
    name{i}=cell2mat(comp.name{i});
end
h = legend(name);
% --- Executes on button press in plt_liq_comp.
function plt_liq_comp_Callback(hObject, eventdata, handles)
% hObject    handle to plt_liq_comp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global results comp
fprintf('***************************************************\n')
disp('Output Data:')
fprintf('============\n');
fprintf('\n Liquid composition: \n ');
fprintf('***************************************************\n')
fprintf('\t \t \t');

for i=1:results.nc
    fprintf('%s \t', string(comp.name{i}));
end
fprintf('\n');
for j=1:results.N
    fprintf('\n stage %d \t',j);
    for i=1:results.nc
        fprintf('%f \t',results.xji(j,i));
    end
end
fprintf('\n***************************************************\n')
figure
h=plot(1:results.N,results.xji(:,:));
grid
xlabel('stage')
ylabel('Liquid composition y')
for i=1:size(comp.name,2)
    name{i}=cell2mat(comp.name{i});
end
h = legend(name);


% --- Executes on button press in energy_d.
function energy_d_Callback(hObject, eventdata, handles)
% hObject    handle to energy_d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global results

h = msgbox(['Condenser Energy = ', sprintf('%1.3e',results.Qj(1)),' kJ/kgmol'],'modal');


% --- Executes on button press in energy_b.
function energy_b_Callback(hObject, eventdata, handles)
% hObject    handle to energy_b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global results
h = msgbox(['Reboiler Energy = ', sprintf('%1.3e',results.Qj(end)),' kJ/kgmol'],'modal');
