function varargout = rotating_beam(varargin)
% ROTATING_BEAM MATLAB code for rotating_beam.fig
%      ROTATING_BEAM, by itself, creates a new ROTATING_BEAM or raises the existing
%      singleton*.
%
%      H = ROTATING_BEAM returns the handle to a new ROTATING_BEAM or the handle to
%      the existing singleton*.
%
%      ROTATING_BEAM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ROTATING_BEAM.M with the given input arguments.
%
%      ROTATING_BEAM('Property','Value',...) creates a new ROTATING_BEAM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before rotating_beam_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to rotating_beam_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help rotating_beam

% Last Modified by GUIDE v2.5 04-Aug-2015 09:38:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @rotating_beam_OpeningFcn, ...
                   'gui_OutputFcn',  @rotating_beam_OutputFcn, ...
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


% --- Executes just before rotating_beam is made visible.
function rotating_beam_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to rotating_beam (see VARARGIN)

% Choose default command line output for rotating_beam
handles.output = hObject;

set(handles.crosssectionlistbox,'Value',1);
set(handles.materiallistbox,'Value',1);

clear_buttons(hObject, eventdata, handles);


clc;

fstr='rotating_beam.png';

bg = imread(fstr);
info = imfinfo(fstr);
w = info.Width;  %An integer indicating the width of the image in pixels
h = info.Height; %An integer indicating the height of the image in pixels
 
axes(handles.axes1);
image(bg);
pos1 = getpixelposition(handles.axes1,true);
 
set(handles.axes1, ...
    'Visible', 'off', ...
    'Units', 'pixels', ...
    'Position', [pos1(1) pos1(2) w h]);
axis off; 



handles.leftBC=1;
handles.rightBC=1;
handles.unit=1;
handles.material=1;
handles.cross_section=1;
 
handles.length=0;
handles.elastic_modulus=0;
handles.mass_density=0;
 
handles.thick=0;
handles.width=0;
 
handles.diameter=0;
handles.wall_thick=0;
 
handles.MOI=0;
handles.area=0;
 
handles.num=100;
handles.n=4;
 
handles.mode_number=1;
 
setappdata(0,'fig_num',1);
 
 
set(handles.edit_neutral_axis,'Visible','off');
set(handles.text_neutral_axis,'Visible','off');
 
set(handles.lengtheditbox,'String','');


set(handles.pushbutton_damping,'Visible','off');
set(handles.pushbutton_force_response,'Visible','off');

% Update handles structure
guidata(hObject, handles);
 
material_change(hObject, eventdata, handles);
geometry_change(hObject, eventdata, handles);



% UIWAIT makes rotating_beam wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = rotating_beam_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in calculatebutton.
function calculatebutton_Callback(hObject, eventdata, handles)
% hObject    handle to calculatebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tpi=2*pi;

omega=str2num(get(handles.edit_rate,'String'));
nrate=get(handles.listbox_rate,'Value');

if(nrate==1) % rpm to rad/sec
    omega=tpi*omega/60;
end
if(nrate==2) % Hz to rad/sec
    omega=tpi*omega;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setappdata(0,'fig_num',1);

handles.cross_section=get(handles.crosssectionlistbox,'Value');
handles.unit=get(handles.unitslistbox,'Value');

E=str2num(get(handles.elasticmoduluseditbox,'String'));
rho=str2num(get(handles.massdensityeditbox,'String'));

nsm=str2num(get(handles.edit_nsm,'String'));

L=str2num(get(handles.lengtheditbox,'String'));

if(L>1.0e-20 && L<1.0e+20)
else
    warndlg('Enter Valid length');
    return;
end

thick=handles.thick;
width=handles.width;
diameter=handles.diameter;
wall_thick=handles.wall_thick;
    
if(handles.unit==1) % English
    rho=rho/386;
    nsm=nsm/386;
else                % metric
    [E]=GPa_to_Pa(E);
    thick=thick/1000;
    width=width/1000;
    diameter=diameter/1000;
    wall_thick=wall_thick/1000;
end
    
setappdata(0,'unit',handles.unit);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(handles.cross_section==4) % other
        set(handles.edit_neutral_axis,'Enable','on');
else
        set(handles.edit_neutral_axis,'Enable','off');        
end
    

    
if(handles.cross_section==1) % rectangular
        [area,MOI,cna]=beam_rectangular_geometry(width,thick); 
        setappdata(0,'chord',width);
end
if(handles.cross_section==2) % pipe
        [area,MOI,cna]=pipe_geometry_wall(diameter,wall_thick);    
        setappdata(0,'chord',diameter);        
end
if(handles.cross_section==3) % solid cylinder
        [area,MOI,cna]=cylinder_geometry(diameter); 
        setappdata(0,'chord',diameter);        
end
if(handles.cross_section==4) % other
        area=handles.area;
        MOI=handles.MOI;
        cna=str2num(get(handles.edit_neutral_axis,'String'));
        chord=sqrt(4*area/pi);
        setappdata(0,'chord',chord);
end


ss=sprintf('%8.4g',cna);
set(handles.edit_neutral_axis,'String',ss); 
setappdata(0,'cna',cna);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out1=sprintf('area=%8.4g MOI=%8.4g',area,MOI);
disp(out1);

rho=rho*area+(nsm/L);   % mass per unit length
    
tmass=rho*L;
    
ne=str2num(get(handles.edit_num_elem,'String'));
setappdata(0,'ne',ne);

deltax=L/ne;

dx=zeros(ne,1);

for i=1:ne
    dx(i)=deltax;
end

xx=zeros((ne+1),1);

for i=2:(ne+1)
   xx(i)=xx(i-1)+deltax; 
end

setappdata(0,'xx',xx);

dof=2*ne+2;
orig_dof=dof;

setappdata(0,'dof',dof);

%
I=MOI;

[klocal] = vibrationdata_local_stiffness_rotating(E,I,ne,dx,omega,L,rho);
[mlocal] = local_mass(rho,ne,dx);

[stiff,mass,dof,dof_status]=rotating_beam_assembly(ne,klocal,mlocal,dof);
          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          

if(handles.unit==1) % English    
    out1=sprintf('\n Total mass = %8.4g lbm ',tmass*386);
else
    out1=sprintf('\n Total mass = %8.4g kg ',tmass);        
end
    
disp(out1);    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[fn,omega,ModeShapes,MST]=Generalized_Eigen(stiff,mass,2);
%

assignin('base','rstiff',stiff);
assignin('base','rmass',mass);

%
disp(' ');
disp('Mode   Freq (Hz) ');

for i=1:dof
    out1 = sprintf('%d  %10.4g ',i,fn(i));
    disp(out1)
    
    if(i==12)
        break;
    end
end

setappdata(0,'fn',fn(1:12));
setappdata(0,'ModeShapes',ModeShapes(:,1:12));
setappdata(0,'dx',dx);
setappdata(0,'L',L);
setappdata(0,'dof_status',dof_status);


set(handles.pushbutton_view_mode_shapes,'Visible','on');
set(handles.uipanel_mode_number,'Visible','on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=length(fn);

ModeShapes_full=zeros(orig_dof,N);

for j=1:N
    
    pf_mode=ModeShapes(:,j); 
    
    for i=1:orig_dof
        if(dof_status(i)==0)
%       
            pf_mode=[pf_mode(1:(i-1));0;pf_mode(i:length(pf_mode))];
%
        end
    end    
 
    max(pf_mode);
    ModeShapes_full(:,j)=pf_mode;
    
end 

setappdata(0,'ModeShapes_full',ModeShapes_full);

set(handles.pushbutton_damping,'Visible','on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


msgbox('Natural Frequencies written to Command Window');

 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on selection change in leftBClistbox.
function leftBClistbox_Callback(hObject, eventdata, handles)
% hObject    handle to leftBClistbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns leftBClistbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from leftBClistbox
clear_buttons(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function leftBClistbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to leftBClistbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in crosssectionlistbox.
function crosssectionlistbox_Callback(hObject, eventdata, handles)
% hObject    handle to crosssectionlistbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns crosssectionlistbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from crosssectionlistbox
clear_buttons(hObject, eventdata, handles);
handles.cross_section=get(hObject,'Value');
guidata(hObject, handles);
geometry_change(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function crosssectionlistbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to crosssectionlistbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function crosssectioneditbox1_Callback(hObject, eventdata, handles)
% hObject    handle to crosssectioneditbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of crosssectioneditbox1 as text
%        str2double(get(hObject,'String')) returns contents of crosssectioneditbox1 as a double

clear_buttons(hObject, eventdata, handles);

handles.thick=0;
handles.diameter=0;
handles.area=0;

string=get(hObject,'String');

if(handles.cross_section==1) % rectangular
    handles.thick=str2num(string);
end
if(handles.cross_section==2) % pipe
    handles.diameter=str2num(string);   
end
if(handles.cross_section==3) % solid cylinder
    handles.diameter=str2num(string);    
end
if(handles.cross_section==4) % other
    handles.area=str2num(string);     
end


guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function crosssectioneditbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to crosssectioneditbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_neutral_axis_Callback(hObject, eventdata, handles)
% hObject    handle to edit_neutral_axis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_neutral_axis as text
%        str2double(get(hObject,'String')) returns contents of edit_neutral_axis as a double


% --- Executes during object creation, after setting all properties.
function edit_neutral_axis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_neutral_axis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function crosssectioneditbox2_Callback(hObject, eventdata, handles)
% hObject    handle to crosssectioneditbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of crosssectioneditbox2 as text
%        str2double(get(hObject,'String')) returns contents of crosssectioneditbox2 as a double

clear_buttons(hObject, eventdata, handles);

string=get(hObject,'String');

if(handles.cross_section==1) % rectangular
    handles.width=str2num(string);
end
if(handles.cross_section==2) % pipe
    handles.wall_thick=str2num(string);   
end
if(handles.cross_section==3) % solid cylinder
    % nothing to do    
end
if(handles.cross_section==4) % other
    handles.MOI=str2num(string);     
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function crosssectioneditbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to crosssectioneditbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in materiallistbox.
function materiallistbox_Callback(hObject, eventdata, handles)
% hObject    handle to materiallistbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns materiallistbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from materiallistbox
clear_buttons(hObject, eventdata, handles);
handles.material=get(hObject,'Value');
material_change(hObject, eventdata, handles);
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function materiallistbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to materiallistbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function elasticmoduluseditbox_Callback(hObject, eventdata, handles)
% hObject    handle to elasticmoduluseditbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of elasticmoduluseditbox as text
%        str2double(get(hObject,'String')) returns contents of elasticmoduluseditbox as a double

clear_buttons(hObject, eventdata, handles);

string=get(hObject,'String');
handles.elastic_modulus=str2num(string);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function elasticmoduluseditbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to elasticmoduluseditbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function massdensityeditbox_Callback(hObject, eventdata, handles)
% hObject    handle to massdensityeditbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of massdensityeditbox as text
%        str2double(get(hObject,'String')) returns contents of massdensityeditbox as a double


clear_buttons(hObject, eventdata, handles);

string=get(hObject,'String');

handles.mass_density=str2num(string);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function massdensityeditbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to massdensityeditbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function lengtheditbox_Callback(hObject, eventdata, handles)
% hObject    handle to lengtheditbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lengtheditbox as text
%        str2double(get(hObject,'String')) returns contents of lengtheditbox as a double
clear_buttons(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function lengtheditbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lengtheditbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in unitslistbox.
function unitslistbox_Callback(hObject, eventdata, handles)
% hObject    handle to unitslistbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns unitslistbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from unitslistbox

clear_buttons(hObject, eventdata, handles);

handles.unit=get(handles.unitslistbox,'Value');

guidata(hObject, handles);

material_change(hObject, eventdata, handles);
geometry_change(hObject, eventdata, handles);


if(handles.unit==1)
    set(handles.lengthlabel,'String','inch');
    set(handles.elasticmoduluslabel,'String','Elastic Modulus (lbf/in^2)');
    set(handles.massdensitylabel,'String','Mass Density (lbm/in^3)'); 
    set(handles.text_nsm,'String','lbm');        
else
    set(handles.lengthlabel,'String','meter');
    set(handles.elasticmoduluslabel,'String','Elastic Modulus (GPa)'); 
    set(handles.massdensitylabel,'String','Mass Density (kg/m^3)'); 
    set(handles.text_nsm,'String','kg');
end

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function unitslistbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to unitslistbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_k2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_k2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_k2 as text
%        str2double(get(hObject,'String')) returns contents of edit_k2 as a double
clear_buttons(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function edit_k2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_k2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_k1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_k1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_k1 as text
%        str2double(get(hObject,'String')) returns contents of edit_k1 as a double
clear_buttons(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function edit_k1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_k1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in rightBClistbox.
function rightBClistbox_Callback(hObject, eventdata, handles)
% hObject    handle to rightBClistbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns rightBClistbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from rightBClistbox
clear_buttons(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function rightBClistbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rightBClistbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function geometry_change(hObject, eventdata, handles)

set(handles.edit_neutral_axis,'Visible','off');
set(handles.edit_neutral_axis,'String','');
set(handles.edit_neutral_axis,'Enable','off');

set(handles.text_neutral_axis,'Visible','off');

        

set(handles.crosssectioneditbox2,'Enable','on');

if(handles.cross_section==4) % other
    set(handles.edit_neutral_axis,'Visible','on');
    set(handles.edit_neutral_axis,'Enable','on');    
    set(handles.text_neutral_axis,'Visible','on');
end

if(handles.unit==1)
    if(handles.cross_section==1) % rectangular
        set(handles.crosssectionlabel1,'String','Thickness (in)');
        set(handles.crosssectionlabel2','String','Width (in)');
    end
    if(handles.cross_section==2) % pipe 
        set(handles.crosssectionlabel1,'String','Outer Diameter (in)');
        set(handles.crosssectionlabel2,'String','Wall Thickness (in)');    
    end
    if(handles.cross_section==3) % solid cylinder
        set(handles.crosssectionlabel1,'String','Diameter (in)');
        set(handles.crosssectionlabel2,'String',''); 
        set(handles.crosssectioneditbox2,'Enable','off');
    end
    if(handles.cross_section==4) % other
        set(handles.crosssectionlabel1,'String','Area (in^2)');
        set(handles.crosssectionlabel2,'String','Area Moment of Inertia(in^4)');
        set(handles.text_neutral_axis,'String','Distance from Neutral Axis (inches)');
    end
else
    if(handles.cross_section==1) % rectangular
        set(handles.crosssectionlabel1,'String','Thickness (mm)');
        set(handles.crosssectionlabel2','String','Width (mm)');
    end
    if(handles.cross_section==2) % pipe 
        set(handles.crosssectionlabel1,'String','Outer Diameter (mm)');
        set(handles.crosssectionlabel2,'String','Wall Thickness (mm)');    
    end
    if(handles.cross_section==3) % solid cylinder
        set(handles.crosssectionlabel1,'String','Diameter (mm)');
        set(handles.crosssectionlabel2,'String',''); 
        set(handles.crosssectioneditbox2,'Enable','off');
    end
    if(handles.cross_section==4) % other
        set(handles.crosssectionlabel1,'String','Area (mm^2)');
        set(handles.crosssectionlabel2,'String','Area Moment of Inertia(mm^4)');
        set(handles.text_neutral_axis,'String','Distance from Neutral Axis (mm)');        
    end
end

set(handles.crosssectioneditbox1,'String',' ');
set(handles.crosssectioneditbox2,'String',' ');

handles.thick=0;
handles.width=0;

handles.diameter=0;
handles.wall_thick=0;

handles.MOI=0;
handles.area=0;


handles.material=get(handles.materiallistbox,'Value');
material_change(hObject, eventdata, handles);

guidata(hObject, handles);


function material_change(hObject, eventdata, handles)

if(handles.unit==1)  % English
    if(handles.material==1) % aluminum
        handles.elastic_modulus=1e+007;
        handles.mass_density=0.1;  
    end  
    if(handles.material==2)  % steel
        handles.elastic_modulus=3e+007;
        handles.mass_density= 0.28;         
    end
    if(handles.material==3)  % copper
        handles.elastic_modulus=1.6e+007;
        handles.mass_density=  0.322;
    end
    if(handles.material==4)  % G10
        handles.elastic_modulus=2.7e+006;
        handles.mass_density=  0.065;
    end
    if(handles.material==5)  % PVC
        handles.elastic_modulus=3.5e+005;
        handles.mass_density=  0.052;
    end
else                 % metric
    if(handles.material==1)  % aluminum
        handles.elastic_modulus=70;
        handles.mass_density=  2700;
    end
    if(handles.material==2)  % steel
        handles.elastic_modulus=205;
        handles.mass_density=  7700;        
    end
    if(handles.material==3)   % copper
        handles.elastic_modulus=110;
        handles.mass_density=  8900;
    end
    if(handles.material==4)  % G10
        handles.elastic_modulus=18.6;
        handles.mass_density=  1800;
    end
    if(handles.material==5)  % PVC
        handles.elastic_modulus=24.1;
        handles.mass_density=  1440;
    end
end

if(handles.material<6)
    ss1=sprintf('%8.4g',handles.elastic_modulus);
    ss2=sprintf('%8.4g',handles.mass_density);

    set(handles.elasticmoduluseditbox,'String',ss1);
    set(handles.massdensityeditbox,'String',ss2);    
end

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function listbox6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox5.
function listbox5_Callback(hObject, eventdata, handles)
% hObject    handle to listbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox5


% --- Executes during object creation, after setting all properties.
function listbox5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox4.
function listbox4_Callback(hObject, eventdata, handles)
% hObject    handle to listbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox4


% --- Executes during object creation, after setting all properties.
function listbox4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox3.
function listbox3_Callback(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox3


% --- Executes during object creation, after setting all properties.
function listbox3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_ne.
function listbox_ne_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_ne (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_ne contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_ne
clear_buttons(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function listbox_ne_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_ne (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function clear_buttons(hObject, eventdata, handles)
%

set(handles.pushbutton_view_mode_shapes,'Visible','off');
set(handles.uipanel_mode_number,'Visible','off');

set(handles.pushbutton_damping,'Visible','off');
set(handles.pushbutton_force_response,'Visible','off');




% --- Executes on button press in pushbutton_return.
function pushbutton_return_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_return (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(rotating_beam);


% --- Executes on button press in pushbutton_view_mode_shapes.
function pushbutton_view_mode_shapes_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_view_mode_shapes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

iu=get(handles.unitslistbox,'Value');

fn=getappdata(0,'fn');
ModeShapes_full=getappdata(0,'ModeShapes_full');

dx=getappdata(0,'dx');
L=getappdata(0,'L');
dof_status=getappdata(0,'dof_status');

dxa=mean(dx);

num=get(handles.listbox_mode_number,'Value');

fig_num=getappdata(0,'fig_num');


orig_dof=length(dof_status);
dof_half=round(orig_dof/2);

%
pf_mode=ModeShapes_full(:,num);  
  
%
pmode=zeros(dof_half);
%        
for i=1:dof_half
    pmode(i)=pf_mode(2*i-1);
end

M=pmode/(max(abs(pmode)));
N=length(M);

x=zeros(N,1);
y=zeros(N,1);
    
for i=2:N
    x(i)=x(i-1)+dxa;
    y(i)=M(i-1);
end    

figure(fig_num);
plot(x,y);
out1=sprintf('Mode %d  fn=%7.3g Hz',num,fn(num));
title(out1);

if(iu==1)
    xlabel('x (in)');
else
    xlabel('x (m)');    
end

fig_num=fig_num+1;
setappdata(0,'fig_num',fig_num);





% --- Executes on selection change in listbox_mode_number.
function listbox_mode_number_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_mode_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_mode_number contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_mode_number


% --- Executes during object creation, after setting all properties.
function listbox_mode_number_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_mode_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_damping.
function pushbutton_damping_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_damping (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


handles.s=generic_damping;    
set(handles.s,'Visible','on'); 


set(handles.pushbutton_force_response,'Visible','on');


% --- Executes on button press in pushbutton_transfer.
function pushbutton_transfer_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_transfer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.s=beam_spring_force_frf;    
set(handles.s,'Visible','on'); 



function edit_nsm_Callback(hObject, eventdata, handles)
% hObject    handle to edit_nsm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_nsm as text
%        str2double(get(hObject,'String')) returns contents of edit_nsm as a double


% --- Executes during object creation, after setting all properties.
function edit_nsm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_nsm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_sine.
function pushbutton_sine_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_sine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.s=beam_spring_sine;    
set(handles.s,'Visible','on'); 



function edit_pml1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_pml1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_pml1 as text
%        str2double(get(hObject,'String')) returns contents of edit_pml1 as a double


% --- Executes during object creation, after setting all properties.
function edit_pml1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_pml1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_pm1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_pm1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_pm1 as text
%        str2double(get(hObject,'String')) returns contents of edit_pm1 as a double


% --- Executes during object creation, after setting all properties.
function edit_pm1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_pm1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_pml2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_pml2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_pml2 as text
%        str2double(get(hObject,'String')) returns contents of edit_pml2 as a double


% --- Executes during object creation, after setting all properties.
function edit_pml2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_pml2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_pml3_Callback(hObject, eventdata, handles)
% hObject    handle to edit_pml3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_pml3 as text
%        str2double(get(hObject,'String')) returns contents of edit_pml3 as a double


% --- Executes during object creation, after setting all properties.
function edit_pml3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_pml3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_pm2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_pm2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_pm2 as text
%        str2double(get(hObject,'String')) returns contents of edit_pm2 as a double


% --- Executes during object creation, after setting all properties.
function edit_pm2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_pm2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_pm3_Callback(hObject, eventdata, handles)
% hObject    handle to edit_pm3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_pm3 as text
%        str2double(get(hObject,'String')) returns contents of edit_pm3 as a double


% --- Executes during object creation, after setting all properties.
function edit_pm3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_pm3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_npm.
function listbox_npm_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_npm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_npm contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_npm


n=get(handles.listbox_npm,'Value');
n=n-1;

if(n==0)
    set(handles.text_pmloc,'Visible','off');
    set(handles.text_pm,'Visible','off');
else
    set(handles.text_pmloc,'Visible','on');
    set(handles.text_pm,'Visible','on');    
end    

set(handles.edit_pm1,'Visible','off');
set(handles.edit_pm2,'Visible','off');
set(handles.edit_pm3,'Visible','off');

set(handles.edit_pml1,'Visible','off');
set(handles.edit_pml2,'Visible','off');
set(handles.edit_pml3,'Visible','off');

if(n>=1)
    set(handles.edit_pm1,'Visible','on');
    set(handles.edit_pml1,'Visible','on');    
end
if(n>=2)
    set(handles.edit_pm2,'Visible','on');
    set(handles.edit_pml2,'Visible','on');     
end
if(n==3)
    set(handles.edit_pm3,'Visible','on');
    set(handles.edit_pml3,'Visible','on');     
end
    






% --- Executes during object creation, after setting all properties.
function listbox_npm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_npm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function uipanel11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function edit_num_elem_Callback(hObject, eventdata, handles)
% hObject    handle to edit_num_elem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_num_elem as text
%        str2double(get(hObject,'String')) returns contents of edit_num_elem as a double


% --- Executes during object creation, after setting all properties.
function edit_num_elem_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_num_elem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_rate_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rate as text
%        str2double(get(hObject,'String')) returns contents of edit_rate as a double


% --- Executes during object creation, after setting all properties.
function edit_rate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_rate.
function listbox_rate_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_rate contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_rate


% --- Executes during object creation, after setting all properties.
function listbox_rate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_damping.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_damping (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.s=generic_damping;    
set(handles.s,'Visible','on'); 


% --- Executes on button press in pushbutton_force_response.
function pushbutton_force_response_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_force_response (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.s=rotating_beam_applied_force;    
set(handles.s,'Visible','on'); 
