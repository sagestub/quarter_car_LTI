function varargout = ri_mdof_modal_arbitrary_force(varargin)
% RI_MDOF_MODAL_ARBITRARY_FORCE MATLAB code for ri_mdof_modal_arbitrary_force.fig
%      RI_MDOF_MODAL_ARBITRARY_FORCE, by itself, creates a new RI_MDOF_MODAL_ARBITRARY_FORCE or raises the existing
%      singleton*.
%
%      H = RI_MDOF_MODAL_ARBITRARY_FORCE returns the handle to a new RI_MDOF_MODAL_ARBITRARY_FORCE or the handle to
%      the existing singleton*.
%
%      RI_MDOF_MODAL_ARBITRARY_FORCE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RI_MDOF_MODAL_ARBITRARY_FORCE.M with the given input arguments.
%
%      RI_MDOF_MODAL_ARBITRARY_FORCE('Property','Value',...) creates a new RI_MDOF_MODAL_ARBITRARY_FORCE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ri_mdof_modal_arbitrary_force_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ri_mdof_modal_arbitrary_force_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ri_mdof_modal_arbitrary_force

% Last Modified by GUIDE v2.5 28-Jul-2017 14:55:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ri_mdof_modal_arbitrary_force_OpeningFcn, ...
                   'gui_OutputFcn',  @ri_mdof_modal_arbitrary_force_OutputFcn, ...
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


% --- Executes just before ri_mdof_modal_arbitrary_force is made visible.
function ri_mdof_modal_arbitrary_force_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ri_mdof_modal_arbitrary_force (see VARARGIN)

% Choose default command line output for ri_mdof_modal_arbitrary_force
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ri_mdof_modal_arbitrary_force wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ri_mdof_modal_arbitrary_force_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_calculate.
function pushbutton_calculate_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_calculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

imethod=get(handles.listbox_method,'Value');

iu=get(handles.listbox_units,'value');
damp=str2num(get(handles.edit_damping,'String'));



FS=get(handles.edit_mass,'String');
M=evalin('base',FS); 

FS=get(handles.edit_stiffness,'String');
K=evalin('base',FS); 

FS=get(handles.edit_force,'String');
F=evalin('base',FS); 


num_modes=str2num(get(handles.edit_num_modes,'String'));

if(isempty(num_modes))
   warndlg('enter number of array '); 
   return; 
end 


if(isempty(M))
   warndlg('enter mass array'); 
   return; 
end  
if(isempty(K))
   warndlg('enter stiffness array'); 
   return; 
end 
if(isempty(F))
   warndlg('enter force array'); 
   return; 
end 
if(isempty(num_modes))
   warndlg('enter number of array '); 
   return; 
end 


szM=size(M);
szK=size(K);
szF=size(F);

num_force_col=szF(2)-1;

NT=szF(1);

sys_dof=szM(1);

if(num_force_col ~= sys_dof)
   warndlg('number of force dof not equal to number of system dof'); 
   return; 
end    

if(num_modes>sys_dof)
   num_modes=sys_dof;
   ns=sprintf('%d',num_modes);
   set(handles.edit_num_modes,'String',ns);
end

ndof=num_modes;

dur=F(NT,1)-F(1,1);

dt=dur/(NT-1);

out1=sprintf(' NT=%d  dur=%8.4g sec  dt=%8.4g sec',NT,dur,dt);
disp(out1);

%
disp(' ');
[fn,omegan,ModeShapes,MST]=Generalized_Eigen(K,M,1);

dampv=ones(ndof,1)*damp;

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t=F(:,1);

F(:,1)=[];

size(F)

[x,v,a,nx,nv,na]=ramp_invariant_force(ModeShapes,F,ndof,sys_dof,omegan,dampv,dt);


% nx=zeros(NT,ndof);
% nv=zeros(NT,ndof);
% na=zeros(NT,ndof);

c1=zeros(num_modes,1);
c2=zeros(num_modes,1);

for j=1:num_modes
    c1(j)=2*dampv(j)/omegan(j);
    c2(j)=1/omegan(j)^2;
end    


if(imethod==2)
    
    if(ndof<sys_dof)

        Kinv=pinv(K);
     
        for i=1:NT
            
            term=zeros(sys_dof,1);
            
            for k=1:sys_dof
                
                for j=1:num_modes

                    term(k)=term(k)+ModeShapes(k,j)*( c1(j)*nv(i,j) + c2(j)*na(i,j) );
               
                end
                   
            end
            
            yyy=Kinv*F(i,:)' -term;
            x(i,:)=yyy';  
            
        end 
        
        ndof=sys_dof;

        %  Velocity & acceleration remains the same per        
        %  Optimization of Large Structural Systems edited by George I. N. Rozvany        
        
    else
       
        msgbox('Mode Acceleration not performed since number of retained modes = total dof ');
        
    end
    
end

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
mdof_plot_legend_3iu(t,x,v,a,sys_dof,iu);

%%%% Leave after plot


clear acceleration;
clear velocity;
clear displacement;
%
if(iu==1)
    a=a/386;
end
%
acceleration=[t a];
velocity=[t v];
displacement=[t x];



setappdata(0,'acceleration',acceleration);
setappdata(0,'velocity',velocity);
setappdata(0,'displacement',displacement);



% --- Executes on button press in pushbutton_return.
function pushbutton_return_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_return (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(ri_mdof_modal_arbitrary_force);


% --- Executes on selection change in listbox_units.
function listbox_units_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_units (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_units contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_units


% --- Executes during object creation, after setting all properties.
function listbox_units_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_units (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_damping_Callback(hObject, eventdata, handles)
% hObject    handle to edit_damping (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_damping as text
%        str2double(get(hObject,'String')) returns contents of edit_damping as a double


% --- Executes during object creation, after setting all properties.
function edit_damping_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_damping (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_mass_Callback(hObject, eventdata, handles)
% hObject    handle to edit_mass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_mass as text
%        str2double(get(hObject,'String')) returns contents of edit_mass as a double


% --- Executes during object creation, after setting all properties.
function edit_mass_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_mass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_stiffness_Callback(hObject, eventdata, handles)
% hObject    handle to edit_stiffness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_stiffness as text
%        str2double(get(hObject,'String')) returns contents of edit_stiffness as a double


% --- Executes during object creation, after setting all properties.
function edit_stiffness_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_stiffness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_num_modes_Callback(hObject, eventdata, handles)
% hObject    handle to edit_num_modes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_num_modes as text
%        str2double(get(hObject,'String')) returns contents of edit_num_modes as a double


% --- Executes during object creation, after setting all properties.
function edit_num_modes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_num_modes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_force_Callback(hObject, eventdata, handles)
% hObject    handle to edit_force (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_force as text
%        str2double(get(hObject,'String')) returns contents of edit_force as a double


% --- Executes during object creation, after setting all properties.
function edit_force_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_force (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_save.
function pushbutton_save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

n=get(handles.listbox_output,'Value');

if(n==1)
    data=getappdata(0,'acceleration');
end
if(n==2)
    data=getappdata(0,'velocity');
end
if(n==3)
    data=getappdata(0,'displacement');
end

output_array=get(handles.edit_output_array,'String');
assignin('base', output_array, data);

h = msgbox('Save Complete'); 


% --- Executes on selection change in listbox_output.
function listbox_output_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_output contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_output


% --- Executes during object creation, after setting all properties.
function listbox_output_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_output_array_Callback(hObject, eventdata, handles)
% hObject    handle to edit_output_array (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_output_array as text
%        str2double(get(hObject,'String')) returns contents of edit_output_array as a double


% --- Executes during object creation, after setting all properties.
function edit_output_array_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_output_array (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_method.
function listbox_method_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_method contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_method


% --- Executes during object creation, after setting all properties.
function listbox_method_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end