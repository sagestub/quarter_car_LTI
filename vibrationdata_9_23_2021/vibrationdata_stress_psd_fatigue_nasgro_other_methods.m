function varargout = vibrationdata_stress_psd_fatigue_nasgro(varargin)
% VIBRATIONDATA_STRESS_PSD_FATIGUE_NASGRO MATLAB code for vibrationdata_stress_psd_fatigue_nasgro.fig
%      VIBRATIONDATA_STRESS_PSD_FATIGUE_NASGRO, by itself, creates a new VIBRATIONDATA_STRESS_PSD_FATIGUE_NASGRO or raises the existing
%      singleton*.
%
%      H = VIBRATIONDATA_STRESS_PSD_FATIGUE_NASGRO returns the handle to a new VIBRATIONDATA_STRESS_PSD_FATIGUE_NASGRO or the handle to
%      the existing singleton*.
%
%      VIBRATIONDATA_STRESS_PSD_FATIGUE_NASGRO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIBRATIONDATA_STRESS_PSD_FATIGUE_NASGRO.M with the given input arguments.
%
%      VIBRATIONDATA_STRESS_PSD_FATIGUE_NASGRO('Property','Value',...) creates a new VIBRATIONDATA_STRESS_PSD_FATIGUE_NASGRO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before vibrationdata_stress_psd_fatigue_nasgro_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to vibrationdata_stress_psd_fatigue_nasgro_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help vibrationdata_stress_psd_fatigue_nasgro

% Last Modified by GUIDE v2.5 29-Feb-2016 11:31:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @vibrationdata_stress_psd_fatigue_nasgro_OpeningFcn, ...
                   'gui_OutputFcn',  @vibrationdata_stress_psd_fatigue_nasgro_OutputFcn, ...
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


% --- Executes just before vibrationdata_stress_psd_fatigue_nasgro is made visible.
function vibrationdata_stress_psd_fatigue_nasgro_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to vibrationdata_stress_psd_fatigue_nasgro (see VARARGIN)

% Choose default command line output for vibrationdata_stress_psd_fatigue_nasgro
handles.output = hObject;

change_unit(hObject, eventdata, handles);
listbox_material_Callback(hObject, eventdata, handles)

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes vibrationdata_stress_psd_fatigue_nasgro wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = vibrationdata_stress_psd_fatigue_nasgro_OutputFcn(hObject, eventdata, handles) 
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

%%
 
A=str2num(get(handles.edit_A,'String'));
B=str2num(get(handles.edit_B,'String'));
C=str2num(get(handles.edit_C,'String'));
P=str2num(get(handles.edit_P,'String'));
 
if(B<0)
    warndlg(' B coefficient should be > 0 '); 
    return;
end

scf=str2num(get(handles.edit_scf,'String'));

FS=get(handles.edit_input_array,'String');
THM=evalin('base',FS); 

disp(' ');
disp('***************************************************');

out1=sprintf('\n PSD filename:  %s \n',FS);
disp(out1);

iu=get(handles.listbox_unit,'Value');

if(iu==1) % psi to ksi
    THM(:,2)=THM(:,2)/1000^2;
end
if(iu==3) % Pa to ksi
    THM(:,2)=THM(:,2)/( 6891^2  );
end
if(iu==4) % MPa to ksi
    THM(:,2)=THM(:,2)/(( 6891e-06  )^2);
end

[f,a,s,rms,df,THM,fi,ai]=fatigue_psd_check(THM,scf);

fig_num=1;

iu=get(handles.listbox_unit,'Value');

fmin=f(1);
fmax=max(f);

xlab='Frequency (Hz)';

if(iu==1)
    ylab='Stress (psi^2/Hz)';
    yu='psi';
    sss=sprintf(' %8.4g psi RMS ',rms);
end
if(iu==2)
    ylab='Stress (ksi^2/Hz)';
    yu='ksi';
    sss=sprintf(' %8.4g ksi RMS ',rms);    
end


if(iu==3)
    ylab='Stress (Pa^2/Hz)';
    yu='Pa';    
    sss=sprintf(' %8.4g Pa RMS ',rms);    
end
if(iu==4)
    ylab='Stress (MPa^2/Hz)';
    yu='MPa';    
    sss=sprintf(' %8.4g MPa RMS ',rms);    
end


t_string=sprintf('Stress PSD   %s',sss);
 
ylab='Stress (ksi^2/Hz)';
out1=sprintf('PSD unit:  %s ',ylab);
disp(out1);
 
 
 
[fig_num]=plot_PSD_function(fig_num,xlab,ylab,t_string,THM,fmin,fmax);

sigma_s=rms;

%
tau=str2num(get(handles.edit_duration,'String'));

 
 
clear length;
n=length(fi);
%

[m0,m1,m2,m4,M2m,Mkp2,m0p75,m1p5,vo,EP,alpha2,e,D1,D2,D3,Q,Rd]=...
                                     spectral_moments_nasgro(n,B,ai,fi,df);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Note that range = (peak-valley)
%
%  maxS is the assumed upper limit of the the histogram.
%
%  The value 6*grms is used as conservative estimate of the upper
%  range maxS for all cycles.  
%
%  The histogram will have 400 bins.  This number is chosen
%  via engineering judgement.
%
%    ds is the bin range width       
%     n is the number of bins
%     N is the cycle count per bin
%     S is the range level for each bin
%
rms=sqrt(m0);

maxS=8*rms;  
%
ds=maxS/400;
%
n=round(maxS/ds);
%
EP_tau=EP*tau;

clear R;
R=str2num(get(handles.edit_R,'String'));

%%%%
%
%   Dirlik
%

out1=sprintf(' D1=%8.4g  D2=%8.4g  D3=%8.4g   ',D1,D2,D3);
disp(out1);

[Dirlik_N,range]=Dirlik_nasgro_probability(n,m0,Q,D1,D2,D3,ds,Rd);

[Dirlik_damage]=Dirlik_nasgro(n,A,B,C,P,R,EP_tau,Dirlik_N,range);

%%%%
%
%  Narrowband
%

D1=0;
D2=0;
D3=1;

[Narrowband_N,range]=Dirlik_nasgro_probability(n,m0,Q,D1,D2,D3,ds,Rd);

[Narrowband_damage]=Dirlik_nasgro(n,A,B,C,P,R,EP_tau,Narrowband_N,range);

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

disp(' '); 
disp(' * * * * ');
disp(' ');
 
 
mlab=getappdata(0,'mlab');
 
out1=sprintf('\n Material: %s ',mlab);
disp(out1);
 
out1=sprintf('\n A=%g  B=%g  C=%g  P=%g ',A,B,C,P);
disp(out1);

out1=sprintf('\n m0=%8.4g \n m1=%8.4g \n m2=%8.4g \n m4=%8.4g ',m0,m1,m2,m4);
disp(out1);
 
out1=sprintf('\n Duration=%8.4g sec \n',tau);
disp(out1);
 

 
out1=sprintf('    v0+=%8.4g Hz, rate of zero upcrossings  ',vo);
disp(out1);
out1=sprintf('     EP=%8.4g Hz, expected peak rate ',EP);
disp(out1);
out1=sprintf(' alpha2=%8.4g ',alpha2);
disp(out1);
out1=sprintf('      e=%8.4g spectral width',e);
disp(out1);
 
 
a=vo*tau;
c=sqrt(2*log(a));
ps=c + 0.5772/c;
 
out1=sprintf('\n Stress Ratio R = %g',R);
disp(out1); 

out1=sprintf(' Stress Concentration Factor scf = %g',scf);
disp(out1);
 
out1=sprintf(' Overall Stress Level with scf = %8.4g ksi rms ',sqrt(m0));
disp(out1);
 
disp(' ');
out1=sprintf(' The maximum expected peak response with scf is: %8.4g ksi (%4.3g sigma) ',ps*rms,ps);
disp(out1);
 
out1=sprintf('\n Number of Rainflow Cycles=%8.4g ',EP_tau);
disp(out1);
 
disp('  ');
disp('    Cumulative Damage ');
out1=sprintf('\n     Dirlik = %8.4g ',Dirlik_damage);
  out2=sprintf(' Narrowband = %8.4g ',Narrowband_damage);
disp(out1);
disp(out2);

disp('  ');
disp('    Damage Rate (1/sec)');
out1=sprintf('\n     Dirlik = %8.4g ',Dirlik_damage/tau);
   out2=sprintf(' Narrowband = %8.4g ',Narrowband_damage/tau);
disp(out1);
disp(out2); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m=B

vo=sqrt(m2/m0);
vp=sqrt(m4/m2);
%
alpha=vo/vp;

alpha_1=m1/sqrt(m0*m2);
alpha_2=m2/sqrt(m0*m4);

delta=sqrt(1-alpha_1);
%
arg=0.5*m+1;
gf=gamma(arg);
%
beta=sqrt(m2*M2m/(m0*Mkp2));
betaem=beta^m;
%
lambda_oc=betaem/alpha;
%
a=0.926-0.033*m;
b=1.587*m-2.323;
ee=sqrt(1-alpha_2);
lambda_wl=a+(1-a)*(1-ee)^b;
%

lambda_ll=M2m^(m/2)/(vo*sigma_s^m);

alpha_0p75=m0p75/sqrt(m0*m1p5);

out5=sprintf('     Wirsching Light = %8.4g ',lambda_wl);
out6=sprintf('          Ortiz Chen = %8.4g ',lambda_oc);
out7=sprintf('      Lutes & Larsen = %8.4g ',lambda_ll);
out8=sprintf('        alpha_0p75^2 = %8.4g ',alpha_0p75^2);

disp(out5);
disp(out6);
disp(out7);
disp(out8);


q=(1-R)^P
1/q

msgbox('Results written to Matlab Command Window');
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%



% --- Executes on button press in pushbutton_return.
function pushbutton_return_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_return (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(vibrationdata_stress_psd_fatigue_nasgro);


% --- Executes on selection change in listbox_unit.
function listbox_unit_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_unit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_unit contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_unit
change_unit(hObject, eventdata, handles); 

% --- Executes during object creation, after setting all properties.
function listbox_unit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_unit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_input_array_Callback(hObject, eventdata, handles)
% hObject    handle to edit_input_array (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_input_array as text
%        str2double(get(hObject,'String')) returns contents of edit_input_array as a double


% --- Executes during object creation, after setting all properties.
function edit_input_array_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_input_array (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_material.
function listbox_material_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_material (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_material contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_material
change_unit(hObject, eventdata, handles);

n=get(handles.listbox_material,'Value');

set(handles.edit_A,'String','');
set(handles.edit_B,'String','');
set(handles.edit_C,'String','');
set(handles.edit_P,'String','');

mlab = 'other';

if(n==1) % Aluminum 6061-T6
    set(handles.edit_A,'String','20.68');
    set(handles.edit_B,'String','9.84');
    set(handles.edit_C,'String','0');
    set(handles.edit_P,'String','0.63');
    mlab = 'Aluminum 6061-T6';
end
if(n==2) % Aluminum 7075-T6
    set(handles.edit_A,'String','18.22');
    set(handles.edit_B,'String','7.77');
    set(handles.edit_C,'String','10.15');
    set(handles.edit_P,'String','0.62');    
    mlab = 'Aluminum 7075-T6';    
end

setappdata(0,'mlab',mlab);


% --- Executes during object creation, after setting all properties.
function listbox_material_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_material (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_m_Callback(hObject, eventdata, handles)
% hObject    handle to edit_m (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_m as text
%        str2double(get(hObject,'String')) returns contents of edit_m as a double


% --- Executes during object creation, after setting all properties.
function edit_m_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_m (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_A_Callback(hObject, eventdata, handles)
% hObject    handle to edit_A (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_A as text
%        str2double(get(hObject,'String')) returns contents of edit_A as a double


% --- Executes during object creation, after setting all properties.
function edit_A_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_A (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_duration_Callback(hObject, eventdata, handles)
% hObject    handle to edit_duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_duration as text
%        str2double(get(hObject,'String')) returns contents of edit_duration as a double


% --- Executes during object creation, after setting all properties.
function edit_duration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_scf_Callback(hObject, eventdata, handles)
% hObject    handle to edit_scf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_scf as text
%        str2double(get(hObject,'String')) returns contents of edit_scf as a double


% --- Executes during object creation, after setting all properties.
function edit_scf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_scf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function change_unit(hObject, eventdata, handles) 
%
%  fix here





% --- Executes on selection change in listbox_mean.
function listbox_mean_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_mean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_mean contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_mean







% --- Executes during object creation, after setting all properties.
function listbox_mean_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_mean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_mean_stress_Callback(hObject, eventdata, handles)
% hObject    handle to edit_mean_stress (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_mean_stress as text
%        str2double(get(hObject,'String')) returns contents of edit_mean_stress as a double


% --- Executes during object creation, after setting all properties.
function edit_mean_stress_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_mean_stress (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_stress_aux_Callback(hObject, eventdata, handles)
% hObject    handle to edit_stress_aux (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_stress_aux as text
%        str2double(get(hObject,'String')) returns contents of edit_stress_aux as a double


% --- Executes during object creation, after setting all properties.
function edit_stress_aux_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_stress_aux (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit_A (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_A as text
%        str2double(get(hObject,'String')) returns contents of edit_A as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_A (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_B_Callback(hObject, eventdata, handles)
% hObject    handle to edit_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_B as text
%        str2double(get(hObject,'String')) returns contents of edit_B as a double


% --- Executes during object creation, after setting all properties.
function edit_B_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_C_Callback(hObject, eventdata, handles)
% hObject    handle to edit_C (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_C as text
%        str2double(get(hObject,'String')) returns contents of edit_C as a double


% --- Executes during object creation, after setting all properties.
function edit_C_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_C (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_P_Callback(hObject, eventdata, handles)
% hObject    handle to edit_P (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_P as text
%        str2double(get(hObject,'String')) returns contents of edit_P as a double


% --- Executes during object creation, after setting all properties.
function edit_P_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_P (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_eq.
function pushbutton_eq_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_eq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

A = imread('nasgro_coefficients.jpg');
figure(555);
imshow(A,'border','tight','InitialMagnification',100);



function edit_R_Callback(hObject, eventdata, handles)
% hObject    handle to edit_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_R as text
%        str2double(get(hObject,'String')) returns contents of edit_R as a double


% --- Executes during object creation, after setting all properties.
function edit_R_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
