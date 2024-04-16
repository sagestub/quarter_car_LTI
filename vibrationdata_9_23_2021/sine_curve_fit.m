function varargout = sine_curve_fit(varargin)
% SINE_CURVE_FIT MATLAB code for sine_curve_fit.fig
%      SINE_CURVE_FIT, by itself, creates a new SINE_CURVE_FIT or raises the existing
%      singleton*.
%
%      H = SINE_CURVE_FIT returns the handle to a new SINE_CURVE_FIT or the handle to
%      the existing singleton*.
%
%      SINE_CURVE_FIT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SINE_CURVE_FIT.M with the given input arguments.
%
%      SINE_CURVE_FIT('Property','Value',...) creates a new SINE_CURVE_FIT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before sine_curve_fit_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to sine_curve_fit_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help sine_curve_fit

% Last Modified by GUIDE v2.5 16-Oct-2018 16:49:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @sine_curve_fit_OpeningFcn, ...
                   'gui_OutputFcn',  @sine_curve_fit_OutputFcn, ...
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


% --- Executes just before sine_curve_fit is made visible.
function sine_curve_fit_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to sine_curve_fit (see VARARGIN)

% Choose default command line output for sine_curve_fit
handles.output = hObject;

set(handles.listbox_method,'Value',1)
set(handles.listbox_type,'Value',1)

listbox_type_Callback(hObject, eventdata, handles);

set(handles.pushbutton_plot_TH,'Enable','on');
set(handles.pushbutton_calculate,'Enable','off');

set(handles.text_start_time,'Enable','off');
set(handles.text_end_time,'Enable','off');
set(handles.edit_start_time,'Enable','off');
set(handles.edit_end_time,'Enable','off');

set(handles.uipanel_save,'Visible','off');

set(handles.listbox_delay,'Value',1);


listbox_method_Callback(hObject, eventdata, handles);

% Find all windows of type figure, which have an empty FileName attribute.
allPlots = findall(0, 'Type', 'figure', 'FileName', []);
delete(allPlots);


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes sine_curve_fit wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = sine_curve_fit_OutputFcn(hObject, eventdata, handles) 
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

% Find all windows of type figure, which have an empty FileName attribute.
allPlots = findall(0, 'Type', 'figure', 'FileName', []);
delete(allPlots);

try
    close(handles.fig2);
catch
end
try
    close(handles.fig10);
catch
end

k=get(handles.listbox_type,'Value');

if(k==1)
    sinefind(hObject, eventdata, handles)
else
    sinefdam(hObject, eventdata, handles)
end

function sinefind(hObject, eventdata, handles)

THM=getappdata(0,'THM');
amp=double(THM(:,2));
 
n=length(amp);
dur=THM(n,1)-THM(1,1);
dt=dur/(n-1);

s_amp=THM(:,2);
s_tim=THM(:,1);

tmi=min(THM(:,1));
tmx=max(THM(:,1));

try
    ts=str2num(get(handles.edit_start_time,'String'));
catch
    ts=tmi;
end

try
    te=str2num(get(handles.edit_end_time,'String'));
catch
    te=tmx;
end


nfr=str2num(get(handles.edit_number,'String'));


%
if(ts>tmx)
    ts=tmi;
    te=tmx;
end
if(te<tmi)
    ts=tmi;
    te=tmx;
end
%

out1=sprintf(' ts=%8.4g  te=%8.4g  tmi=%8.4g  dt=%8.4g  ',ts,te,tmi,dt);
disp(out1);

n1=fix((ts-tmi)/dt);
n2=fix((te-tmi)/dt);
%
if(n1<1)
    n1=1;
end
%
if(n2>n)
    n2=n;
end
%
if(n1>n2)
    n2=n1;
end
%

out1=sprintf(' n1=%d  n2=%d  ',n1,n2);
disp(out1);
%

if(n1==n2)
    warndlg(' n1=n2 ');
    return;
end

%
tim=s_tim(n1:n2)';
t=tim;
%
amp=s_amp(n1:n2)';

mr=get(handles.listbox_remove_mean,'Value');

if(mr==1)
    amp=amp-mean(amp);
end

a=amp;
    
amp_orig=amp;

%
np=length(tim);


if(np<=0)
    warndlg(' Number of points is zero ');
    return;
end

%
out1 = sprintf(' np  = %d  \n',np);
disp(out1)
%
%
dur=max(t)-min(t);

out1=sprintf('\n dur = %8.4g sec \n',dur);
disp(out1);

if(dur<=1.0e-20)
    warndlg(' Duration = 0 ');
    return;
end





[syn,residual,Ar,Br,omeganr]=sine_find_function(dur,a,amp_orig,t,dt,nfr);

figure(999)
plot(tim,amp_orig,tim,syn);

%
tmax=1.0e+90;

LM=length(s_tim);


for i=1:LM
    terr=abs(s_tim(i)-t(1));
    if(terr<tmax)
        tmax=terr;
    end
end

s_tim=fix_size(s_tim);
s_amp=fix_size(s_amp);

t=fix_size(t);
syn=fix_size(syn);


tim=fix_size(tim);
residual=fix_size(residual);

handles.fig2=figure(2);
plot(tim,residual);
xlabel('Time (sec)');
hold on
title('Residual');


%
figure(3)
plot(s_tim,s_amp,'b');
xlabel('Time (sec)');
hold on
plot(t,syn,'r');
xlabel(' Time (sec) ');
legend ('Input Data','Synthesis');
tmin=tmi;
tmax=tmx;


max1=max(s_amp);
min1=min(s_amp);

max2=max(syn);
min2=min(syn);


ymax=1.5*max( [max1 max2]);
ymin=1.5*min( [min1 min2]);
axis([tmin,tmax,ymin,ymax]);
zoom on;
hold off;


synth_th=[t syn ];
residual_th=[tim residual];

setappdata(0,'synth_th',synth_th);
setappdata(0,'residual_th',residual_th);



tp=2*pi;
%
%      y=Ar*cos(omeganr*tt) + Br*sin(omeganr*tt); 
%

mag=zeros(nfr,1);
phase=zeros(nfr,1);


disp(' ');
disp(' Results');
disp(' ');
disp(' Case   Amplitude   fn(Hz)   Phase(rad)  ');
%

for ie=1:nfr
    mag(ie)=norm([Ar(ie) Br(ie)]);
    phase(ie)=atan2(Ar(ie),Br(ie));
    out4 = sprintf('  %d  %10.4f  %10.4f  %10.4f  ',ie,mag(ie),omeganr(ie)/tp,phase(ie));
    disp(out4)      
end    


omeganr=fix_size(omeganr);

aaa=[mag omeganr/tp phase];
bbb=sortrows(aaa,2);
setappdata(0,'curve_fit_table',bbb);


set(handles.uipanel_save,'Visible','on');


guidata(hObject, handles); % leave here

msgbox('Results Written to Matlab Command Window'); 


function sinefdam(hObject, eventdata, handles)

tic


THM=getappdata(0,'THM');


ndel=get(handles.listbox_delay,'Value');

tp=2*pi;

ts=str2num(get(handles.edit_start_time,'String'));
te=str2num(get(handles.edit_end_time,'String'));

nfr=str2num(get(handles.edit_number,'String'));


md1=str2num(get(handles.edit_max_damp,'String'));
md2=str2num(get(handles.edit_min_damp,'String'));


s_amp=THM(:,2);
s_tim=THM(:,1);


mr=get(handles.listbox_remove_mean,'Value');

if(mr==1)
    s_amp=s_amp-mean(s_amp);
end

%
n = length(s_amp); 
%
tmx=max(s_tim);
tmi=min(s_tim);
%
out1 = sprintf('\n  %d samples \n',n);
disp(out1)
%
dt=(tmx-tmi)/(n-1);
sr=1./dt;
%
out1 = sprintf(' SR  = %g samples/sec    dt = %g sec \n',sr,dt);
out3 = sprintf(' start  = %g sec    end = %g sec            ',tmi,tmx);
disp(out1)
disp(out3)
%
if(ts>tmx)
    ts=tmx;
end
if(te<tmi)
    ts=tmi;
end
%
n1=fix((ts-tmi)/dt);
n2=fix((te-tmi)/dt);
%
if(n1<1)
    n1=1;
end
if(n2>n)
    n2=n;
end    
%
if(n1>n2)
    temp=n1;
    n1=n2;
    n2=temp;
end
%
amp=s_amp(n1:n2)';
tim=s_tim(n1:n2)';

amp=amp-mean(amp);

amp_orig=amp;


np=length(tim);
n=np;
num2=np;

%
out1 = sprintf(' np  = %d  \n',np);
disp(out1)
%
%
t=tim;
av=amp;
%
clear THM;
%
dur=t(num2)-t(1);
%
dt=dur/(num2-1);


if(ndel==1)
    [syn,residual,Ar,Br,omeganr,dampr]=damped_sine_find_function(dur,av,amp_orig,t,dt,nfr,md1,md2);
else
    [syn,residual,Ar,Br,omeganr,dampr,delayr]=damped_sine_find_function_delay(dur,av,amp_orig,t,dt,nfr,md1,md2);    
end


% % % % % %

fx=omeganr/tp;

clear xx;
clear yy;
    
handles.fig10=figure(10);


sx1r=zeros(nfr,1);

for i=1:nfr
    
    sx1r(i)=sqrt( Ar(i)^2 + Br(i)^2 );
    
    xx=[fx(i) fx(i)];
    yy=[0  sx1r(i)];
    
    line(xx,yy,'Color','b');
end

xlabel('Frequency(Hz)');
ylabel('Amplitude');
grid on;
title('Damped Sine Curve-fit Spectra');


%

tmax=1.0e+90;

LM=length(s_tim);

for i=1:LM
    terr=abs(s_tim(i)-t(1));
    if(terr<tmax)
        tmax=terr;
    end
end

s_tim=fix_size(s_tim);
s_amp=fix_size(s_amp);

t=fix_size(t);
syn=fix_size(syn);

tim=fix_size(tim);
residual=fix_size(residual);
 


handles.fig2=figure(2);
plot(tim,residual);
xlabel('Time (sec)');
hold on
title('Residual')


%
figure(3);
plot(s_tim,s_amp,'b');
xlabel('Time (sec)');
hold on
plot(t,syn,'r');
xlabel(' Time (sec) ');
legend ('Input Data','Synthesis');
tmin=tmi;
tmax=tmx;

ymax=1.5*max(s_amp);
ymin=1.5*min(s_amp);
axis([tmin,tmax,ymin,ymax]);
zoom on;
hold off;



synth_th=[t syn ];
residual_th=[tim residual];

setappdata(0,'synth_th',synth_th);
setappdata(0,'residual_th',residual_th);


mag=zeros(nfr,1);
phase=zeros(nfr,1);

omeganr=fix_size(omeganr);


for ie=1:nfr    
    mag(ie)=norm([Ar(ie) Br(ie)]);
    phase(ie)=atan2(Ar(ie),Br(ie));
end        

disp(' ');
disp(' Results');
disp(' ');

if(ndel==1)
    
    disp(' Case   Amplitude   fn(Hz)   Phase(rad)  damp ');

    for ie=1:nfr
        out4 = sprintf('  %d  %10.4f  %10.4f  %10.4f  %8.4f',ie,mag(ie),omeganr(ie)/tp,phase(ie),dampr(ie));
        disp(out4)      
    end  
    
    aaa=[mag omeganr/tp phase dampr];
    
else
    
    disp(' Case   Amplitude   fn(Hz)   Phase(rad)  damp  delay(sec)');    
 
    for ie=1:nfr
        out4 = sprintf('  %d  %10.4f  %10.4f  %10.4f  %8.4f  %8.4f',ie,mag(ie),omeganr(ie)/tp,phase(ie),dampr(ie),delayr(ie));
        disp(out4)      
    end       
    
    aaa=[mag omeganr/tp phase dampr delayr];
    
end

bbb=sortrows(aaa,2);
setappdata(0,'curve_fit_table',bbb);



disp(' ');
disp(' Results Sorted by Frequency');
disp(' ');

if(ndel==1)

    disp(' Case   Amplitude   fn(Hz)   Phase(rad)  damp');

    for ie=1:nfr
        out4 = sprintf('  %d  %10.4f  %10.4f  %10.4f  %8.4f',ie,bbb(ie,1),bbb(ie,2),bbb(ie,3),bbb(ie,4));
        disp(out4)      
    end   
    
    aaa=[mag omeganr/tp phase dampr];
    
else
    
    disp(' Case   Amplitude   fn(Hz)   Phase(rad)  damp   delay(sec)');

    for ie=1:nfr
        out4 = sprintf('  %d  %10.4f  %10.4f  %10.4f  %8.4f  %8.4f',ie,bbb(ie,1),bbb(ie,2),bbb(ie,3),bbb(ie,4),bbb(ie,5));
        disp(out4)      
    end   
    
    aaa=[mag omeganr/tp phase dampr delayr];
    
end


bbb=sortrows(aaa,2);
setappdata(0,'curve_fit_table',bbb);


set(handles.uipanel_save,'Visible','on');

disp(' ');
toc

msgbox('Results Written to Matlab Command Window'); 


% --- Executes on button press in pushbutton_return.
function pushbutton_return_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_return (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(sine_curve_fit)


% --- Executes on selection change in listbox_method.
function listbox_method_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_method contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_method
n=get(handles.listbox_method,'Value');

set(handles.pushbutton_calculate,'Enable','off');

set(handles.text_input_array_name,'Visible','on');
set(handles.edit_input_array,'Visible','on');

set(handles.edit_input_array,'String',' ');

if(n==1)
   set(handles.edit_input_array,'enable','on') 
else
   set(handles.text_input_array_name,'Visible','off');
   set(handles.edit_input_array,'Visible','off');
    
   set(handles.edit_input_array,'enable','off')
   
   [filename, pathname] = uigetfile('*.*');
   filename = fullfile(pathname, filename);
   fid = fopen(filename,'r');
   THM = fscanf(fid,'%g %g',[2 inf]);
   THM=THM';
    
   setappdata(0,'THM',THM);
end

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


% --- Executes on selection change in listbox_type.
function listbox_type_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_type


n=get(handles.listbox_type,'Value');

set(handles.text_max_damp,'Visible','off');
set(handles.edit_max_damp,'Visible','off');

set(handles.text_min_damp,'Visible','off');
set(handles.edit_min_damp,'Visible','off');

set(handles.text_zd,'Visible','off');
set(handles.listbox_delay,'Visible','off');


if(n==1)
    set(handles.text_number,'String','Number of Sines');
else
    set(handles.text_number,'String','Number of Damped Sines'); 
    set(handles.text_max_damp,'Visible','on');
    set(handles.edit_max_damp,'Visible','on');    
    set(handles.text_min_damp,'Visible','on');
    set(handles.edit_min_damp,'Visible','on');
    
    set(handles.text_zd,'Visible','on');
    set(handles.listbox_delay,'Visible','on');    
end



% --- Executes during object creation, after setting all properties.
function listbox_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_plot_TH.
function pushbutton_plot_TH_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_plot_TH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

k=get(handles.listbox_method,'Value');
 
if(k==1)
  try  
    FS=get(handles.edit_input_array,'String');
    THM=evalin('base',FS);
  catch
    warndlg('File does not exist');  
    return  
  end  
else
  THM=getappdata(0,'THM');    
end


mr=get(handles.listbox_remove_mean,'Value');

if(mr==1)
    THM(:,2)=THM(:,2)-mean(THM(:,2));
end

setappdata(0,'THM',THM);

figure(1);
plot(THM(:,1),THM(:,2));
title('Time History');
xlabel('Time (sec)');
grid on;


set(handles.pushbutton_calculate,'Enable','on');

set(handles.text_start_time,'Enable','on');
set(handles.text_end_time,'Enable','on');
set(handles.edit_start_time,'Enable','on');
set(handles.edit_end_time,'Enable','on');


tmi=min(THM(:,1));
tmx=max(THM(:,1));

s1=sprintf(' %8.4g ',tmi);
s2=sprintf(' %8.4g ',tmx);


hhh=get(handles.edit_start_time,'String');

if( isempty(hhh) )
    set(handles.edit_start_time,'String',s1);
end

set(handles.edit_end_time,'String',s2);



function edit_number_Callback(hObject, eventdata, handles)
% hObject    handle to edit_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_number as text
%        str2double(get(hObject,'String')) returns contents of edit_number as a double


% --- Executes during object creation, after setting all properties.
function edit_number_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_trials_Callback(hObject, eventdata, handles)
% hObject    handle to edit_trials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_trials as text
%        str2double(get(hObject,'String')) returns contents of edit_trials as a double


% --- Executes during object creation, after setting all properties.
function edit_trials_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_trials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_start_time_Callback(hObject, eventdata, handles)
% hObject    handle to edit_start_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_start_time as text
%        str2double(get(hObject,'String')) returns contents of edit_start_time as a double


% --- Executes during object creation, after setting all properties.
function edit_start_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_start_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_end_time_Callback(hObject, eventdata, handles)
% hObject    handle to edit_end_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_end_time as text
%        str2double(get(hObject,'String')) returns contents of edit_end_time as a double


% --- Executes during object creation, after setting all properties.
function edit_end_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_end_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_export.
function pushbutton_export_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

n=get(handles.listbox_save,'Value');

if(n==1)
    data=getappdata(0,'synth_th');
end    
if(n==2)
    data=getappdata(0,'residual_th');    
end  
if(n==3)
    data=getappdata(0,'curve_fit_table');    
end  
    
    
output_name=get(handles.edit_output_array,'String');
assignin('base', output_name, data);

msgbox('Save Complete');


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


% --- Executes on selection change in listbox_save.
function listbox_save_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_save contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_save


% --- Executes during object creation, after setting all properties.
function listbox_save_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_max_damp_Callback(hObject, eventdata, handles)
% hObject    handle to edit_max_damp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_max_damp as text
%        str2double(get(hObject,'String')) returns contents of edit_max_damp as a double


% --- Executes during object creation, after setting all properties.
function edit_max_damp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_max_damp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function pushbutton_export_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on key press with focus on edit_input_array and none of its controls.
function edit_input_array_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to edit_input_array (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

set(handles.pushbutton_calculate,'Enable','off');


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over edit_input_array.
function edit_input_array_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to edit_input_array (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.pushbutton_calculate,'Enable','on');



function edit_min_damp_Callback(hObject, eventdata, handles)
% hObject    handle to edit_min_damp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_min_damp as text
%        str2double(get(hObject,'String')) returns contents of edit_min_damp as a double


% --- Executes during object creation, after setting all properties.
function edit_min_damp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_min_damp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_delay.
function listbox_delay_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_delay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_delay contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_delay


% --- Executes during object creation, after setting all properties.
function listbox_delay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_delay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function pushbutton_calculate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_calculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in pushbutton_close_all_plots.
function pushbutton_close_all_plots_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_close_all_plots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%% evalin('base', 'close all')
setappdata(0,'fig_num',1);


Figures = findobj( 'Type', 'Figure' , '-not' , 'Tag' , get( handles.output , 'Tag' ) );
NFigures = length( Figures );
for nFigures = 1 : NFigures
  close( Figures( nFigures ) );
end

msgbox('Plots Cleared'); 


% --- Executes on selection change in listbox_remove_mean.
function listbox_remove_mean_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_remove_mean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_remove_mean contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_remove_mean


% --- Executes during object creation, after setting all properties.
function listbox_remove_mean_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_remove_mean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
