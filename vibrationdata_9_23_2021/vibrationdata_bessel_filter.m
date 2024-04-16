function varargout = vibrationdata_bessel_filter(varargin)
% VIBRATIONDATA_BESSEL_FILTER MATLAB code for vibrationdata_bessel_filter.fig
%      VIBRATIONDATA_BESSEL_FILTER, by itself, creates a new VIBRATIONDATA_BESSEL_FILTER or raises the existing
%      singleton*.
%
%      H = VIBRATIONDATA_BESSEL_FILTER returns the handle to a new VIBRATIONDATA_BESSEL_FILTER or the handle to
%      the existing singleton*.
%
%      VIBRATIONDATA_BESSEL_FILTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIBRATIONDATA_BESSEL_FILTER.M with the given input arguments.
%
%      VIBRATIONDATA_BESSEL_FILTER('Property','Value',...) creates a new VIBRATIONDATA_BESSEL_FILTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before vibrationdata_bessel_filter_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to vibrationdata_bessel_filter_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help vibrationdata_bessel_filter

% Last Modified by GUIDE v2.5 21-Oct-2016 16:13:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @vibrationdata_bessel_filter_OpeningFcn, ...
                   'gui_OutputFcn',  @vibrationdata_bessel_filter_OutputFcn, ...
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


% --- Executes just before vibrationdata_bessel_filter is made visible.
function vibrationdata_bessel_filter_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to vibrationdata_bessel_filter (see VARARGIN)

% Choose default command line output for vibrationdata_bessel_filter
handles.output = hObject;

set(handles.edit_results,'Enable','off');

set(handles.listbox_method,'Value',1);


set(handles.pushbutton_save,'Enable','off');
set(handles.edit_output_array,'Enable','off');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes vibrationdata_bessel_filter wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = vibrationdata_bessel_filter_OutputFcn(hObject, eventdata, handles) 
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

set(handles.edit_results,'Enable','on');

YS=get(handles.edit_ylabel,'String');

k=get(handles.listbox_method,'Value');
 
if(k==1)
  FS=get(handles.edit_input_array,'String');
  THM=evalin('base',FS);   
else
  THM=getappdata(0,'THM');
end

y=double(THM(:,2));

t=THM(:,1);

rms_input=sqrt(mean(y)^2+std(y)^2);

n=length(y);

dur=THM(n,1)-THM(1,1);

dt=dur/(n-1);

fc=str2num(get(handles.edit_fc,'String'));

if isempty(fc)
   warndlg('Enter Lowpass Frequency'); 
   return;
end

scale=1.3617;
OM=tan(pi*fc*dt/scale);
%
OM2=OM^2;
%
den=1+3*OM+3*OM2;
%
b0=3*OM2/den;
b1=2*b0;
b2=b0;
%
a1=2*(-1+3*OM2)/den;
a2=(1-3*OM+3*OM2)/den;
%
out1=sprintf('\n OM=%8.4g ',OM);
out2=sprintf(' b0=%8.4g  b1=%8.4g  b2=%8.4g ',b0,b1,b2);
out3=sprintf(' a1=%8.4g  a2=%8.4g ',a1,a2);
disp(out1);
disp(out2);
disp(out3);
%
forward=[ b0,  b1,  b2 ];
back   =[     1, a1, a2 ];
yf=filter(forward,back,y);
%

t=fix_size(t);
yf=fix_size(yf);

filtered_data=[t yf];
%
figure(1);
plot(t,THM(:,2));
title('Input Time History');
xlabel(' Time(sec) ')
ylabel(YS)
grid on;

figure(2);
plot(t,yf);
out1=sprintf('Lowpass Filtered Time History %g Hz',fc);
title(out1);
xlabel(' Time(sec) ')
ylabel(YS)
grid on;


res=THM(:,2)-yf;
residual_data=[t res];

figure(3);
plot(t,res);
out1=sprintf('Residual Time History %g Hz',fc);
title(out1);
xlabel(' Time(sec) ')
ylabel(YS)
grid on;

rms_output=sqrt(mean(yf)^2+std(yf)^2);

setappdata(0,'filtered_data',filtered_data);   
setappdata(0,'residual_data',residual_data);  

s1=sprintf(' Input');
s2=sprintf('\n %8.4g RMS',rms_input);
s3=sprintf('\n\n Filtered Data');
s4=sprintf('\n %8.4g RMS',rms_output);

string=strcat(s1,s2,s3,s4);

set(handles.edit_results,'String',string)

set(handles.pushbutton_save,'Enable','on');
set(handles.edit_output_array,'Enable','on');

% --- Executes on button press in pushbutton_return.
function pushbutton_return_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_return (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(vibrationdata_bessel_filter)

% --- Executes on selection change in listbox_method.
function listbox_method_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_method contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_method
n=get(hObject,'Value');

set(handles.text_input_array_name,'Visible','on');
set(handles.edit_input_array,'Visible','on');

set(handles.pushbutton_save,'Enable','off');
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



function edit_ylabel_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ylabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ylabel as text
%        str2double(get(hObject,'String')) returns contents of edit_ylabel as a double


% --- Executes during object creation, after setting all properties.
function edit_ylabel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ylabel (see GCBO)
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

n=get(handles.listbox_type,'Value');

if(n==1)
    data=getappdata(0,'filtered_data');
else
    data=getappdata(0,'residual_data');    
end
    
output_name=get(handles.edit_output_array,'String');
assignin('base', output_name, data);

h = msgbox('Save Complete'); 


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



function edit_fc_Callback(hObject, eventdata, handles)
% hObject    handle to edit_fc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_fc as text
%        str2double(get(hObject,'String')) returns contents of edit_fc as a double
set(handles.edit_results,'Enable','off');

% --- Executes during object creation, after setting all properties.
function edit_fc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function text12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function edit_results_Callback(hObject, eventdata, handles)
% hObject    handle to edit_results (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_results as text
%        str2double(get(hObject,'String')) returns contents of edit_results as a double


% --- Executes during object creation, after setting all properties.
function edit_results_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_results (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_transfer.
function pushbutton_transfer_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_transfer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fc=str2num(get(handles.edit_fc,'String'));
L=2;  % two-poles, second order

dzero=(factorial(2*L))/((2^L)*factorial(L));
%
nn=5000;
df=fc/100;
%
fH=zeros(nn,1);
H=zeros(nn,1);
B=zeros(L+1,1);
%
for i=1:nn
    ff=i*df;
    fH(i)=ff;
    s=complex(0,(ff/fc));
    B(1)=1;
    B(2)=s+1;
    for j=3:L+1
        B(j)=(2*L-1)*B(j-1) + s^2*B(j-2);
    end
    H(i)=dzero/B(L+1);
%  
end
%
lpf=fc;
out1=sprintf(' Bessel Lowpass Filter order=%d fc=%g Hz ',L,lpf);
%
%  Phase
%
figure(3);
subplot(3,1,1);

theta=angle(H)*180/pi;

plot(fH,theta);
ylim([-180,0]);
set(gca,'MinorGridLineStyle','none','GridLineStyle',':','XScale','log',...
                    'YScale','lin','ytick',[-180,-90,0]);
set(gca,'MinorGridLineStyle','none','GridLineStyle',':','XScale','log','YScale','lin')
ylabel(' Phase(deg) ');
xlabel(' Frequency (Hz) ');
title(out1);
grid on;
%
%  Magnitude
%
subplot(3,1,[2 3]);
plot(fH,abs(H));
set(gca,'MinorGridLineStyle','none','GridLineStyle',':','XScale','log','YScale','log');
grid on;
xlabel(' Frequency (Hz) ');
ylabel(' Magnitude ');

%
figure(4);
plot(fH,real(H),fH,imag(H));
set(gca,'MinorGridLineStyle','none','GridLineStyle',':','XScale','log','YScale','lin')
ylabel(' Amplitude ');
xlabel(' Frequency (Hz) ');
title(out1);
legend ('real','imaginary');  
grid on;


% --- Executes on selection change in listbox_type.
function listbox_type_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_type


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
