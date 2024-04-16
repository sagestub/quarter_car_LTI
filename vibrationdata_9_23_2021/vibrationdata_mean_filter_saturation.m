function varargout = vibrationdata_mean_filter_saturation(varargin)
% VIBRATIONDATA_MEAN_FILTER_SATURATION MATLAB code for vibrationdata_mean_filter_saturation.fig
%      VIBRATIONDATA_MEAN_FILTER_SATURATION, by itself, creates a new VIBRATIONDATA_MEAN_FILTER_SATURATION or raises the existing
%      singleton*.
%
%      H = VIBRATIONDATA_MEAN_FILTER_SATURATION returns the handle to a new VIBRATIONDATA_MEAN_FILTER_SATURATION or the handle to
%      the existing singleton*.
%
%      VIBRATIONDATA_MEAN_FILTER_SATURATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIBRATIONDATA_MEAN_FILTER_SATURATION.M with the given input arguments.
%
%      VIBRATIONDATA_MEAN_FILTER_SATURATION('Property','Value',...) creates a new VIBRATIONDATA_MEAN_FILTER_SATURATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before vibrationdata_mean_filter_saturation_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to vibrationdata_mean_filter_saturation_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help vibrationdata_mean_filter_saturation

% Last Modified by GUIDE v2.5 06-Aug-2014 14:31:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @vibrationdata_mean_filter_saturation_OpeningFcn, ...
                   'gui_OutputFcn',  @vibrationdata_mean_filter_saturation_OutputFcn, ...
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


% --- Executes just before vibrationdata_mean_filter_saturation is made visible.
function vibrationdata_mean_filter_saturation_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to vibrationdata_mean_filter_saturation (see VARARGIN)

% Choose default command line output for vibrationdata_mean_filter_saturation
handles.output = hObject;

set(handles.edit_Q,'String','10');

set(handles.listbox_method,'Value',1);
set(handles.listbox_saturation_removal,'Value',1);
set(handles.listbox_passes,'Value',1);

set(handles.pushbutton_save,'Enable','off');
set(handles.edit_output_array,'Enable','off');


set(handles.uipanel_Saturation_Removal,'Visible','off');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes vibrationdata_mean_filter_saturation wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = vibrationdata_mean_filter_saturation_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_return.
function pushbutton_return_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_return (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(vibrationdata_mean_filter_saturation);


% --- Executes on button press in pushbutton_calculate.
function pushbutton_calculate_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_calculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fig_num=1;

iu=get(handles.listbox_unit,'Value');

k=get(handles.listbox_method,'Value');

iflag=1;
 
if(k==1)
     try  
         FS=get(handles.edit_input_array,'String');
         THM=evalin('base',FS);  
         iflag=1;
     catch
         iflag=0; 
         warndlg('Input Array does not exist.  Try again.')
     end
else
  THM=getappdata(0,'THM');
end
 
if(iflag==0)
    return;
end 

setappdata(0,'t',THM(:,1)); 

%
figure(fig_num);
fig_num=fig_num+1;
plot(THM(:,1),THM(:,2));
title('Base Input Acceleration Time History');
xlabel(' Time(sec) ')
if(iu==1 || iu==2)
    ylabel('Accel (G)')
else
    ylabel('Accel (m/sec^2)')    
end
grid on;



y=double(THM(:,2));

n=length(y);
last=n;

dur=THM(n,1)-THM(1,1);

dt=dur/(n-1);
sr=1/dt;

Q=str2num(get(handles.edit_Q,'String'));
setappdata(0,'Q',Q);  

fstart=str2num(get(handles.edit_start_frequency,'String'));

damp=1/(2*Q);

fn(1)=fstart;

fmax=sr/8;

oct=2^(1/12);

num=1;

while(1)
    
  num=num+1;
  
  fn(num)=fn(num-1)*oct;
  
  if(fn(num)>fmax)
      break;
  end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Initialize coefficients
%
[a1,a2,b1,b2,b3,rd_a1,rd_a2,rd_b1,rd_b2,rd_b3]=...
                                   srs_coefficients(fn,damp,dt);
%
%
%  SRS engine
%
clear length;
b=y;
nspec=length(fn);
good_freq=1.0e+20;
slope_goal=6.;
max_ref=zeros(nspec,1);
min_ref=zeros(nspec,1);
%
[xmax,xmin,dmax,slope_p,slope_n,ref_error,ppp,nnn,nnnoct]=...
    vibrationdata_mr_SRS(dt,nspec,last,b,a1,a2,b1,b2,b3,fn,good_freq,max_ref,min_ref,slope_goal);
%
%
max_ref=xmax;
min_ref=xmin;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

%
figure(fig_num);
fig_num=fig_num+1;
plot(fn,xmax,fn,xmin);
legend('positive','negative');
out1=sprintf('Acceleration SRS Q=%g',Q);
title(out1);
set(gca,'MinorGridLineStyle',':','GridLineStyle',':','XScale','log','YScale','log');
xlabel('Natural Frequency (Hz)');
if(iu==1 || iu==2)
    ylabel('Peak Accel (G)')
else
    ylabel('Peak Accel (m/sec^2)')    
end
plot_fmax=str2num(get(handles.edit_plot_fmax,'String'));
xlim([fn(1) plot_fmax]);
grid on;

fn=fix_size(fn);
xmax=fix_size(xmax);
xmin=fix_size(xmin);

setappdata(0,'raw_srs',[fn xmax xmin]);

m=get(handles.listbox_saturation_removal,'Value');

if(m==1)
    
    setappdata(0,'fn',fn)    
    setappdata(0,'damp',damp);    
    setappdata(0,'dt',dt); 
    setappdata(0,'a1',a1); 
    setappdata(0,'a2',a2);    
    setappdata(0,'b1',b1); 
    setappdata(0,'b2',b2);    
    setappdata(0,'b3',b3);  
    
    setappdata(0,'nspec',nspec); 
    setappdata(0,'last',last);     
    setappdata(0,'b_raw',b);
    
    setappdata(0,'max_ref',max_ref);  
    setappdata(0,'min_ref',min_ref);
    
    setappdata(0,'fig_num',fig_num);
    
    
    set(handles.uipanel_Saturation_Removal,'Visible','on');

    for j=nspec:-1:2    
%
%        printf(" f=%8.4g Hz  xmax=%8.4g  xmin=%8.4g \n",f(j),xmax(j),xmin(j)); 
%
        if( xmax(j)>0 && abs(xmin(j))>0)
%       
            aaa=abs(20.*log10(  xmax(j)/abs(xmin(j)  ) ) );
            bbb=abs(20.*log10(xmax(j-1)/abs(xmin(j-1)) ) );
%
%           printf("      aaa=%8.4g  bbb=%8.4g \n",aaa,bbb); 
%
            if(aaa>=4. && bbb>=4.)
%
                out1=sprintf('%8.4g ',fn(j-1));
                set(handles.edit_good_frequency,'String',out1);
                break;
            end
        end
    end
    
    
%
end


% --- Executes on button press in pushbutton_save.
function pushbutton_save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

n=get(handles.listbox_output_type,'Value');

if(n==1)
    data=getappdata(0,'raw_srs'); 
end
if(n==2)
    data=getappdata(0,'cleaned_srs');
end
if(n==3)
    data=getappdata(0,'cleaned_th'); 
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


% --- Executes on selection change in listbox_output_type.
function listbox_output_type_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_output_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_output_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_output_type


% --- Executes during object creation, after setting all properties.
function listbox_output_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_output_type (see GCBO)
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


% --- Executes on selection change in listbox_method.
function listbox_method_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_method contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_method

n=get(hObject,'Value');

set(handles.pushbutton_save,'Enable','off');
set(handles.edit_input_array,'String',' ');

set(handles.text_input_array_name,'Visible','on');
set(handles.edit_input_array,'Visible','on');

if(n==1)
   set(handles.edit_input_array,'enable','on') 
else
   set(handles.edit_input_array,'enable','off')
   set(handles.text_input_array_name,'Visible','off');
   set(handles.edit_input_array,'Visible','off');
   
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


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in listbox_unit.
function listbox_unit_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_unit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_unit contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_unit


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



function edit_start_frequency_Callback(hObject, eventdata, handles)
% hObject    handle to edit_start_frequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_start_frequency as text
%        str2double(get(hObject,'String')) returns contents of edit_start_frequency as a double


% --- Executes during object creation, after setting all properties.
function edit_start_frequency_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_start_frequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Q_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Q (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Q as text
%        str2double(get(hObject,'String')) returns contents of edit_Q as a double


% --- Executes during object creation, after setting all properties.
function edit_Q_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Q (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_plot_fmax_Callback(hObject, eventdata, handles)
% hObject    handle to edit_plot_fmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_plot_fmax as text
%        str2double(get(hObject,'String')) returns contents of edit_plot_fmax as a double


% --- Executes during object creation, after setting all properties.
function edit_plot_fmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_plot_fmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_saturation_removal.
function listbox_saturation_removal_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_saturation_removal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_saturation_removal contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_saturation_removal


% --- Executes during object creation, after setting all properties.
function listbox_saturation_removal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_saturation_removal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_slope_Callback(hObject, eventdata, handles)
% hObject    handle to edit_slope (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_slope as text
%        str2double(get(hObject,'String')) returns contents of edit_slope as a double


% --- Executes during object creation, after setting all properties.
function edit_slope_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_slope (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_passes.
function listbox_passes_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_passes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_passes contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_passes


% --- Executes during object creation, after setting all properties.
function listbox_passes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_passes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_calculate_saturation_removal.
function pushbutton_calculate_saturation_removal_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_calculate_saturation_removal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

    good_freq=str2num(get(handles.edit_good_frequency,'String'));

    fig_num=getappdata(0,'fig_num');

    t=getappdata(0,'t'); 
    fn=getappdata(0,'fn');    
    damp=getappdata(0,'damp');  
    Q=getappdata(0,'Q');  
    dt=getappdata(0,'dt'); 
    a1=getappdata(0,'a1'); 
    a2=getappdata(0,'a2');    
    b1=getappdata(0,'b1'); 
    b2=getappdata(0,'b2');    
    b3=getappdata(0,'b3');  
    
    nspec=getappdata(0,'nspec'); 
    last=getappdata(0,'last');     
    b=getappdata(0,'b_raw');
    a=b;
    clear b;
    
    max_ref=getappdata(0,'max_ref');  
    min_ref=getappdata(0,'min_ref');
    
    slope_goal=str2num(get(handles.edit_slope,'String'));
     
    max_window=str2num(get(handles.edit_window_size,'String'));
    
    max_passes=get(handles.listbox_passes,'value');
    
    if(max_passes>10000)
		max_passes=9999;
    end
    if(max_window>10000)
		max_window=9999;
    end    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

record=1.0e+90;
%
disp(' ');
%
nmp = 1;
%
ahold=a;
%
nt=max_passes*max_window;
%
progressbar;
for ix=1:max_passes    
for iq=3:2:max_window
%
    progressbar(iq*ix/nt);
%			 
    window_size=iq;
%
    number_passes=ix;
%      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
    w=window_size;
	k=fix(double(w-1)/2.);
	np=number_passes;
	npr_hold=np; % number of passes
	kr_hold=k;   % window size
%
    mf=zeros(last,1);
    
    a=ahold;
    
	for m=1:np	
		for i=1:last
%
			ave=0.;
			n=0;
%
			for j=(i-k):(i+k)
				if(j>=1 && j<=last )
					ave=ave+a(j);
					n=n+1;
                end		
            end
			if(n>1)			
				mf(i)=mf(i)+ave/double(n);
            end
        end
        a=ahold-mf;
    end
%
    b=a;
    ams=std(a);
    bms=std(b);
%
	arms=sqrt(ams/double(last));
	brms=sqrt(bms/double(last));
%
	if(arms<1.0e-20 || brms<1.0e-20 )
%
		disp('  Error 2 ');
		out1=sprintf('  arms=%8.4g  brms=%8.4g \n\n',arms,brms);
        disp(out1);
%
    end
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  SRS engine
%
%
[xmax,xmin,dmax,slope_p,slope_n,ref_error,ppp,nnn,nnnoct]=...
    vibrationdata_mr_SRS(dt,nspec,last,b,a1,a2,b1,b2,b3,fn,good_freq,max_ref,min_ref,slope_goal);
%

	dbc=dmax;   % want to minimize
%
	slope_pc=slope_p;
	slope_nc=slope_n;
%
	total_db = ref_error + dbc + ( slope_pc + slope_nc );
%
	if( total_db<record )
%
		record=total_db;
%
        mfr=mf;
        br=b;
%
		ps=xmax;
		pn=abs(xmin);
%
		npr=number_passes;
		kr=window_size;
%
		slope_pos=20.*log10(ppp)/nnnoct;
		slope_neg=20.*log10(nnn)/nnnoct;
%
 		out1=sprintf(' %ld  err=%8.4g dB  np=%ld  ws=%ld  +slope=%6.3g dB  -slope=%6.3g dB ',iq,record,number_passes,window_size,slope_pos,slope_neg);
        disp(out1);
    end
%
end
end    
%
pause(0.5);
progressbar(1);


out1=sprintf(' \n Best case: \n\n  Window size=%ld  Passes=%ld \n',kr,npr);   
disp(out1);
%    
clear n;
clear v;
%
n=last;
v(1)=b(1)*dt/2;
for(i=2:(n-1))
    v(i)=v(i-1)+b(i)*dt;
end
v(n)=v(n-1)+b(n)*dt/2;
%
t=fix_size(t);
b=fix_size(b);
mfr=fix_size(mfr);
%
saturation=[t mfr];
cleaned=[t b];

%
setappdata(0,'cleaned_th',cleaned);
%
iu=get(handles.listbox_unit,'Value');
%
figure(fig_num);
fig_num=fig_num+1;
plot(t,mfr)
grid on;
title(' Saturation Estimate ');
xlabel(' Time(sec) ');

if(iu==1 || iu==2)
    ylabel('Accel (G)')
else
    ylabel('Accel (m/sec^2)')    
end

%
figure(fig_num);
fig_num=fig_num+1;
plot(t,b)
grid on;
title(' Acceleration Cleaned Signal ');
xlabel(' Time(sec) ');
if(iu==1 || iu==2)
    ylabel('Accel (G)')
else
    ylabel('Accel (m/sec^2)')    
end

%
if(iu==1)
    v=v*386;
end
if(iu==2)
    v=v*9.81;
end
%
figure(fig_num);
fig_num=fig_num+1;
plot(t,v)
grid on;
title(' Velocity Cleaned Signal ');
xlabel(' Time(sec) ');
%
if(iu==1)
    ylabel('Vel (in/sec)')
else
    ylabel('Vel (m/sec)')    
end

%
figure(fig_num);
fig_num=fig_num+1;


sz=max(size(xmax));
fff=fn(1:sz);

fff=fix_size(fff);
ps=fix_size(ps);
pn=fix_size(pn);

setappdata(0,'cleaned_srs',[fff ps pn]);

plot(fff,ps,fff,abs(pn),'-.');
legend ('positive','negative');
out1=sprintf(' SRS Q=%g    Cleaned signal ',Q);   
title(out1);
xlabel(' Natural Frequency(sec) ');
if(iu==1 || iu==2)
    ylabel('Peak Accel (G)')
else
    ylabel('Peak Accel (m/sec^2)')    
end

grid;
set(gca,'MinorGridLineStyle',':','GridLineStyle',':','XScale','log','YScale','log');
%

set(handles.pushbutton_save,'Enable','on');
set(handles.edit_output_array,'Enable','on');


function edit_window_size_Callback(hObject, eventdata, handles)
% hObject    handle to edit_window_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_window_size as text
%        str2double(get(hObject,'String')) returns contents of edit_window_size as a double


% --- Executes during object creation, after setting all properties.
function edit_window_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_window_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_good_frequency_Callback(hObject, eventdata, handles)
% hObject    handle to edit_good_frequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_good_frequency as text
%        str2double(get(hObject,'String')) returns contents of edit_good_frequency as a double


% --- Executes during object creation, after setting all properties.
function edit_good_frequency_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_good_frequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
