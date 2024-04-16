function varargout = vibrationdata_psd(varargin)
% VIBRATIONDATA_PSD MATLAB code for vibrationdata_psd.fig
%      VIBRATIONDATA_PSD, by itself, creates a new VIBRATIONDATA_PSD or raises the existing
%      singleton*.
%
%      H = VIBRATIONDATA_PSD returns the handle to a new VIBRATIONDATA_PSD or the handle to
%      the existing singleton*.
%
%      VIBRATIONDATA_PSD('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIBRATIONDATA_PSD.M with the given input arguments.
%
%      VIBRATIONDATA_PSD('Property','Value',...) creates a new VIBRATIONDATA_PSD or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before vibrationdata_psd_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to vibrationdata_psd_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help vibrationdata_psd

% Last Modified by GUIDE v2.5 06-Feb-2015 17:15:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @vibrationdata_psd_OpeningFcn, ...
                   'gui_OutputFcn',  @vibrationdata_psd_OutputFcn, ...
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


% --- Executes just before vibrationdata_psd is made visible.
function vibrationdata_psd_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to vibrationdata_psd (see VARARGIN)

% Choose default command line output for vibrationdata_psd
handles.output = hObject;


    
set(handles.listbox_method,'Value',1);

set(handles.pushbutton_calculate,'Enable','off');
set(handles.pushbutton_view_options,'Enable','off');
set(handles.listbox_numrows,'Enable','off');
set(handles.listbox_numrows,'Visible','off');
set(handles.text_select_option,'Visible','off');

set(handles.edit_input_array,'Visible','on');
set(handles.text_input_array_name,'Visible','on');

set(handles.uitable_advise,'Visible','off');

set(handles.pushbutton_spec,'Visible','off');

set(handles.listbox_export_plot,'Value',2);

set(handles.listbox_crms,'Value',2);

listbox_export_plot_Callback(hObject, eventdata, handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes vibrationdata_psd wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = vibrationdata_psd_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


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


% --- Executes on selection change in listbox_method.
function listbox_method_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_method contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_method


set(handles.pushbutton_calculate,'Enable','off');
set(handles.pushbutton_view_options,'Enable','off');
set(handles.listbox_numrows,'Enable','off');
set(handles.listbox_numrows,'Visible','off');
set(handles.text_select_option,'Visible','off');
set(handles.uitable_advise,'Visible','off');

n=get(hObject,'Value');

set(handles.pushbutton_view_options,'Enable','on');


set(handles.edit_output_array,'Enable','off');


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



function edit_ylabel_input_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ylabel_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ylabel_input as text
%        str2double(get(hObject,'String')) returns contents of edit_ylabel_input as a double


% --- Executes during object creation, after setting all properties.
function edit_ylabel_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ylabel_input (see GCBO)
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


% --- Executes on selection change in listbox_mean_removal.
function listbox_mean_removal_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_mean_removal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_mean_removal contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_mean_removal


% --- Executes during object creation, after setting all properties.
function listbox_mean_removal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_mean_removal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_window.
function listbox_window_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_window contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_window


% --- Executes during object creation, after setting all properties.
function listbox_window_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_view_options.
function pushbutton_view_options_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_view_options (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB

set(handles.pushbutton_spec,'Visible','off');

set(handles.pushbutton_calculate,'Enable','on');
set(handles.listbox_numrows,'Enable','on');
set(handles.uitable_advise,'Visible','on');
set(handles.listbox_numrows,'Visible','on');
set(handles.text_select_option,'Visible','on');

set(handles.listbox_numrows,'String',' ');

k=get(handles.listbox_method,'Value');

if(k==1)
    try  
         FS=get(handles.edit_input_array,'String');
         THM=evalin('base',FS);  
         iflag=1;
    catch
         iflag=0; 
    end
    
else
    try  
       THM=getappdata(0,'THM');
       iflag=1;
    catch
       iflag=0; 
    end
end

if(iflag==0)
    warndlg('Input Array does not exist.  Try again.');
    return
end

t=THM(:,1);
y=double(THM(:,2));

n=length(t);

dur=t(n)-t(1);
dt=dur/(n-1);

NC=0;
for i=1:1000
%    
    nmp = 2^(i-1);
%   
    if(nmp <= n )
        ss(i) = 2^(i-1);
        seg(i) =n/ss(i);
        i_seg(i) = fix(seg(i));
        NC=NC+1;
    else
        break;
    end
end

disp(' ')
out4 = sprintf(' Number of   Samples per   Time per        df               ');
out5 = sprintf(' Segments     Segment      Segment(sec)   (Hz)     dof     ');
%
disp(out4)
disp(out5)
%
k=1;
for i=1:NC
    j=NC+1-i;
    if j>0
        if( i_seg(j)>0 )
%           str = int2str(i_seg(j));
            tseg=dt*ss(j);
            ddf=1./tseg;
            out4 = sprintf(' %8d  %8d    %11.5f    %9.4f   %d',i_seg(j),ss(j),tseg,ddf,2*i_seg(j));
            disp(out4)
            data(k,:)=[i_seg(j),ss(j),tseg,ddf,2*i_seg(j)];
            k=k+1;
        end
    end
    if(i==12)
        break;
    end
end
%

max_num_rows=k-1;


for i=1:max_num_rows
    handles.number(i)=i;
end

set(handles.listbox_numrows,'String',handles.number);


cn={'No. of Segments ','Samples/Segments','Time/Segment (sec)','df (Hz)','dof'};

set(handles.uitable_advise,'Data',data,'ColumnWidth','auto','ColumnName',cn);

setappdata(0,'advise_data',data);


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_return.
function pushbutton_return_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_return (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(vibrationdata_psd);

% --- Executes on button press in pushbutton_save.
function pushbutton_save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

data=getappdata(0,'PSD');

output_name=get(handles.edit_output_array,'String');
assignin('base', output_name, data);

h = msgbox('Save Complete'); 

set(handles.pushbutton_convert_octave,'Enable','on');



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


% --- Executes on selection change in listbox_numrows.
function listbox_numrows_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_numrows (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_numrows contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_numrows


% --- Executes during object creation, after setting all properties.
function listbox_numrows_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_numrows (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_calculate.
function pushbutton_calculate_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_calculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.pushbutton_spec,'Visible','off');

set(handles.pushbutton_save,'Enable','On');    
set(handles.edit_output_array,'Enable','On'); 

iflag=0;

k=get(handles.listbox_method,'Value');
 
if(k==1)

    iflag=1;

    try  
         FS=get(handles.edit_input_array,'String');
         THM=evalin('base',FS);  
         iflag=1;
    catch
         iflag=0; 
         warndlg('Input Array does not exist.  Try again.');
         return;
    end
    
else
    try  
        THM=getappdata(0,'THM');
        iflag=1;
    catch
         iflag=0; 
         warndlg('Input Array does not exist.  Try again.');
         return;
    end
end

if(iflag==0)
    return
end

t=THM(:,1);
amp=THM(:,2);

n=length(t);

dur=t(n)-t(1);

dt=dur/(n-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%(' Number of   Samples per   Time per        df               ');
%%(' Segments     Segment      Segment(sec)   (Hz)     dof     ');

advise_data=getappdata(0,'advise_data');

q=get(handles.listbox_numrows,'Value');

NW=advise_data(q,1);  % Number of Segments

mmm = 2^fix(log(n/NW)/log(2));
%
df=1./(mmm*dt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mr_choice=get(handles.listbox_mean_removal,'Value');

h_choice =get(handles.listbox_window,'Value');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[freq,full,rms]=psd_core(mmm,NW,mr_choice,h_choice,amp,df);

    disp(' ');
    out4 = sprintf(' Overall RMS = %10.3g ',rms);
    out5 = sprintf(' Three Sigma = %10.3g ',3*rms);
    disp(out4)
    disp(out5)
    disp(' ');
%

   clear power_spectral_density;
%
    freq=fix_size(freq);
    full=fix_size(full);
%
    power_spectral_density=[freq full]; 
%
    [~,fmax]=find_max(power_spectral_density);
%
    out5 = sprintf(' Peak occurs at %10.5g Hz \n',fmax);
    disp(out5)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

YS=get(handles.edit_ylabel_input,'String');


figure(1);
plot(THM(:,1),THM(:,2));
title('Input Time History');
xlabel(' Time(sec) ');
ylabel(YS);
grid on;

m=get(handles.listbox_type,'Value');


t_string=sprintf('Power Spectral Density %6.3g %sRMS Overall ',YS,rms);

k = strfind(YS,'/');


setappdata(0,'irms',rms);


if(m==1)
    out2=sprintf('Accel (%s^2/Hz)',YS);
    if( k>=1)
      out2=sprintf('Accel ((%s)^2/Hz)',YS);        
    end
end
if(m==2)
    out2=sprintf('Vel ((%s)^2/Hz)',YS);
    if( k>=1)
      out2=sprintf('Vel ((%s)^2/Hz)',YS);        
    end    
end
if(m==3)
    out2=sprintf('Disp (%s^2/Hz)',YS);
end
if(m==4)
    out2=sprintf('Force (%s^2/Hz)',YS);
end
if(m==5)
    out2=sprintf('Pressure (%s^2/Hz)',YS);
    if( k>=1)
      out2=sprintf('Pressure ((%s)^2/Hz)',YS);        
    end     
end
if(m==6)
    out2=sprintf('(%s^2/Hz)',YS);
end

fig_num=2;
xlab='Frequency (Hz)';
ylab=out2;

sfmin=get(handles.edit_fmin,'String');
sfmax=get(handles.edit_fmax,'String');


if isempty(sfmin)
    string=sprintf('%8.4g',df);
    set(handles.edit_fmin,'String',string); 
    sfmin=get(handles.edit_fmin,'String');    
end

if isempty(sfmax)
    sr=1/dt;
    nyf=sr/2;    
    string=sprintf('%8.4g',max(freq));
    set(handles.edit_fmax,'String',string);
    sfmax=get(handles.edit_fmax,'String');    
end

fmin=str2num(sfmin);
fmax=str2num(sfmax);


if(rms>=1000)

 t_string=...
    sprintf(' Power Spectral Density  Overall Level = %7.4g %sRMS ',rms,YS); 

else
    
 t_string=...
    sprintf(' Power Spectral Density  Overall Level = %6.3g %sRMS ',rms,YS);    

end

    [fig_num,h2]=plot_PSD_function(fig_num,xlab,ylab,t_string,power_spectral_density,fmin,fmax);

nps=get(handles.listbox_export_plot,'Value');    

if(nps==1)
    pname=get(handles.edit_plot_name','String');
    set(gca,'Fontsize',12);
    print(h2,pname,'-dpng','-r300');    
    disp(' ');
    disp(' Plot file exported to hard drive as: ');
    out1=sprintf(' %s.png',pname);
    disp(out1);    
end    
    
setappdata(0,'PSD',power_spectral_density);

setappdata(0,'xlab',xlab);
setappdata(0,'ylab',ylab);
setappdata(0,'t_string',t_string);
setappdata(0,'fmin',fmin);
setappdata(0,'fmax',fmax);

ncrms=get(handles.listbox_crms,'Value');

if(ncrms==1)

    [crms]=cumulative_rms(power_spectral_density);

    figure(3);
    plot(crms(:,1),crms(:,2));
    set(gca,'MinorGridLineStyle',':','GridLineStyle',':','XScale','log',...
                     'YScale','lin','XminorTick','off','YminorTick','off');
%
    set(gca,'XGrid','on','GridLineStyle',':');
    set(gca,'YGrid','on','GridLineStyle',':');
    xlabel('Frequency (Hz)');
    title('Cumulative RMS');
    out1=sprintf(' %sRMS ',YS); 
    ylabel(out1);
    xlim([fmin,fmax]);

end

set(handles.pushbutton_spec,'Visible','on');

setappdata(0,'fig_num',fig_num);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes on key press with focus on edit_input_array and none of its controls.
function edit_input_array_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to edit_input_array (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

set(handles.pushbutton_spec,'Visible','off');

set(handles.uitable_advise,'Visible','off');
set(handles.pushbutton_calculate,'Enable','off');
set(handles.pushbutton_view_options,'Enable','off');
set(handles.listbox_numrows,'Enable','off');
set(handles.listbox_numrows,'Visible','off');
set(handles.text_select_option,'Visible','off');


set(handles.pushbutton_view_options,'Enable','on');

set(handles.pushbutton_spec,'Visible','on');

% --- Executes on key press with focus on listbox_numrows and none of its controls.
function listbox_numrows_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to listbox_numrows (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)



function edit_fmax_Callback(hObject, eventdata, handles)
% hObject    handle to edit_fmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_fmax as text
%        str2double(get(hObject,'String')) returns contents of edit_fmax as a double


% --- Executes during object creation, after setting all properties.
function edit_fmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_fmin_Callback(hObject, eventdata, handles)
% hObject    handle to edit_fmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_fmin as text
%        str2double(get(hObject,'String')) returns contents of edit_fmin as a double


% --- Executes during object creation, after setting all properties.
function edit_fmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_convert_octave.
function pushbutton_convert_octave_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_convert_octave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.s= vibrationdata_psd_oct;   
set(handles.s,'Visible','on'); 


% --- Executes on button press in pushbutton_spec.
function pushbutton_spec_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_spec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


prompt = {'Enter spec array name'};
dlg_title = 'Input';
num_lines = 1;

try  
    temp=inputdlg(prompt,dlg_title,num_lines);
    FS=temp{:};
    THM=evalin('base',FS);  
catch 
    warndlg('Input Array does not exist.  Try again.');
    return;
end

[aslope,srms] = calculate_PSD_slopes(THM(:,1),THM(:,2));

ppp=getappdata(0,'PSD');
qqq=THM; 

irms=getappdata(0,'irms');

x_label=getappdata(0,'xlab');
y_label=getappdata(0,'ylab');
% t_string=getappdata(0,'t_string');
fmin=getappdata(0,'fmin');
fmax=getappdata(0,'fmax');

yu=get(handles.edit_ylabel_input,'String');

fig_num=getappdata(0,'fig_num');

leg_b=sprintf('input %6.3g %sRMS',irms,yu);
leg_a=sprintf('spec  %6.3g %sRMS',srms,yu);

t_string='Power Spectral Density';

nps=get(handles.listbox_export_plot,'Value');    

if(nps==1)
    pname=get(handles.edit_plot_name,'String');
else
    pname=' ';
end    



[fig_num]=plot_PSD_two_f(fig_num,x_label,y_label,t_string,qqq,ppp,leg_a,leg_b,fmin,fmax,pname,nps);
    


% --- Executes on selection change in listbox_export_plot.
function listbox_export_plot_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_export_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_export_plot contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_export_plot

n=get(handles.listbox_export_plot,'Value');

if(n==1)
    set(handles.text_plot_name,'Visible','on');
    set(handles.edit_plot_name,'Visible','on');
else
    set(handles.text_plot_name,'Visible','off');
    set(handles.edit_plot_name,'Visible','off');   
end


% --- Executes during object creation, after setting all properties.
function listbox_export_plot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_export_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_plot_name_Callback(hObject, eventdata, handles)
% hObject    handle to edit_plot_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_plot_name as text
%        str2double(get(hObject,'String')) returns contents of edit_plot_name as a double


% --- Executes during object creation, after setting all properties.
function edit_plot_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_plot_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_crms.
function listbox_crms_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_crms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_crms contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_crms


% --- Executes during object creation, after setting all properties.
function listbox_crms_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_crms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
