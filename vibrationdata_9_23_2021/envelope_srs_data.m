function varargout = envelope_srs_data(varargin)
% ENVELOPE_SRS_DATA MATLAB code for envelope_srs_data.fig
%      ENVELOPE_SRS_DATA, by itself, creates a new ENVELOPE_SRS_DATA or raises the existing
%      singleton*.
%
%      H = ENVELOPE_SRS_DATA returns the handle to a new ENVELOPE_SRS_DATA or the handle to
%      the existing singleton*.
%
%      ENVELOPE_SRS_DATA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ENVELOPE_SRS_DATA.M with the given input arguments.
%
%      ENVELOPE_SRS_DATA('Property','Value',...) creates a new ENVELOPE_SRS_DATA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before envelope_srs_data_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to envelope_srs_data_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help envelope_srs_data

% Last Modified by GUIDE v2.5 12-Oct-2017 15:37:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @envelope_srs_data_OpeningFcn, ...
                   'gui_OutputFcn',  @envelope_srs_data_OutputFcn, ...
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


% --- Executes just before envelope_srs_data is made visible.
function envelope_srs_data_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to envelope_srs_data (see VARARGIN)

% Choose default command line output for envelope_srs_data
handles.output = hObject;

set(handles.uipanel_save,'Visible','off');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes envelope_srs_data wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = envelope_srs_data_OutputFcn(hObject, eventdata, handles) 
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

delete(envelope_srs_data);


% --- Executes on button press in pushbutton_calculate.
function pushbutton_calculate_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_calculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fig_num=1;

try
    FS=get(handles.edit_input_array,'String');
    THF=evalin('base',FS);
catch
    warndlg('Input file not found');
    return; 
end

dB=str2num(get(handles.edit_margin,'String'));

mscale=10^(dB/20);


THF=sortrows(THF,1);

ioct=6;

fr=THF(:,1);
 r=THF(:,2);

[f,srs_in]=SRS_specification_interpolation_nw(fr,r,ioct);


omegan=2*pi*f;

nf=length(f);

fmin=f(1);
fmax=f(nf);

amin=min(srs_in);
amax=max(srs_in);


ccc_record=1.0e+99;
rmax=1.0e+99;

for i=1:2000
   
    if(i<20 || rand()<0.5)
    
        ijk=round(nf*rand());

        
        if(ijk<=1)
            ijk=2;
        end
        if(ijk>(nf-1))
            ijk=nf-1;
        end
        
        fint=f(ijk);
        aa=amax*(rand())^3;
    
    else
        
        fint=frec*(0.95+0.1*rand());
          aa=aarec*(0.95+0.1*rand());
        
    end  

        
    if(aa<amin)
        aa=amin;
    end
    
    fr=[fmin fint fmax];
    r=[ aa   amax  amax];
   
    [f,srs_trial]=SRS_specification_interpolation_nw(fr,r,ioct);    
    
% scale

     sc=ones(nf,1);
     


     for j=1:nf
         
         if(srs_trial(j)<srs_in(j))
            
             sc(j)= srs_in(j)/srs_trial(j);
             
         end
         
     end
     

     max_sc=max(sc);
     
     srs_trial=srs_trial*max_sc;
             r=r*max_sc;

% calculate PV

     PV=zeros(nf,1);

     for j=1:nf
         
         PV(j)=386*srs_trial(j)/omegan(j);
         
     end
     
     PV_max=max(PV);
     
     ccc=r(1)*PV_max*max(r);
     
     if(ccc<ccc_record)
        
         ccc_record=ccc;

      
         fr=fix_size(fr);
          r=fix_size(r);
         
         r=r*mscale; 
          
         spec=[fr r];
         
         out1=sprintf('\n i=%d  PV_max=%8.4g  ',i,PV_max);
         disp(out1);
         
         spec  
         
         frec=fint;
         aarec=aa;
         
     end

% check for ccc_record
    
    
end

x_label='Natural Frequency (Hz)';
y_label='Peak Accel (G)';
t_string='Shock Response Spectrum Q=10';

ppp1=spec;
ppp2=THF;

leg1='envelope';
leg2='measure data';

md=5;

[fig_num,h2]=plot_loglog_function_md_two_h2(fig_num,x_label,...
               y_label,t_string,ppp1,ppp2,leg1,leg2,fmin,fmax,md);

spec           
           
setappdata(0,'new_spec',spec);

set(handles.uipanel_save,'Visible','on');

out1=sprintf('\n mscale = %6.3g \n',mscale);
disp(out1);


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



function edit_margin_Callback(hObject, eventdata, handles)
% hObject    handle to edit_margin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_margin as text
%        str2double(get(hObject,'String')) returns contents of edit_margin as a double


% --- Executes during object creation, after setting all properties.
function edit_margin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_margin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function pushbutton_calculate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_calculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function edit_output_array_name_Callback(hObject, eventdata, handles)
% hObject    handle to edit_output_array_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_output_array_name as text
%        str2double(get(hObject,'String')) returns contents of edit_output_array_name as a double


% --- Executes during object creation, after setting all properties.
function edit_output_array_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_output_array_name (see GCBO)
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



data=getappdata(0,'new_spec');

data

output_name=get(handles.edit_output_array_name,'String');
assignin('base', output_name, data);

h = msgbox('Save Complete');


% --- Executes on key press with focus on edit_input_array and none of its controls.
function edit_input_array_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to edit_input_array (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

set(handles.uipanel_save,'Visible','off');


% --- Executes on key press with focus on edit_margin and none of its controls.
function edit_margin_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to edit_margin (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

set(handles.uipanel_save,'Visible','off');
