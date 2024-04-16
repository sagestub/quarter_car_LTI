function varargout = nastran_applied_force_frf(varargin)
% NASTRAN_APPLIED_FORCE_FRF MATLAB code for nastran_applied_force_frf.fig
%      NASTRAN_APPLIED_FORCE_FRF, by itself, creates a new NASTRAN_APPLIED_FORCE_FRF or raises the existing
%      singleton*.
%
%      H = NASTRAN_APPLIED_FORCE_FRF returns the handle to a new NASTRAN_APPLIED_FORCE_FRF or the handle to
%      the existing singleton*.
%
%      NASTRAN_APPLIED_FORCE_FRF('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NASTRAN_APPLIED_FORCE_FRF.M with the given input arguments.
%
%      NASTRAN_APPLIED_FORCE_FRF('Property','Value',...) creates a new NASTRAN_APPLIED_FORCE_FRF or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before nastran_applied_force_frf_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to nastran_applied_force_frf_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help nastran_applied_force_frf

% Last Modified by GUIDE v2.5 11-Feb-2019 17:19:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @nastran_applied_force_frf_OpeningFcn, ...
                   'gui_OutputFcn',  @nastran_applied_force_frf_OutputFcn, ...
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


% --- Executes just before nastran_applied_force_frf is made visible.
function nastran_applied_force_frf_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to nastran_applied_force_frf (see VARARGIN)

% Choose default command line output for nastran_applied_force_frf
handles.output = hObject;


set(handles.listbox_units,'Value',2);


num_elem_quad4stress =0; 
num_elem_tria3stress=0;

num_elem_quad4strain =0; 
num_elem_tria3strain=0;


num_node_disp=0;
num_node_velox=0;
num_node_accel=0;

setappdata(0,'num_node_disp',num_node_disp);
setappdata(0,'num_node_velox',num_node_velox);
setappdata(0,'num_node_accel',num_node_accel);

elem_quad4stress=0;
elem_tria3stress=0;

elem_quad4strain=0;
elem_tria3strain=0;


setappdata(0,'num_elem_quad4stress',num_elem_quad4stress );
setappdata(0,'num_elem_tria3stress',num_elem_tria3stress);

setappdata(0,'elem_quad4stress',elem_quad4stress);
setappdata(0,'elem_tria3stress',elem_tria3stress);

setappdata(0,'num_elem_quad4strain',num_elem_quad4strain);
setappdata(0,'num_elem_tria3strain',num_elem_tria3strain);

setappdata(0,'elem_quad4strain',elem_quad4strain);
setappdata(0,'elem_tria3strain',elem_tria3strain);




% Update handles structure
guidata(hObject, handles);

% UIWAIT makes nastran_applied_force_frf wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = nastran_applied_force_frf_OutputFcn(hObject, eventdata, handles) 
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

delete(nastran_applied_force_frf);


% --- Executes on button press in pushbutton_read.
function pushbutton_read_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_read (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

   
   fig_num=1;
   setappdata(0,'fig_num',fig_num);
 
   
   
   iu=get(handles.listbox_units,'Value');
   setappdata(0,'iu',iu);
   
   setappdata(0,'metric_displacement',0);
   setappdata(0,'metric_velocity',0);
   setappdata(0,'metric_acceleration',0);
   
   setappdata(0,'metric_quad4stress',0);
   setappdata(0,'metric_tria3stress',0);
   
   setappdata(0,'metric_quad4strain',0);
   setappdata(0,'metric_tria3strain',0);   
      
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   disp('  ');
   disp(' * * * * * * * * * * * * * * ');
   disp('  ');
 
   [filename, pathname] = uigetfile('*.f06');
   filename = fullfile(pathname, filename);
 
   try
      fid = fopen(filename,'r');
   catch
      warndlg('File not opened'); 
      return; 
   end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
   
   disp(' Read f06 file... ');
   disp('  ');   
 
   try
       sarray = textscan(fid,'%s','Delimiter','\n');
       fclose(fid);
       setappdata(0,'sarray',sarray);
    
   catch
       disp(' ');
       disp(' sarray error');
       return;
   end
        
    kv=cellfun(@length,sarray);
 
    out1=sprintf(' %d number of lines',kv);
    disp(out1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   rb_displacement = get(handles.radiobutton_displacement, 'Value');
   rb_velocity = get(handles.radiobutton_velocity, 'Value');
   rb_acceleration = get(handles.radiobutton_acceleration, 'Value');
   
   rb_plate_quad4_stress = get(handles.radiobutton_plate_quad4_stress, 'Value');
   rb_plate_quad4_strain = get(handles.radiobutton_plate_quad4_strain, 'Value');   
   
   rb_plate_tria3_stress = get(handles.radiobutton_plate_tria3_stress, 'Value');
   rb_plate_tria3_strain = get(handles.radiobutton_plate_tria3_strain, 'Value');
   
   
   rb_beam_stress = get(handles.radiobutton_beam_stress, 'Value');   
   
   setappdata(0,'rb_displacement',rb_displacement);
   setappdata(0,'rb_velocity',rb_velocity);
   setappdata(0,'rb_acceleration',rb_acceleration);
   
   setappdata(0,'rb_plate_quad4_stress',rb_plate_quad4_stress);
   setappdata(0,'rb_plate_quad4_strain',rb_plate_quad4_strain);   
   
   setappdata(0,'rb_plate_tria3_stress',rb_plate_tria3_stress);   
   setappdata(0,'rb_plate_tria3_strain',rb_plate_tria3_strain);     
   
   
   setappdata(0,'rb_beam_stress',rb_beam_stress);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

find_number_frequency_points(hObject, eventdata, handles);




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


% --- Executes on button press in radiobutton_displacement.
function radiobutton_displacement_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_displacement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_displacement


% --- Executes on button press in radiobutton_velocity.
function radiobutton_velocity_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_velocity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_velocity


% --- Executes on button press in radiobutton_acceleration.
function radiobutton_acceleration_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_acceleration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_acceleration


% --- Executes on button press in radiobutton_plate_quad4_stress.
function radiobutton_plate_quad4_stress_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_plate_quad4_stress (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_plate_quad4_stress


% --- Executes on button press in radiobutton_plate_quad4_strain.
function radiobutton_plate_quad4_strain_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_plate_quad4_strain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_plate_quad4_strain


% --- Executes on button press in radiobutton_plate_tria3_stress.
function radiobutton_plate_tria3_stress_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_plate_tria3_stress (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_plate_tria3_stress


% --- Executes on button press in radiobutton_plate_tria3_strain.
function radiobutton_plate_tria3_strain_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_plate_tria3_strain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_plate_tria3_strain


% --- Executes on button press in radiobutton_beam_stress.
function radiobutton_beam_stress_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_beam_stress (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_beam_stress

function find_number_frequency_points(hObject, eventdata, handles)
       
disp('  ');
disp(' Find number of frequency points... ');
disp('  ');     
    
sarray=getappdata(0,'sarray'); 

rb_displacement=getappdata(0,'rb_displacement');
rb_velocity=getappdata(0,'rb_velocity');
rb_acceleration=getappdata(0,'rb_acceleration');

rb_plate_quad4_stress=getappdata(0,'rb_plate_quad4_stress');
rb_plate_quad4_strain=getappdata(0,'rb_plate_quad4_strain');


rb_plate_tria3_stress=getappdata(0,'rb_plate_tria3_stress');
rb_plate_tria3_strain=getappdata(0,'rb_plate_tria3_strain');


rb_beam_stress=getappdata(0,'rb_beam_stress');


num_node_disp=0;
num_node_velox=0;
num_node_accel=0;

num_elem_quad4stress =0;
num_elem_quad4strain =0;

num_elem_tria3stress=0;   
num_elem_tria3strain=0; 

num_elem_beamstress=0;
   
    if(rb_displacement==1)    
        idd = strfind(sarray{:}, 'C O M P L E X   D I S P L A C E M E N T   V E C T O R');
        idd = find(not(cellfun('isempty', idd)));
        num_node_disp = length(idd);
    end
    
    if(rb_velocity==1)
        idv = strfind(sarray{:}, 'C O M P L E X   V E L O C I T Y   V E C T O R');
        idv = find(not(cellfun('isempty', idv)));
        num_node_velox = length(idv);
    end
    
    if(rb_acceleration==1)
        ida = strfind(sarray{:}, 'C O M P L E X   A C C E L E R A T I O N   V E C T O R');
        ida = find(not(cellfun('isempty', ida)));
        num_node_accel = length(ida);
        setappdata(0,'num_node_accel',num_node_accel);
        
        if(num_node_accel==0)
            warndlg('Acceleration required');
            return;
        end
        
    end
    
    if(rb_plate_quad4_stress==1)
        iquad4stress = strfind(sarray{:}, 'S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )');
        iquad4stress = find(not(cellfun('isempty', iquad4stress)));
        num_elem_quad4stress = length(iquad4stress);   
    end
     if(rb_plate_quad4_strain==1)
        iquad4strain = strfind(sarray{:}, 'S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )');
        iquad4strain = find(not(cellfun('isempty', iquad4strain)));
        num_elem_quad4strain = length(iquad4strain);   
    end   
    
    
    
    if(rb_plate_tria3_stress==1)
        itria3stress = strfind(sarray{:},'S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )');
        itria3stress = find(not(cellfun('isempty', itria3stress)));
        num_elem_tria3stress = length(itria3stress);   
    end
    if(rb_plate_tria3_strain==1)
        itria3strain = strfind(sarray{:},'S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )');
        itria3strain = find(not(cellfun('isempty', itria3strain)));
        num_elem_tria3strain = length(itria3strain);   
    end
    

    if(rb_beam_stress==1)    
        ibeam_stress = strfind(sarray{:}, 'C O M P L E X   S T R E S S E S   I N   B E A M   E L E M E N T S');
        ibeam_stress = find(not(cellfun('isempty', ibeam_stress)));
        num_elem_beamstress = length(ibeam_stress);
    end    
    
    
    
    
    idf = strfind(sarray{:}, 'FREQUENCY =');
    idf = find(not(cellfun('isempty', idf)));
%%%    num_freq = length(idf); 

    
    num_freq=1;

    tb=0;
    
    i=1;
    
    ijk=1;
    
    while(1)
    
        j=idf(i);   
        ss=sarray{1}{j};        
        strs = strsplit(ss,' ');

        if( strcmp(char(strs(1)),'FREQUENCY'))

            tb=str2double(strs(3));
            fmin=tb; 
            fn(ijk)=tb;
            ione=i+1;
            break;
        end    
            
        i=i+1;    
        
    end
    
    for i=(ione-2):-1:1
        idf(i)=[];
    end
    
    ijk=ijk+1;
    
    
    idfL=length(idf);
    
    out1=sprintf(' idfL = %d ',idfL);
    disp(out1);
    
    progressbar;
    
    for i=2:idfL
        
        progressbar(i/idfL);
        
        j=idf(i);
        
        ss=sarray{1}{j};
        strs = strsplit(ss,' ');
        tt=str2double(strs(3));
        
        if(tt>tb)
            tb=tt;
            fn(ijk)=tb;
            ijk=ijk+1;
        end
    end

    
    pause(0.2);
    progressbar(1);
        
    metric_displacement=0;
    metric_velocity=0;
    metric_acceleration=0;
    
    metric_quad4stress=0;
    metric_quad4strain=0;    
    
    metric_tria3stress=0;
    metric_tria3strain=0;    
    
    metric_beamstress=0;
    
    if(num_node_disp>0)
        metric_displacement=1;
    end
    if(num_node_velox>0)
        metric_velocity=1;          
    end   
    if(num_node_accel>0)
        metric_acceleration=1;          
    end
    
    if(num_elem_quad4stress >0)
        metric_quad4stress=1;          
    end 
    if(num_elem_quad4strain >0)
        metric_quad4strain=1;          
    end     
    
    if(num_elem_tria3stress>0)
        metric_tria3stress=1;          
    end  
    if(num_elem_tria3strain>0)
        metric_tria3strain=1;          
    end      
    
    if(num_elem_beamstress>0)
        metric_beamstress=1;          
    end     
    
    
    if(rb_displacement==0 || num_node_disp==0)
        metric_displacement=0;      
    end 
    if(rb_velocity==0 || num_node_velox==0)
        metric_velocity=0;           
    end 
    if(rb_acceleration==0 || num_node_accel==0)
        metric_acceleration=0;       
    end 
    
    if(rb_plate_quad4_stress==0 || num_elem_quad4stress ==0)
        metric_quad4stress=0;      
    end     
    if(rb_plate_quad4_strain==0 || num_elem_quad4strain ==0)
        metric_quad4strain=0;      
    end         
    
    
    if(rb_plate_tria3_stress==0 || num_elem_tria3stress==0)
        metric_tria3stress=0;      
    end  
    if(rb_plate_tria3_strain==0 || num_elem_tria3strain==0)
        metric_tria3strain=0;      
    end  
    
    
    if(rb_beam_stress==0 || num_elem_beamstress==0)
        metric_beamstress=0;      
    end     
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
   out1=sprintf(' num_node_disp=%d  num_node_velox=%d  num_node_accel=%d ',num_node_disp,num_node_velox,num_node_accel);
   disp(out1);
   
   out1=sprintf(' num_elem_quad4stress =%d num_elem_quad4strain=%d ',num_elem_quad4stress,num_elem_quad4strain);
   disp(out1);
   
   out1=sprintf(' num_elem_tria3stress =%d num_elem_tria3strain=%d ',num_elem_tria3stress,num_elem_tria3strain);
   disp(out1);
   
   
   out1=sprintf(' num_elem_beamstress =%d ',num_elem_beamstress);
   disp(out1)
   
   zflag=0;
   setappdata(0,'zflag',zflag);
   
   nf=length(fn);
   
   if(nf==0)
       out1=sprintf('f06 file must contain displacement, velocity, acceleration, \n quad4 stress, tria3 stress or beam stress');
       warndlg(out1);
       zflag=1;
       setappdata(0,'zflag',zflag);
       return;
   end
   
    
   disp(' ');
   disp(' Form frequency vector...');
   disp(' ');   
   
   fn=fix_size(fn);
      
   setappdata(0,'fn',fn);
   
   out1=sprintf(' fmin=%9.5g  fmax=%9.5g  nf=%9.5g',fn(1),fn(nf),nf);
   disp(out1);
   
   if(rb_displacement==1)
        setappdata(0,'idd',idd);   
   end 
   if(rb_velocity==1)
        setappdata(0,'idv',idv);          
   end 
   if(rb_acceleration==1)
        setappdata(0,'ida',ida);
   end 
   if(rb_plate_quad4_stress==1)
        setappdata(0,'iquad4stress',iquad4stress);
   end     
   if(rb_plate_quad4_strain==1)
        setappdata(0,'iquad4strain',iquad4strain);
   end    
   
   if(rb_plate_tria3_stress==1)
        setappdata(0,'itria3stress',itria3stress);      
   end
   if(rb_plate_tria3_strain==1)
        setappdata(0,'itria3strain',itria3strain);      
   end   
   
   
   if(rb_beam_stress==1)
        setappdata(0,'ibeam_stress',ibeam_stress);      
   end      
   
   
   try
        setappdata(0,'idf',idf);
   catch
        warndlg('Frequency column error');
        return;
   end
        
   
   disp(' find frequency')
   
       
   setappdata(0,'num_node_disp',num_node_disp);
   setappdata(0,'num_node_velox',num_node_velox);
   setappdata(0,'num_node_accel',num_node_accel);
   
   setappdata(0,'num_elem_quad4stress',num_elem_quad4stress );
   setappdata(0,'num_elem_quad4strain',num_elem_quad4strain );   
   
   setappdata(0,'num_elem_tria3stress',num_elem_tria3stress);
   setappdata(0,'num_elem_tria3strain',num_elem_tria3strain);  
   
   setappdata(0,'num_elem_beamstress',num_elem_beamstress);  
   setappdata(0,'num_freq',num_freq);

   setappdata(0,'metric_displacement',metric_displacement);
   setappdata(0,'metric_velocity',metric_velocity);
   setappdata(0,'metric_acceleration',metric_acceleration);
   
   setappdata(0,'metric_quad4stress',metric_quad4stress);
   setappdata(0,'metric_quad4strain',metric_quad4strain);   
   
   setappdata(0,'metric_tria3stress',metric_tria3stress);
   setappdata(0,'metric_tria3strain',metric_tria3strain);   
   
   setappdata(0,'metric_beamstress',metric_beamstress);
   
   setappdata(0,'nf',nf);   
   setappdata(0,'fn',fn);  
  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

   if(metric_acceleration==1 && rb_acceleration==1)

        find_acceleration_nodes(hObject, eventdata, handles);    
   
        nflag=getappdata(0,'nflag');
        
        if(nflag==1)
            return;
        end
        
   end    

%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
   if(metric_displacement==1 && rb_displacement==1)
       
        find_displacement_nodes(hObject, eventdata, handles);    
        
        qflag=getappdata(0,'qflag');
        
        if(qflag==1)
            return;
        end
        
   end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   if(metric_velocity==1 && rb_velocity==1)
        find_velocity_nodes(hObject, eventdata, handles);    
   end    

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
   if(metric_quad4stress==1 && rb_plate_quad4_stress==1)
        find_quad4stress_elem(hObject, eventdata, handles);          
   end  
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
   if(metric_quad4strain==1 && rb_plate_quad4_strain==1)
        find_quad4strain_elem(hObject, eventdata, handles);          
   end   
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   if(metric_tria3stress==1 && rb_plate_tria3_stress==1)
        find_tria3stress_elem(hObject, eventdata, handles);  
   end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   if(metric_tria3strain==1 && rb_plate_tria3_strain==1)
        find_tria3strain_elem(hObject, eventdata, handles);  
   end   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

   if(metric_beamstress==1 && rb_beam_stress==1)
        find_beamstress_elem(hObject, eventdata, handles);      
   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  


function find_velocity_nodes(hObject, eventdata, handles)

    sarray=getappdata(0,'sarray');

    idv=getappdata(0,'idv');  
    
    nf=getappdata(0,'nf');
    fn=getappdata(0,'fn');
    
    disp('    ');
    disp(' Find velocity response nodes ');
    disp('    ');
        
    j=1;
    
    jk=1;
        
    iflag=0;
        
    for i=idv(1):idv(2)
                        
        ss=sarray{1}{i};
            
        k=strfind(ss, 'POINT ID.');
            
        if(~isempty(k))
                
                while(1) 
                    ss=sarray{1}{i+j};
                  
                    kk=strfind(ss,'FRF');
                    
                    if(~isempty(kk))
                        iflag=1;
                        break;
                    else
                        strs = strsplit(ss,' ');
                        if( strcmp(char(strs(3)),'G'))
                            node_velox(jk)= str2double(strs(2));
                            jk=jk+1;
                        end
                    end
                   
                    j=j+1;
                    
                end
                
                if(iflag==1)
                    break;
                end                    
        end
    end
   
   num_node_velox=length(node_velox);
   
   for i=1:num_node_velox
       out1=sprintf('  %d',node_velox(i));
       disp(out1);
   end    
   
   setappdata(0,'num_node_velox',num_node_velox);
   setappdata(0,'node_velox',node_velox);
    
   
   qflag=0;
   setappdata(0,'qflag',qflag);
   
disp('    ');
disp(' Form velocity FRF ');
disp('    ');
        
nm=max(node_velox);
        
node_index=zeros(nm,1);
        
for i=1:num_node_velox
    j=node_velox(i);
    node_index(j)=i;
end
                
progressbar;
            
nndd=nf;
        
TM1=zeros(nf,num_node_velox);
TM2=zeros(nf,num_node_velox);
TM3=zeros(nf,num_node_velox);
RM1=zeros(nf,num_node_velox);
RM2=zeros(nf,num_node_velox);
RM3=zeros(nf,num_node_velox);        
                 
TP1=zeros(nf,num_node_velox);
TP2=zeros(nf,num_node_velox);
TP3=zeros(nf,num_node_velox);
RP1=zeros(nf,num_node_velox);
RP2=zeros(nf,num_node_velox);
RP3=zeros(nf,num_node_velox);   
 
TC1=zeros(nf,num_node_velox);
TC2=zeros(nf,num_node_velox);
TC3=zeros(nf,num_node_velox);
RC1=zeros(nf,num_node_velox);
RC2=zeros(nf,num_node_velox);
RC3=zeros(nf,num_node_velox);  
        
 
 
C=pi/180;
 
for j=1:nndd
                
    progressbar(j/nndd); 
            
    k1=idv(j);
                
    k=k1-1;
                
    while(1)
                    
        k=k+1;
                                        
        sst=sarray{1}{k};             
                    
        kt=strfind(sst, 'NASTRAN');
                    
        if(~isempty(kt))
            break;
        end
                    
        ka=strfind(sst, 'POINT ID.');
            
        if(~isempty(ka))
                        
            for ijk=1:2:2*num_node_velox 
                sst=sarray{1}{k+ijk};
                strs = strsplit(sst,' ');
                iv=node_index(str2num(char(strs(2))));
                        
                TM1(j,iv)=str2double(char(strs(4)));
                TM2(j,iv)=str2double(char(strs(5)));
                TM3(j,iv)=str2double(char(strs(6)));  
                RM1(j,iv)=str2double(char(strs(7)));
                RM2(j,iv)=str2double(char(strs(8)));
                RM3(j,iv)=str2double(char(strs(9)));   
                
                sst=sarray{1}{k+ijk+1};
                strs = strsplit(sst,' ');
                
                TP1(j,iv)=C*str2double(char(strs(1)));
                TP2(j,iv)=C*str2double(char(strs(2)));
                TP3(j,iv)=C*str2double(char(strs(3)));  
                RP1(j,iv)=C*str2double(char(strs(4)));
                RP2(j,iv)=C*str2double(char(strs(5)));
                RP3(j,iv)=C*str2double(char(strs(6)));   
         
                TC1(j,iv)=TM1(j,iv)*( cos(TP1(j,iv)) + 1i*sin(TP1(j,iv))  );
                TC2(j,iv)=TM2(j,iv)*( cos(TP2(j,iv)) + 1i*sin(TP2(j,iv))  );
                TC3(j,iv)=TM3(j,iv)*( cos(TP3(j,iv)) + 1i*sin(TP3(j,iv))  );
                RC1(j,iv)=RM1(j,iv)*( cos(RP1(j,iv)) + 1i*sin(RP1(j,iv))  );
                RC2(j,iv)=RM2(j,iv)*( cos(RP2(j,iv)) + 1i*sin(RP2(j,iv))  );
                RC3(j,iv)=RM3(j,iv)*( cos(RP3(j,iv)) + 1i*sin(RP3(j,iv))  );                  
                
                
            end
                        
            break;
                           
        end
    end
end
 
%%
pause(0.2);
progressbar(1);       
        
disp('  ');
disp(' Writing velocity arrays'); 
disp('  ');        

progressbar;
                 
for i=1:num_node_velox
 
    node=node_velox(i);
        
    output_TM1=sprintf('velox_%d_mag_frf_T1',node);
    output_TM2=sprintf('velox_%d_mag_frf_T2',node);            
    output_TM3=sprintf('velox_%d_mag_frf_T3',node);         
    output_RM1=sprintf('velox_%d_mag_frf_R1',node);
    output_RM2=sprintf('velox_%d_mag_frf_R2',node);            
    output_RM3=sprintf('velox_%d_mag_frf_R3',node); 
        
    assignin('base', output_TM1, [fn TM1(:,i)]);            
    assignin('base', output_TM2, [fn TM2(:,i)]); 
    assignin('base', output_TM3, [fn TM3(:,i)]);             
    assignin('base', output_RM1, [fn RM1(:,i)]);            
    assignin('base', output_RM2, [fn RM2(:,i)]); 
    assignin('base', output_RM3, [fn RM3(:,i)]);
    
end
 
progressbar;
