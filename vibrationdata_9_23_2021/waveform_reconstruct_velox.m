%
 clear THM; 
 clear aaa;
 clear amp;
 clear amp_original;
 clear ddd;
 clear delay;
 clear disp;
 clear duration;
 clear error;
 clear f;
 clear fl;
 clear fmax;
 clear freq;
 clear fu;
 clear iunit;
 clear mm;
 clear nfr;
 clear nhs;
 clear nt;
 clear num2;
 clear out1;
 clear p;
 clear residual;
 clear sd;
 clear sr;
 clear t;
 clear tp;
 clear tpi; 
 clear velox;
 clear vlast;
 clear vvv;
 clear x1r;
 clear x2r;
 clear x3r;
 clear x4r;
 clear y;
 clear ffmax;
 clear ffmin;
%
disp(' ');
disp(' wavelet_reconstruct.m   version 1.1   June 27, 2009');
disp(' By Tom Irvine   Email:  tomirvine@aol.com ');
%
disp(' ');
disp(' This program synthesizes a wavelet series to approximate ');
disp(' an input time history. ');
%
disp(' ');
disp(' The input file must be: time(sec) and amplitude(units) ');
%
MAX=50000;
MAXP=1000;
%
tp=2.*pi;
tpi=tp;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
disp(' ');
iunit=input(' Enter acceleration unit:   1= G   2= m/sec^2  ');	
%
if(iunit==1)
    iscale=386.;
else
    iscale=1.;
end
%
disp(' ')
disp(' Select file input method ');
disp('   1=external ASCII file ');
disp('   2=file preloaded into Matlab ');
disp('   3=Excel file');
file_choice = input('');
%
if(file_choice==1)
        [filename, pathname] = uigetfile('*.*');
        filename = fullfile(pathname, filename);
        fid = fopen(filename,'r');
        THM = fscanf(fid,'%g %g',[2 inf]);
        fclose(fid);
        THM=THM';
end
if(file_choice==2)
        THM = input(' Enter the matrix name:  ');
end
if(file_choice==3)
        [filename, pathname] = uigetfile('*.*');
        xfile = fullfile(pathname, filename);
%        
        THM = xlsread(xfile);
%         
end
%
t=double(THM(:,1));
y=double(THM(:,2));
amp=y;
amp_original=y;
%
disp(' ')
ffmin=input(' Enter minimum frequency (Hz) ');
disp(' ')
ffmax=input(' Enter maximum frequency (Hz) ');
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
ik=0;
%      
x1r=zeros(MAXP,1);	
x2r=zeros(MAXP,1);	
x3r=zeros(MAXP,1);	
x4r=zeros(MAXP,1);	
%
disp(' ');		
nt=input(' Enter the number of trials per frequency. ');
disp(' ');		
nfr=input(' Enter the number of frequencies. ');
%
if(nfr>MAXP)
    nfr=MAXP;
end
%		   
fmax=0.;
%
residual=y;
%
num2=max(size(y));
%
clear y;
%
out1=sprintf(' number of input point= %d ',num2)
disp(out1);
%		
duration=t(num2)-t(1);
%
sr=duration/(num2-1);
sr=1./sr;
%
out1=sprintf(' sample rate = %10.4g \n',sr);
disp(out1);
%
fl=3./duration;
fu=sr/5.;
%
clear y;
%
for(ie=1:nfr)
%
    out1=sprintf(' frequency case %d ',ie);
    disp(out1);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
    disp(' ref 1');
    [x1r,x2r,x3r,x4r]=wgen(num2,t,residual,duration,sr,x1r,x2r,x3r,x4r,fl,fu,nt,ie,ffmin,ffmax);
    disp(' ref 2');
%            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
    if(x2r(ie)<ffmin*tp)
        x2r(ie)=ffmin*tp;
        x1r(ie)=1.0e-20;
        x3r(ie)=3;
        x4r(ie)=0;
    end
%
	for(i=1:num2)
%			
	    tt=t(i);
%
		t1=x4r(ie) + t(1);
		t2=t1 + tp*x3r(ie)/(2.*x2r(ie));
%
		y=0.;
%
		if( tt>= t1 && tt <= t2)
%				
			arg=x2r(ie)*(tt-t1);  
			y=x1r(ie)*sin(arg/double(x3r(ie)))*sin(arg);			       
%
		    residual(i)=residual(i)-y;
        end
    end       
    ave=mean(residual); 
    sd=std(residual);
%
    out1=sprintf(' ave=%12.4g  sd=%12.4g \n\n',ave,sd);
    disp(out1);       
%			 
end    
%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
clear aaa;
clear vvv;
clear ddd;
clear mm;
clear w_diff;
%
aaa=zeros(num2,1); 
vvv=zeros(num2,1); 
ddd=zeros(num2,1); 
w_diff=zeros(num2,1); 
%
		disp(' ');
%
		kk=0;
%
		for(i=1:nfr)
%
            if(x2r(i)<ffmin*tp)
                x2r(i)=ffmin*tp;
                x1r(i)=1.0e-20;
                x3r(i)=3;
                x4r(i)=0;
            end
%
			out1=sprintf(' amp=%10.4f   freq=%10.3f Hz   nhs=%d   delay=%10.4f ',x1r(i),x2r(i)/tp,x3r(i),x4r(i));
            disp(out1);
%
            wavelet_table(i,1)=i;
            wavelet_table(i,2)=x2r(i)/tpi;           
            wavelet_table(i,3)=x1r(i);
            wavelet_table(i,4)=x3r(i); 
            wavelet_table(i,5)=x4r(i);              
%
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		tf=0.;
%
		vlast=zeros(num2,1);
%
%%        num2
%%        nfr
%%        iscale
		for(k=1:num2)
%
				tt=t(k);
%              
				for(j=1:nfr)
%				
                    w=0.;
					v=0.;
					d=0.;
%
					v1=0.;
					v2=0.;
%
					d1=0.;
					d2=0.;
%
					t1=x4r(j)+t(1);
%
                    if(x2r(j)>=1.0e-10)
%
                    else
                        x2r(j)=ffmin*tp;
                        x1r(j)=1.0e-20;
                        x3r(j)=3;
                        x4r(j)=0;
                        t1=x4r(j)+t(1);  
                    end    
                    t2=tp*x3r(j)/(2.*x2r(j))+t1; 
%
					if( tt>=t1  && tt <=t2  )
%					
						arg=x2r(j)*(tt-t1);  
%
						w=  x1r(j)*sin(arg/double(x3r(j)))*sin(arg);
%
                           wd=cos(arg/double(x3r(j)))*sin(arg)/double(x3r(j)); 
                        wd=wd+sin(arg/double(x3r(j)))*cos(arg);
                        wd=wd*x1r(j)*x2r(j);
%
                        aa=x2r(j)/double(x3r(j));
						bb=x2r(j);
%
						alpha1=aa+bb;
					    alpha2=aa-bb;
%
						te=tt-t1;
%
						v1= -sin(alpha1*tf)/(2.*alpha1) + sin(alpha2*tf)/(2.*alpha2);
				        v2= -sin(alpha1*te)/(2.*alpha1) + sin(alpha2*te)/(2.*alpha2);
%
						v=(v2-v1)*iscale*x1r(j);
%
						vlast(j)=v;
%
                    end
%
                    w_diff(k)=w_diff(k)+wd;
					vvv(k)=vvv(k)+w; 
					ddd(k)=ddd(k)+v;
%
					amp(k)=amp(k)-w;
%
					if(x3r(j)<1)
                        printf(' error x3r ');
                        break;
                    end
%
                end
%
%%               out1=sprintf(' %d t=%8.4g  %8.4g  %8.4g  %8.4g ',k,t(k),aaa(k),vvv(k),ddd(k));
%%               disp(out1);
        end
%%
        aaa=w_diff;
        error=[t,amp];
        accel=[t,aaa];
        velox=[t,vvv];
        sz=size(ddd);
        if(sz(2)>sz(1))
            ddd=ddd';
        end
        disp=[t,ddd];
%
        figure(1);
        plot(t,ddd);     
        title(' Displacement Time History ');
        xlabel(' Time(sec)');
        if(iunit==1)
            ylabel(' Displacment(inch)');
        else
            ylabel(' Displacement(m)');          
        end
        grid on;
 %
        figure(2);
        plot(t,amp_original,t,vvv);
        legend ('raw data','synthesis');       
        title(' Velocity Time History ');
        xlabel(' Time(sec)');
        if(iunit==1)
            ylabel(' Velocity(in/sec)');
        else
            ylabel(' Velocity(m/sec)');          
        end
        grid on;
%        
        figure(3);
        plot(t,aaa);     
        title(' Acceleration Time History ');
        xlabel(' Time(sec)');
        if(iunit==1)
            ylabel(' Accel(G)');
        else
            ylabel(' Accel(m/sec^2)');          
        end
        grid on;
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
 clear THM; 