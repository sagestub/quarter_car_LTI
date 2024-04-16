%
%   quake_srs_synth_phase_1.m  ver 1.7  by Tom Irvine
%
function[ramp,error_i,wavelet_low,wavelet_up,alpha,beta,t]=...
    quake_srs_synth_phase_1(f,dt,damp,amp,NHS,td,last_wavelet,iter,...
                            fn,srs_spec,dlimit,error_limit,iu,uewf,dur)
                                         
    
tpi=2*pi;                        
                        
%%

% Need to do this here.  Input data may have different dt than output

nt=round(dur/dt);

t=linspace(0,nt*dt,nt+1); 


alpha=zeros(last_wavelet,1);
wavelet_low=zeros(last_wavelet,1);
wavelet_up =zeros(last_wavelet,1);

 
beta=tpi*f;
 
for i=1:last_wavelet
    alpha(i)=beta(i)/double(NHS(i));
    ud=NHS(i)/(2*f(i)); 
    
    wavelet_low(i)=floor(td(i)/dt);
    
    if(wavelet_low(i)<=0)
        wavelet_low(i)=1;       
    end    
    
    wavelet_up(i)=wavelet_low(i)+floor(ud/dt);   
%    

    if(wavelet_up(i)>nt)
        wavelet_up(i)=nt;       
    end   
    
end


%%
                                         
                                         
ic=0;
                                     
ramp=amp;

nfn=length(fn);

ih=round(iter/3);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Initialize coefficients
%
[a1,a2,b1,b2,b3,rd_a1,rd_a2,rd_b1,rd_b2,rd_b3]=...
                                   srs_coefficients(fn,damp,dt);    

drec=1.0e+90;    

error_max=1.0e+90;                               
         

for i=1:length(f)
    if(f(i)>max(fn))
        amp(i)=0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

progressbar;
%
for i=1:iter
    progressbar(i/iter);
%
    [accel]=generate_time_history_wavelet_table_q(f,amp,NHS,td,t,...
                                         wavelet_low,wavelet_up,alpha,beta);
    
    xabs=zeros(nfn,1);
    xmax=zeros(nfn,1);
    xmin=zeros(nfn,1);
        
    
    for j=1:nfn
       
        forward=[ b1(j),  b2(j),  b3(j) ];    
        back   =[     1, -a1(j), -a2(j) ];    
%    
        resp=filter(forward,back,accel);
%
        xmax(j)=max(resp);
        xmin(j)=abs(min(resp));
        
        xabs(j)= max([xmax(j) xmin(j)]);
    end
 

    
    error=0;
    error_i=0;
    
    e2=0;
    
    fmin=0.;
    
    for j=1:length(srs_spec)
        if(srs_spec(j)>0)
            errora=abs(20*log10(srs_spec(j)/xmin(j)));
            errorb=abs(20*log10(srs_spec(j)/xmax(j)));
            
            if(xmin(j)>srs_spec(j))
                errora=uewf*errora;
            end
            if(xmax(j)>srs_spec(j))
                errorb=uewf*errorb;
            end           
            
            
            if(errora>error)
                error=errora;
                fmax=fn(j);
            end
            if(errorb>error)
                error=errorb;
                fmin=fn(j);
            end
            
            error_i=error;
            
            e2=e2+errora+errorb;
        end    
    end 
    
    error=error*e2;

  
    
    [v]=integrate_function(accel,dt);
    [d]=integrate_function(v,dt);

    
    dmin=abs(min(d));
    dmax=abs(max(d));
    
      
    dratio=max([dmax dmin]);
    
    if(iu==1)
        dratio=dratio*386;
    else
        dratio=dratio*9.81*1000;
    end   
    
    
    if(error<error_max && (dratio<dlimit || dratio <drec));
        
        if(dratio<drec)
            drec=dratio;
        end
        
        error_max=error;
        ramp=amp;
        raccel=accel;
        rxmax=xmax;
        rxmin=xmin;
        

        
        out1=sprintf('%d  %8.4g  %8.4g  %8.4g %d ',i,error,error_i,dratio,ic);
        disp(out1);
        
        if(error_i<error_limit)
            progressbar(1);
            break;
        end
    end
    
%
    
    if(i<iter)
        if(i> round(iter*0.03) || i>100)

            if(rand()<0.5)
                if(rand()<0.6 || i > ih)
                    amp=ramp;
                    jk=1+floor(last_wavelet*rand());
                    amp(jk)=ramp(jk)*(0.95+0.1*rand());
                    ic=1;
                else    
                    if(rand()<0.7)
                        for j=1:last_wavelet
                            amp(j)=ramp(j)*(0.99+0.02*rand());
                        end
                        ic=2;
                    else
                        for j=1:last_wavelet
                            amp(j)=ramp(j)*(0.95+0.10*rand());
                        end  
                        ic=3;
                    end 
                end
            else
                if(rand()<0.5)
                    ferr=1.0e+20;
                    for jj=1:length(f)
                        diff=abs(f(jj)-fmax);
                        if(diff<ferr)
                            jk=jj;
                            ferr=diff;
                        end
                    end
                    amp(jk)=ramp(jk)*(0.95+0.1*rand());
                    ic=4;
                else
                    ferr=1.0e+20;
                    for jj=1:length(f)
                        diff=abs(f(jj)-fmin);
                        if(diff<ferr)
                            jk=jj;
                            ferr=diff;
                        end
                    end
                    amp(jk)=ramp(jk)*(0.95+0.1*rand());
                    ic=5;                   
                end
            end
        else
            scale=0.5+1.0*rand();
            amp=ramp*scale;
        end    
    end
%
end

pause(0.5);