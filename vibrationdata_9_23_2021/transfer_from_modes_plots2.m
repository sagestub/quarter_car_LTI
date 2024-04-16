%
%  transfer_from_modes_plots2.m  ver 1.3  by Tom Irvine
%
function[fig_num]=...
    transfer_from_modes_plots2(iu,iam,i,k,fig_num,minf,maxf,ymin,ymax,freq,PPP,PHA)
%
    ymax=10^(ceil(0.1+log10(ymax)));
%
    ymin=10^(floor(log10(ymin)));
%
    if(ymin<ymax*0.0001)
          ymin=ymax*0.0001;
    end

%
    freq=fix_size(freq);
%

    try

        PPP=fix_size(PPP); 
%
        af=[freq PPP]; 
        [xmax,max_freq]=find_max(af);
    
    end
%

    try

        figure(fig_num);
        fig_num=fig_num+1;
        subplot(3,1,1);    
        plot(freq,PHA);
        axis([minf,maxf,-180,180]); 
        set(gca,'ytick',[-180 0 180]);    
        if(iam==1)
            out1=sprintf('Receptance FRF  H %d %d',i,k);
        end  
        if(iam==2)
            out1=sprintf('Mobility FRF  H %d %d',i,k);
        end  
        if(iam==3)
            out1=sprintf('Accelerance FRF  H %d %d',i,k);
        end    
        if(iam==6)
            out1=sprintf('Bending Stress FRF  H %d %d',i,k);
        end     
        title(out1);     
        grid on;
        ylabel('Phase(deg)');
        set(gca,'MinorGridLineStyle',':','GridLineStyle',':','XScale','log');
    
        [xtt,xTT,iflag]=xtick_label(minf,maxf);

        if(iflag==1)
            set(gca,'xtick',xtt);
            set(gca,'XTickLabel',xTT);
            fmin=min(xtt);
            fmax=max(xtt);
            xlim([fmin,fmax]);    
        end   

    
    %
        subplot(3,1,[2 3]);       
        plot(freq,PPP);
        grid on;
        axis([minf,maxf,ymin,ymax]);
        set(gca,'MinorGridLineStyle',':','GridLineStyle',':','XScale','log','YScale','log');
        xlabel('Frequency (Hz)');    
    %
    %
        if(iam==1)
            disp(' ');
            disp(' * * * * Maximum Values * * * * ');  
            if(iu==1)
                ylabel('Disp/Force (in/lbf)');
                out1=sprintf('\n (Disp/Force) = %8.4g (in/lbf) at %8.4g Hz ',xmax,max_freq);  
            else
                ylabel('Disp/Force (m/N) ');
                out1=sprintf('\n (Disp/Force) = %8.4g  (m/N) at %8.4g Hz ',xmax,max_freq);
            end      
        end
%
        if(iam==2)
            if(iu==1)
                ylabel('Vel/Force (in/sec/lbf)');
                out1=sprintf('\n (Vel/Force) = %8.4g (in/sec/lbf) at %8.4g Hz ',xmax,max_freq);        
            else
                ylabel('Vel/Force (m/sec/N) ');
                out1=sprintf('\n (Vel/Force) = %8.4g (m/sec/N) at %8.4g Hz ',xmax,max_freq); 
            end
        end
%
        if(iam==3)
            if(iu==1)
                ylabel('Accel/Force (G/lbf)');
                out1=sprintf('\n (Accel/Force) = %8.4g (G/lbf) at %8.4g Hz ',xmax,max_freq); 
            else
                ylabel('Accel/Force (m/sec^2/N) '); 
                out1=sprintf('\n (Accel/Force) = %8.4g (m/sec^2/N) at %8.4g Hz ',xmax,max_freq);        
            end
        end
%

        try

            if(iam==6)
                if(iu==1)
                    ylabel('Stress/Force (psi/lbf)');
                    out1=sprintf('\n (Stress/Force) = %8.4g (psi/lbf) at %8.4g Hz ',xmax,max_freq); 
                else
                    ylabel('Stress/Force (Pa/N) '); 
                    out1=sprintf('\n (Stress/Force = %8.4g (Pa/N) at %8.4g Hz ',xmax,max_freq);        
                end
            end    
    
        end

        disp(out1);  

        if(iflag==1)
            set(gca,'xtick',xtt);
            set(gca,'XTickLabel',xTT);
            fmin=min(xtt);
            fmax=max(xtt);
            xlim([fmin,fmax]);    
        end 
    
    end