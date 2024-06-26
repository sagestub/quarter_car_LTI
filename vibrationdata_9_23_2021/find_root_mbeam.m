%
%  beam_type =  1 clamped_pinned_free
%               2 clamped_pinned_clamped
%
%
function[beta_L]=find_root_mbeam(x,a,b,beam_type)
%
k=1;
%
for(i=1:length(x))
%
    iflag=0;
    if(i>1 & k>=2 & x(i)<= 1.01*beta_L(k-1))
        beta=1.01*x(i);
    else
        beta=x(i);
    end    
    beta_a=beta*a;
    beta_b=beta*b;
%
    if(beam_type==1)
        [y(i)]=cpf_C(beta_a,beta_b);
    end
    if(beam_type==2)   
        [y(i)]=cpc_C(beta_a,beta_b);        
    end
%    
    YS(i)=y(i);
    XS(i)=x(i);
%
    if(i>1)
        if(y(i)==0)
            beta_L(k)=x(i);
            if(k>=5)
                break;
            end
            k=k+1;           
        end
        if(y(i-1)*y(i)<0)
%
            x1=x(i-1);
            x2=x(i);
            y1=y(i-1);
            y2=y(i);
%  
            delta=x2-x1;
            if(abs(delta)<1.0e-08)
                  beta_L(k)=(x2+x1)/2;
                  k=k+1;
                  iflag=1;
                  break;
            end
%               
            m=(y2-y1)/delta;
            r=y1-m*x1;
            xx=-r/m;
%
            if(iflag==0)
%
                for j=1:40  % secant method
%
                    beta=xx;
                    beta_a=beta*a;
                    beta_b=beta*b;
 %                   
                    if(beam_type==1)
                        [yy]=cpf_C(beta_a,beta_b);
                    end
                    if(beam_type==2)   
                        [yy]=cpc_C(beta_a,beta_b);        
                    end  
%
                    if(y(i)==0)
                        break;         
                    end              
%
                    if((yy*y1)<0)
                        x2=beta;
                        y2=yy; 
                    else
                        x1=beta;
                        y1=yy;                   
                    end
%
                    delta=x2-x1;
                    if(abs(delta)<1.0e-08)
                        beta_L(k)=(x2+x1)/2;
                        k=k+1;
                        iflag=1;
                        break;
                    end
                    m=(y2-y1)/delta;
                    r=y1-m*x1;
                    xx=-r/m;
%
                end % secant loop
%
                if(iflag==0)
                    beta=beta_a/a;
                    beta_L(k)=xx;
                    if(k>=5)
                        break;
                    end
                    k=k+1; 
%           
                end
            end % if(iflag==0)
%            
        end % if(y(i-1)*y(i)<0)
%        
    end % if i>1
%
end  % length