%   
%    calculate_PSD_slopes_alt.m,  ver 1.7,  October 25, 2018
%
function [s,grms] = calculate_PSD_slopes_alt(A)
%
f=A(:,1);
a=A(:,2);
%
MAX = 12000;
%
ra=0.;
grms=0.;
iflag=0;
%
s=zeros(1,MAX);
%
if(f(1) < .0001)
    f(1)=[];
    a(1)=[];
end
%
nn=length(f)-1;
%
for  i=1:nn
%
    if(  f(i) <=0 )
        disp(' frequency error ')
        out=sprintf(' f(%d) = %6.2f ',i,f(i));
        disp(out)
        iflag=1;
    end
    if(  a(i) <=0 )
        disp(' amplitude error ')
        out=sprintf(' a(%d) = %6.2f ',i,a(i));
        disp(out)
        iflag=1;
    end  
    if(  f(i+1) < f(i) )
        disp(' frequency error ')
        iflag=1;
    end  
    if(  iflag==1)
        break;
    end
%    
    if(f(i+1)~=f(i))
        s(i)=log10( a(i+1)/ a(i) )/log10( f(i+1)/f(i) );
    else
        s(i)=NaN;
    end
%   
 end
 %
 %disp(' RMS calculation ');
 %
 if( iflag==0)
    for i=1:nn
        if(abs(s(i))<1000)
            if(s(i) < -1.0001 ||  s(i) > -0.9999 )
                ra = ra + ( a(i+1) * f(i+1)- a(i)*f(i))/( s(i)+1.);
            else
                ra = ra + a(i)*f(i)*log( f(i+1)/f(i));
            end
        end
    end
    grms=sqrt(ra);
end
         
