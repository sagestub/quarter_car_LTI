%
%   interpolate_PSD_arbitary_frequency.m  ver 1.0   February 11, 2015
%
function[fi,ai] = interpolate_PSD_arbitary_frequency(f,a,new_freq)
%
    if(f(1) < .0001)
        f(1)=[];
        a(1)=[];
    end
%
    m=length(f);
%
%   recalculate slope
%
    for i=1:m-1
        s(i)=log(  a(i+1) / a(i)  )/log( f(i+1) / f(i) );
    end    
%
    np = length(new_freq);
%
    fi=zeros(np,1);
    ai=zeros(np,1);
%
    for  i=1:np 
%       
		fi(i)=new_freq(i); 
%
        for j=1:(m-1)
%
            if(fi(i)==f(j))
                ai(i)=a(i);
                break;
            end
%
			if( ( fi(i) >= f(j) ) && ( fi(i) <= f(j+1) )  )
				ai(i)=a(j)*( ( fi(i) / f(j) )^ s(j) );
				break;
            end
        end
%               
    end
    nn=length(fi);
    ai(nn)=a(m);