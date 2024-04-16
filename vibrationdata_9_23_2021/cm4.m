%
%   cm4.m  ver 1.0  December 6, 2013
%
function[m4]=cm4(apsd,fn,nnn)
%
    ra=0.;
%
    for i=1:nnn
        apsd(i)=apsd(i)*fn(i)^4;
    end
%      
	for i=1:(nnn-1)
%
        s=log( apsd(i+1)/apsd(i) )/log( fn(i+1)/fn(i) );
%
        if(s < -1.0001 ||  s > -0.9999 )
            ra=ra+ ( apsd(i+1) * fn(i+1)- apsd(i)*fn(i))/( s+1.);
        else
            ra=ra+ apsd(i)*fn(i)*log( fn(i+1)/fn(i));
        end
%       
	end
%    
    m4=ra;