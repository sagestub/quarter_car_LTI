
% peak_vrs.m  ver 1.0 by Tom Irvine

function[accel_vrs,psvel_vrs,reldisp_vrs]=...
       peak_vrs(fn,a_vrs,pv_vrs,rd_vrs,dur,spec_minf,spec_maxf)
   
j=1;
 
for i=1:length(fn)
%    
    if(fn(i)>=spec_minf && fn(i)<=spec_maxf)
%
        C=sqrt(2*log(fn(i)*dur));
        ms=C + (0.5772/C);
%
           a_vrs_peak=ms*a_vrs(i);
           pv_vrs_peak=ms*pv_vrs(i);        
           rd_vrs_peak=ms*rd_vrs(i);
           
           
           accel_vrs(j,:)  =[ fn(i)  a_vrs(i)   3*a_vrs(i)   a_vrs_peak];
           psvel_vrs(j,:)  =[ fn(i) pv_vrs(i)  3*pv_vrs(i)  pv_vrs_peak];
           reldisp_vrs(j,:)=[ fn(i) rd_vrs(i)  3*rd_vrs(i)  rd_vrs_peak];
           
           j=j+1;
    end   
%        
end   
