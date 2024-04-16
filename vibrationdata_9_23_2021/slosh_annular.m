disp(' ');
disp(' slosh_annular.m   ver 1.1   July 16, 2010 ');
disp(' by Tom Irvine ');
%
disp(' ');
disp(' This program calculates the slosh frequency in ');
disp(' a annular cylinder basin.');
disp(' ');
disp(' Assume');
disp(' 1. The liquid is homogeneous, inviscid, irrotational, & incompressible');  
disp(' 2. The boundaries of the basin are rigid ');   
disp(' 3. Small wave amplitudes, linear behavior');   
disp(' 4. The influence of the surrounding atmosphere is negligible');   
disp(' 5. The influence of surface tension is negligible');   
%
clear f;
clear ff;
clear a;
clear b;
clear B;
clear D;
clear DD;
%
tpi=2*pi;
% 
disp(' ');
disp(' Enter acceleration of gravity (in/sec^2) ');
g=input(' ');
disp(' Enter the liquid surface height (inch)');
h=input(' ');
disp(' Enter the outer diameter (inch)');
D1=input(' ');
disp(' Enter the inner diameter (inch)');
D2=input(' ');
%
R2=D2/2;
R1=D1/2;
%
x=R2/R1;
if(x>1)
    x=1/x;
end
if(R1<R2)
    R1=R2;
end
%
out1=sprintf('%8.4g %8.4g %8.4g %8.4g %8.4g',g,D2,D1,h,x);
disp(out1);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
term1=inline('(besselj(i-1,x)-besselj(i+1,x))/2','i','x');
term2=inline('(bessely(i-1,x)-bessely(i+1,x))/2','i','x');  
%
lambda=zeros(10,10);
for(i=0:6)
   n=1;
   for(kj=1:500)
       k=kj*0.1;
       x=k;
       xba=x*R2/R1;
       CF=term1(i,x)*term2(i,xba)-term1(i,xba)*term2(i,x);
 %
 %      out1=sprintf('kj=%d  x=%8.4g  t1=%8.4g  t2=%8.4g  t3=%8.4g  t4=%8.4g  CF=%8.4g',kj,x,term1,term2,term3,term4,CF);
 %      disp(out1);
 %
       if(CF==0)
           lambda(i+1,n)=k;
           n=n+1;
           if(n==4)
               break;
           end
       else
           if(kj>=2)
              if((CF*CFb)<0.)
                 for(m=1:20)
                     nk=kb-CFb*(k-kb)/(CF-CFb);
                     x=nk;
                     xba=x*R2/R1;
                     CFn=term1(i,x)*term2(i,xba)-term1(i,xba)*term2(i,x);
                     %
                     if(CFn*CF>0.)
                         CF=CFn;
                         k=nk;
                     else
                         CFb=CFn;
                         kb=nk;                         
                     end
                 end
                 lambda(i+1,n)=nk;
                 n=n+1;
                 if(n==11)
                    break;
                 end  
              end
           end
       end
       CFb=CF;
       kb=k;
   end
end   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
kk=1;
for(i=1:6)
     for(j=1:3)
         lam=lambda(i,j)/R1;
         f(i,j)=(1/tpi)*sqrt((lam*g)*tanh(lam*h));
         ii=i-1;
         if(i==1)
             jj=j;
         else
             jj=j-1;
         end
%%         out1=sprintf(' f=%8.4g Hz  lambda=%10.6g  i=%d j=%d',f(i,j),lambda(i,j),ii,jj);
%%         disp(out1);
         ff(kk,1)=f(i,j);
         ff(kk,2)=lambda(i,j);
         ff(kk,3)=ii;
         ff(kk,4)=jj;     
         kk=kk+1;
     end
end    
%
B = sortrows(ff,1);
%
kk=kk-1;
disp(' ');
disp(' i=number of nodal diameters ');
disp(' j=number of nodal circles ');
disp(' ');
disp(' Freq(Hz)      lambda        i   j  ');
for(i=1:kk)
    out1=sprintf(' %8.4g \t %10.6g \t %d \t %d',B(i,1),B(i,2),round(B(i,3)),round(B(i,4)));
    disp(out1);
end