%
%  srs_coefficients.m  ver 1.5  October 5, 2012
%
function[a1,a2,b1,b2,b3,rd_a1,rd_a2,rd_b1,rd_b2,rd_b3]=...
                                      srs_coefficients(f,damp,dt)    
%
%%% disp(' ')
%%% disp(' Initialize coefficients')
%
%%% disp(' ')
%%% disp(' Select algorithm: ')
%%% disp(' 1=Kelly-Richman  2=Smallwood ');
%%% ialgorithm=input(' ');
%
%  Set to Smallwood for sake of MDOF programs which use this function
%
ialgorithm=2;
%
tpi=2*pi;
%
num_fn=max(size(f));
%
a1=zeros(num_fn,1);
a2=zeros(num_fn,1);
b1=zeros(num_fn,1);
b2=zeros(num_fn,1);   
b3=zeros(num_fn,1);   
%
rd_a1=zeros(num_fn,1);
rd_a2=zeros(num_fn,1);
rd_b1=zeros(num_fn,1);
rd_b2=zeros(num_fn,1);   
rd_b3=zeros(num_fn,1); 
%
num_damp=length(damp);
%
for j=1:num_fn
%
        omega=(tpi*f(j));
%
        if(num_damp==1)
            ddd=damp;
        else
            ddd=damp(j);
        end
%
        omegad=(omega*sqrt(1.-ddd^2));
%            
        cosd=(cos(omegad*dt));
        sind=(sin(omegad*dt));
        domegadt=(ddd*omega*dt);
%
        rd_a1(j)=2.*exp(-domegadt)*cosd; 
        rd_a2(j)=-exp(-2.*domegadt);
        rd_b1(j)=0.;
        rd_b2(j)=-(dt/omegad)*exp(-domegadt)*sind;
        rd_b3(j)=0;
%
        if(ialgorithm==1)
%
            a1(j)=(2.*exp(-domegadt)*cosd);
            a2(j)=(-exp(-2.*domegadt));
            b1(j)=(2.*domegadt);
            b2(j)=(omega*dt*exp(-domegadt));
            b2(j)=b2(j)*(( (omega/omegad)*(1-2*(ddd^2))*sind -2*ddd*cosd ));
            b3(j)=0.;
%           
        else
%
            E=0;
            K=0;
            C=0;
            S=0;
            Sp=0;
%
            E=(exp(-ddd*omega*dt));
            K=(omegad*dt);
            C=(E*cos(K));
            S=(E*sin(K));
%
            Sp=S/K;
%
            a1(j)=(2*C);
            a2(j)=(-(E^2));
%
            b1(j)=(1.-Sp);
            b2(j)=(2.*(Sp-C));
            b3(j)=((E^2)-Sp);
        end
%
end