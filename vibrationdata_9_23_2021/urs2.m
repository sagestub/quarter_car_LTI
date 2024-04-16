%
%   urs_peak_test.m  ver 1.1  by Tom Irvine
%

clear fn;
clear n;
clear nx;
clear r;
clear rl;
clear c;
clear T;
clear aa

mu=0;
sigma=1;

sr=10000;
dt=1/sr;
f=1000;

T=1200;

fn=800;

[r]=maximax_peak(fn,T);

arg=fn*T;

nt=round(T*sr);

damp=0.05;

num=100;



nx=zeros(num,1);



num=1000;

ra=zeros(num,1);
re=zeros(num,1);
aa=zeros(num,1);


for i=1:num
    

        aa(i)=i/1000;

        ra(i)=r*sqrt(1-(log(aa(i))/log(arg))); 
        
        [ps]=ERS_peak(arg,aa(i));   
        
        re(i)=ps;
end 


ry=zeros(200,1);
yy=zeros(200,1);

for i=1:200
    ry(i)=0.05*i;
    yy(i)=sum(peak>ry(i))/length(peak);
    
    out1=sprintf(' yy=%8.4g  ra=%8.4g  ',yy(i),ra(i));
    disp(out1);
end

figure(340);
plot(re,aa,ry,yy);
grid on
xlim([4.8 6]);

p_sum1=[re aa];
p_sum2=[ry yy];