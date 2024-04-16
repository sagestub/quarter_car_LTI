%
%   xtick_label.m  ver 2.9  by Tom Irvine
%
function[xtt,xTT,iflag]=xtick_label(fmin,fmax)


try

xtt=[];
xTT={};
iflag=0;



%
%  start 0.001
%
if(fmin>=0.001 && fmin<0.005)

    if(fmax<=0.2)
        xtt=[0.001 0.01 0.02 ];
        xTT={'0.001';'0.01';'0.02'};
        iflag=1; 
        return;
    end
    if(fmax<=0.1)
        xtt=[0.001 0.01 0.1 ];
        xTT={'0.001';'0.01';'0.1'};
        iflag=1; 
        return;
    end
    if(fmax<=0.2)
        xtt=[0.001 0.01 0.1 0.2];
        xTT={'0.001';'0.01';'0.1';'0.2'};
        iflag=1; 
        return;
    end
    if(fmax<=0.5)
        xtt=[0.001 0.01 0.1 0.5];
        xTT={'0.001';'0.01';'0.1';'0.5'};
        iflag=1; 
        return;
    end
    if(fmax<=1)
        xtt=[0.001 0.01 0.1 1];
        xTT={'0.001';'0.01';'0.1';'1'};
        iflag=1; 
        return;
    end
    if(fmax<=5)
        xtt=[0.001 0.01 0.1 1 5];
        xTT={'0.001';'0.01';'0.1';'1';'5'};
        iflag=1; 
        return;
    end    
end


%
%  start 0.01
%
if(fmin>=0.01 && fmin<0.03)

    if(fmax<=0.2)
        xtt=[0.01 0.1 0.2 ];
        xTT={'0.01';'0.1';'0.2'};
        iflag=1; 
        return;
    end
    if(fmax<=1)
        xtt=[0.01 0.1 1 ];
        xTT={'0.01';'0.1';'1'};
        iflag=1; 
        return;
    end
    if(fmax<=2)
        xtt=[0.01 0.1 1 2];
        xTT={'0.01';'0.1';'1';'2'};
        iflag=1; 
        return;
    end
    if(fmax<=5)
        xtt=[0.01 0.1 1 5];
        xTT={'0.01';'0.1';'1';'5'};
        iflag=1; 
        return;
    end

end


%
%  start 0.1
%
if(fmin>=0.1 && fmin<0.3)

    if(fmax<=2)
        xtt=[0.1 1 2 ];
        xTT={'0.1';'1';'2'};
        iflag=1; 
        return;
    end
    if(fmax<=10)
        xtt=[0.1 1 10 ];
        xTT={'0.1';'1';'10'};
        iflag=1; 
        return;
    end
    if(fmax<=20)
        xtt=[0.1 1 10 20];
        xTT={'0.1';'1';'10';'20'};
        iflag=1; 
        return;
    end
    if(fmax<=50)
        xtt=[0.1 1 10 50];
        xTT={'0.1';'1';'10';'50'};
        iflag=1; 
        return;
    end

end

%
%
%  start 0.3
%
if(fmin>=0.3 && fmin<1)

    if(fmax<=20)
        xtt=[0.3 1 10 20];
        iflag=1; 
        return;
    end
    if(fmax<=50)
        xtt=[0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 2 3 4 5 6 7 8 9 10 20 30 40 50];
        xTT={'0.3';'';'';'';'';'';'';'1';'';'';'';'';'';'';'';'';'10';'';'';'';'50'};
        iflag=1; 
        return;
    end
    if(fmax<=100)
        xtt=[0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 2 3 4 5 6 7 8 9 10 20 30 40 50 60 70 80 90 100];
        xTT={'0.3';'';'';'';'';'';'';'1';'';'';'';'';'';'';'';'';'10';'';'';'';'';'';'';'';'';'100'};
        iflag=1; 
        return;
    end
end

%
%
%  start 1
%

if(fmin>=1 && fmin<2)

    if(fmax<=20)
        xtt=[1 10 20];
        xTT={'1';'10';'20'};
        iflag=1; 
        return;
    end
    if(fmax<=50)
        xtt=[1 10 20 30 40 50];
        xTT={'1';'10';'20';'30';'40';'50'};
        iflag=1; 
        return;
    end  
    if(fmax<=100)
        xtt=[1 2 3 4 5 6 7 8 9 10 20 30 40 50 60 70 80 90 100];
        xTT={'1';'';'';'';'';'';'';'';'';'10';'';'';'';'';'';'';'';'';'100'};
        iflag=1; 
        return;
    end  

end



%
%  start 2
%


if(fmin>=2 && fmin<5)

    if(fmax<=80)
        xtt=[2 3 4 5 6 7 8 9 10 20 30 40 50 60 70 80];
        xTT={'2';'';'';'';'';'';'';'';'10';'';'';'';'';'';'';'80'};
        iflag=1; 
        return;
    end
end


%
%  start 5
%

if(fmin>=5 && fmin<10)

    if(fmax<51)
        xtt=[5 6 7 8 9 10 20 30 40 50];
        xTT={'5';'';'';'';'';'10';'';'';'';'50'};
        iflag=1; 
        return;
    end
    if(fmax<101)
        xtt=[5 6 7 8 9 10 20 30 40 50 60 70 80 90 100];
        xTT={'5';'';'';'';'';'10';'';'';'';'';'';'';'';'';'100'};
        iflag=1; 
        return;
    end
    if(fmax<510)
        xtt=[5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
        xTT={'5';'';'';'';'';'10';'';'';'';'';'';'';'';'';'100';'';'';'';'500'};
        iflag=1; 
        return;
    end
    if(fmax<1010)
        xtt=[5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500 600 700 800 900 1000];
        xTT={'5';'';'';'';'';'10';'';'';'';'';'';'';'';'';'100';'';'';'';'500';'';'';'';'';'1000'};
        iflag=1; 
        return;
    end
    if(fmax<=2050)
        xtt=[5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500 600 700 800 900 1000 2000];
        xTT={'5';'';'';'';'';'10';'';'';'';'';'';'';'';'';'100';'';'';'';'500';'';'';'';'';'1000';'2000'};
        iflag=1; 
        return;
    end

end


%
%  start 10
%

if(fmin>=10 && fmin<20)

    if(fmax<51)
        xtt=[10 20 30 40 50];
        xTT={'10';'20';'30';'40';'50'};
        iflag=1; 
        return;
    end
    if(fmax<=105)
        xtt=[10 20 30 40 50 60 70 80 90 100];
        xTT={'10';'';'';'';'';'';'';'';'';'100'};
        iflag=1; 
        return;    
    end
    if(fmax<=205)
        xtt=[10 20 30 40 50 60 70 80 90 100 200];
        xTT={'10';'';'';'';'';'';'';'';'';'100';'200'};
        iflag=1; 
        return;    
    end
    if(fmax<=505)
        xtt=[10 20 30 40 50 60 70 80 90 100 200 300 400 500];
        xTT={'10';'';'';'';'';'';'';'';'';'100';'';'';'';'500'};
        iflag=1; 
        return;
    end
    if(fmax<=1020)
        xtt=[10 20 30 40 50 60 70 80 90 100 200 300 400 500 600 700 800 900 1000];
        xTT={'10';'';'';'';'';'';'';'';'';'100';'';'';'';'';'';'';'';'';'1000'};
        iflag=1; 
        return;
    end
    if(fmax<=2050)
        xtt=[10 20 30 40 50 60 70 80 90 100 200 300 400 500 600 700 800 900 1000 2000];
        xTT={'10';'';'';'';'';'';'';'';'';'100';'';'';'';'';'';'';'';'';'1000';'2000'};
        iflag=1; 
        return;
    end
    if(fmax<=3050)
        xtt=[10 20 30 40 50 60 70 80 90 100 200 300 400 500 600 700 800 900 1000 2000 3000];
        xTT={'10';'';'';'';'';'';'';'';'';'100';'';'';'';'';'';'';'';'';'1000';'';'3000'};
        iflag=1; 
        return;
    end
    if(fmax<=4050)
        xtt=[10 20 30 40 50 60 70 80 90 100 200 300 400 500 600 700 800 900 1000 2000 3000 4000];
        xTT={'10';'';'';'';'';'';'';'';'';'100';'';'';'';'';'';'';'';'';'1000';'';'';'4000'};
        iflag=1; 
        return;
    end
    if(fmax<=5100)
        xtt=[10 20 30 40 50 60 70 80 90 100 200 300 400 500 600 700 800 900 1000 2000 3000 4000 5000];
        xTT={'10';'';'';'';'';'';'';'';'';'100';'';'';'';'';'';'';'';'';'1000';'';'';'';'5000'};
        iflag=1; 
        return;
    end
%
    if(fmax<=10010)
        xtt=[10 20 30 40 50 60 70 80 90 100 200 300 400 500 600 700 800 900 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000];
        xTT={'10';'';'';'';'';'';'';'';'';'100';'';'';'';'';'';'';'';'';'1000';'';'';'';'';'';'';'';'';'10K'};
        iflag=1; 
        return;
    end
%
    if(fmax<=20010)
        xtt=[10 20 30 40 50 60 70 80 90 100 200 300 400 500 600 700 800 900 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000 20000];
        xTT={'10';'';'';'';'';'';'';'';'';'100';'';'';'';'';'';'';'';'';'1000';'';'';'';'';'';'';'';'';'10K';'20K'};
        iflag=1; 
        return;
    end

end


%
%  start 20
%

if(fmin>=20 && fmin<50)
    
    if(fmax<=100)
        xtt=[20 30 40 50 60 70 80 90 100];
        xTT={'20';'';'';'';'';'';'';'';'100'};
        iflag=1; 
        return;
    end 
    if(fmax<=200)
        xtt=[20 30 40 50 60 70 80 90 100 200];
        xTT={'20';'';'';'';'';'';'';'';'100';'200'};
        iflag=1; 
        return;
    end       
    if(fmax<=1010)
        xtt=[20 30 40 50 60 70 80 90 100 200 300 400 500 600 700 800 900 1000];
        xTT={'20';'';'';'';'';'';'';'';'100';'';'';'';'';'';'';'';'';'1000'};
        iflag=1; 
        return;
    end    
    if(fmax<=2010)
        xtt=[20 30 40 50 60 70 80 90 100 200 300 400 500 600 700 800 900 1000 2000];
        xTT={'20';'';'';'';'';'';'';'';'100';'';'';'';'';'';'';'';'';'1000';'2000'};
        iflag=1; 
        return;
    end
    if(fmax<=3000)
        xtt=[20 30 40 50 60 70 80 90 100 200 300 400 500 600 700 800 900 1000 2000 3000];
        xTT={'20';'';'';'';'';'';'';'';'100';'';'';'';'';'';'';'';'';'1K';'2K';'3K'};
        iflag=1; 
        return;
    end
    if(fmax<=4000)
        xtt=[20 30 40 50 60 70 80 90 100 200 300 400 500 600 700 800 900 1000 2000 3000 4000];
        xTT={'20';'';'';'';'';'';'';'';'100';'';'';'';'';'';'';'';'';'1000';'';'';'4000'};
        iflag=1; 
        return;
    end
    if(fmax<=5010)
        xtt=[20 30 40 50 60 70 80 90 100 200 300 400 500 600 700 800 900 1000 2000 3000 4000 5000];
        xTT={'20';'';'';'';'';'';'';'';'100';'';'';'';'';'';'';'';'';'1000';'';'';'';'5000'};
        iflag=1; 
        return;
    end
    if(fmax<=10010)
        xtt=[20 30 40 50 60 70 80 90 100 200 300 400 500 600 700 800 900 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000];
        xTT={'20';'';'';'';'';'';'';'';'100';'';'';'';'';'';'';'';'';'1000';'';'';'';'';'';'';'';'';'10K'};
        iflag=1; 
        return;
    end
    if(fmax<=15010)
        xtt=[20 30 40 50 60 70 80 90 100 200 300 400 500 600 700 800 900 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000 15000];
        xTT={'20';'';'';'';'';'';'';'';'100';'';'';'';'';'';'';'';'';'1000';'';'';'';'';'';'';'';'';'10K';'15K'};
        iflag=1; 
        return;
    end
    if(fmax<=20010)
        xtt=[20 30 40 50 60 70 80 90 100 200 300 400 500 600 700 800 900 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000 20000];
        xTT={'20';'';'';'';'';'';'';'';'100';'';'';'';'';'';'';'';'';'1000';'';'';'';'';'';'';'';'';'10K';'20K'};
        iflag=1; 
        return;
    end

end


%
%  start 50
%

if(fmin>=50 && fmin<100)

    if(fmax<=2010)
        xtt=[50 60 70 80 90 100 200 300 400 500 600 700 800 900 1000 2000];
        xTT={'50';'';'';'';'';'100';'';'';'';'';'';'';'';'';'1K';'2K'};
        iflag=1; 
        return;
    end    
    if(fmax<=3010)
        xtt=[50 60 70 80 90 100 200 300 400 500 600 700 800 900 1000 2000 3000];
        xTT={'50';'';'';'';'';'100';'';'';'';'';'';'';'';'';'1K';'2K';'3K'};
        iflag=1; 
        return;
    end
    if(fmax<=4010)
        xtt=[50 60 70 80 90 100 200 300 400 500 600 700 800 900 1000 2000 3000 4000];
        xTT={'50';'';'';'';'';'100';'';'';'';'';'';'';'';'';'1K';'2K';'3K','4K'};
        iflag=1; 
        return;
    end

end

%
%  start 100
%

if(fmin>=100 && fmin<200)

    if(fmax<=2010)
        xtt=[100 200 300 400 500 600 700 800 900 1000 2000];
        xTT={'100';'';'';'';'';'';'';'';'';'1000';'2000'};
        iflag=1; 
        return;
    end

    if(fmax<=3010)
        xtt=[100 200 300 400 500 600 700 800 900 1000 2000 3000];
        xTT={'100';'200';'';'';'500';'';'';'';'';'1000';'2000';'3000'};
        iflag=1; 
        return;
    end

    if(fmax<=4010)
        xtt=[100 200 300 400 500 600 700 800 900 1000 2000 3000 4000];
        xTT={'100';'';'';'';'';'';'';'';'';'1000';'2000';'3000';'4000'};
        iflag=1; 
        return;
    end

    if(fmax<=5010)
        xtt=[100 200 300 400 500 600 700 800 900 1000 2000 3000 4000 5000];
        xTT={'100';'';'';'';'';'';'';'';'';'1000';'2000';'3000';'4000';'5000'};
        iflag=1; 
        return;
    end

    if(fmax<=10010)
        xtt=[100 200 300 400 500 600 700 800 900 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000];
        xTT={'100';'';'';'';'';'';'';'';'';'1000';'';'';'';'';'';'';'';'';'10K'};
        iflag=1; 
        return;
    end

end

%
%  start 200
%

if(fmin>=200 && fmin<=2000)

    if(fmax<=2010)
        xtt=[200 300 400 500 600 700 800 900 1000 2000];
        xTT={'200';'';'';'';'';'';'';'';'1000';'2000'};
        iflag=1; 
        return;
    end
    if(fmax<=3010)
        xtt=[200 300 400 500 600 700 800 900 1000 2000 3000];
        xTT={'200';'';'';'';'';'';'';'';'1000';'2000';'3000'};
        iflag=1; 
        return;
    end   
    if(fmax<=10010)
        xtt=[200 300 400 500 600 700 800 900 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000];
        xTT={'200';'';'';'';'';'';'';'';'1000';'';'';'';'';'';'';'';'';'10K'};
        iflag=1; 
        return;
    end     
    
end


%
%

catch
   
   out1=sprintf(' fmin=%9.5g  fmax=%9.5g  ',fmin,fmax);
   disp(out1);
   
   warndlg('xtick_label error ');
   return;
end

