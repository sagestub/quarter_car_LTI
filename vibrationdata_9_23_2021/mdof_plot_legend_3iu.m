%
%  mdof_plot_legend_3iu.m  ver 1.1   Tom Irvine
%
function[]=mdof_plot_legend_3iu(t,x,y,accel,num,iu)
%
    [colororder,color_rows]=line_colors();
%
    try
        close 1
    end
    figure(1);
    hold('all');
%
    j=1;
    for i = 1:num
        if(j>color_rows)
            j=j-color_rows;
        end
        out1=sprintf('dof %d',i);
        plot(t,x(:,i),'Color', colororder(j,:),'DisplayName',out1);
        legend('-DynamicLegend');
        j=j+1;
    end
    hold off;
    grid on;
    xlabel('Time(sec)');
    title('Displacement');
    if(iu==1)
        ylabel('Disp(in)');
    else
        ylabel('Disp(m)');    
    end
    plot_legend(num);
    hold off;
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

    try
        close 2
    end
    figure(2);
    hold('all');
%
    j=1;
    for i = 1:num
        if(j>color_rows)
            j=j-color_rows;
        end
        out1=sprintf('dof %d',i);
        plot(t,y(:,i),'Color', colororder(j,:),'DisplayName',out1);
        legend('-DynamicLegend');
        j=j+1;
    end
    hold off;
    grid on;
    xlabel('Time(sec)');
    title('Velocity');
    if(iu==1)
        ylabel('Vel(in/sec)');
    else
        ylabel('Vel(m/sec)');    
    end
    plot_legend(num);
%
    hold off;
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

    try
        close 3
    end
    figure(3);
    N=length(t);
%
    if(iu<=2)
        accel=accel/386;
    end
%
%%
    hold('all');
%
    j=1;
    for i = 1:num
        if(j>color_rows)
            j=j-color_rows;
        end
        out1=sprintf('dof %d',i);
        plot(t,accel(:,i),'Color', colororder(j,:),'DisplayName',out1);
        legend('-DynamicLegend');
        j=j+1;
    end
    hold off;
%%
    xlabel('Time(sec)');
    title('Acceleration');
    if(iu==1)
        ylabel('Accel(G)');
    else
        ylabel('Accel(m/sec^2)');    
    end
    grid on;
    hold off;
%