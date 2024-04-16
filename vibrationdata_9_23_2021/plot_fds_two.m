%
%  plot_fds_two.m  ver 1.2  by Tom Irvine
%
function[fig_num]=...
         plot_fds_two(fig_num,x_label,y_label,ppp,qqq,leg_a,leg_b,Q,bex,iu,nmetric)
%
f=ppp(:,1);
a=ppp(:,2);
%
figure(fig_num);
fig_num=fig_num+1;
%
plot(ppp(:,1),ppp(:,2),qqq(:,1),qqq(:,2));
%
legend(leg_a,leg_b);
%
set(gca,'MinorGridLineStyle','none','GridLineStyle',':','XScale','log','YScale','log','XminorTick','off','YminorTick','off');
%
fmax=max(f);
fmin=min(f);


if(fmin<1)
        fmin=1;
end

[xtt,xTT,iflag]=xtick_label(fmin,fmax);


%

[y_label,t_string]=fds_ylabel(Q,bex,nmetric,iu);


out=sprintf(t_string);
title(out);

ylabel(y_label);   
xlabel(x_label);

%
set(gca,'XGrid','on','GridLineStyle',':');
set(gca,'YGrid','on','GridLineStyle',':');
%

if(iflag==1)
        set(gca,'xtick',xtt);
        set(gca,'XTickLabel',xTT);
        fmin=min(xtt);
        fmax=max(xtt);
end   

%
xlim([fmin,fmax]);
%