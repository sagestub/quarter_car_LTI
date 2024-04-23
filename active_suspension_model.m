
clear all
close all
clc

% Constants
lb2kg = 1/2.2; % conversion from lbs to kg

% Vehicle Profile
M = 2825.6;          % vehicle curb weight in (kg)
M_fsus = 47.18;      % front suspension unsprung mass (minus tire/wheel) (kg)
M_rsus = 51.03;      % rear suspension unsprung mass (minus tire/wheel) (kg)
M_whl = 39;          % tire/wheel mass (kg)
MR_f = 0.6;          % front spring motion ratio
MR_r = 1;            % rear spring motion ratio
Kf = 123510/MR_f;    % front spring rate (N/m)
Kr = 58620/MR_r;     % rear spring rate (N/m)
Bf = 1200/MR_f;     % rear approximated damping rate (N-s/m)
Br = 1150/MR_r;     % rear approximated damping rate (N-s/m)
Kt = 317400;         % tire spring rate (N/m)
Bt = 0;              % tire damping rate (N-s/m)


% model front 1/4 car
m1 = M/4;                % 1/4 bus body mass                      (kg)
m2 = M_rsus+M_whl  ;     % suspension mass                        (kg)
k1 = Kr;                 % spring constant of suspension system   (N/m)
k2 = Kt;                 % spring constant of wheel and tire      (N/m)
b1 = Br;                 % damping constant of suspension system  (N.s/m)
b2 = Bt;                 % damping constant of wheel and tire     (N.s/m)
%U                      %control force

s = tf('s');
delta = ((m1*s^2+b1*s+k1)*(m2*s^2+(b1+b2)*s+(k1+k2))-(b1*s+k1)*(b1*s+k1));

G_u = ((m1+m2)*s^2+b2*2+k2)/delta; %control input transfer function
G_w = (-m1*b2*s^3+m1*k2*s^2)/delta; %disturbance input transfer function

%% Plot transfer function characteristics
figure()
% step(G_u)


figure()
% step(G_w)


figure()
hold on;
options = bodeoptions;
options.FreqUnits = 'Hz';
bodemag(G_u,options);
bodemag(G_w,options);
title('Frequency response to open-loop inputs')
legend('control input response','road disturbance response')

%% open-loop state-space model
% inputs: U (active forcer), W (road profile)
% states:
% x(0) = position of sprung mass (body)
% x(1) = velocity of sprung mass (body)
% x(2) = position of unsprung mass (wheel)
% x(3) = velocity of unsprung mass (wheel)
A=[0 1   0  0
  -k1/m1 -b1/m1 k1/m1 b1/m1
   0 0 0 1
   k1/m2 b1/m2 -(k1+k2)/m2 -b1/m2];

B=[0     0 
   1/m1  0
   0     0
   -1/m2 k2/m2];

C=[1 0 0 0
   0 1 0 0
   0 0 1 0];

D=[0 0
   0 0
   0 0];
sys=ss(A,B,C,D);
controllability = ctrb(sys);
rank_ctrb = rank(controllability)

observability = obsv(sys);
rank_obsv = rank(observability)

cl = 0;
% plot open-loop poles
% figure()
% p = eig(A)
% plot(real(p), imag(p), 'r*'); % plot poles of open-loop system
% title('Open-loop System Poles')

%% Create closed-loop state-space model
% inputs: U (active forcer), W (road profile)
% states:
% x(0) = position of sprung mass (body)
% x(1) = velocity of sprung mass (body)
% x(2) = position of unsprung mass (wheel)
% x(3) = velocity of unsprung mass (wheel)

% R needs to match number of inputs
% Q needs to match number of states
R_k = 1
Q_k = [1e5 0 0 0
       0 1e5 0 0
       0 0 1e5 0
       0 0 0 0];
fprintf('Q to R = %.2s',Q_k(1,1)/R_k)
% [M, K, L] = icare(A,B,Q_k,R_k) 
[K, ~, ~] = lqr(ss(A,B(:,1),C,D(:,1)),Q_k,R_k)

sys = ss(A-B(:,1)*K(1,:),B(:,2),C,D(:,2));
cl = 1;


%% run model on road
model = "iso";

if model == "iso"
    newRoad = 0;
    spd_kph = 60;           % kilometers per hour test speed
    if newRoad % load road profile if none already created
        dist = 250;             % meters test distance
        dx = 0.1;               % meters distance incrementwlen = dist/dx;         % calculate length of road vector
        w = func_roadElevationProfile(6, dist,dx,'figure',false,'fignum',6); %road profile vector
        w = w-w(1);
    else % use saved road profile
        load("Class6_roadProfile.mat");
    end
    
    if cl
        input = [w];
    else
    u = zeros(1,numel(w));  % create empty input vector 
    input = [u; w];
    end
    
    x0 = [0 0 0 0]; % initial conditions of state variables
    
    spd_mps = spd_kph/3.6;  % convert nominal steady-state speed of simulation to meters/sec
    t_end = dist/spd_mps;   % final time value
    dt = t_end/size(w,2);   % calculate time step increment
    t = 0:dt:dist/spd_mps;  % create time vector
    t = t(1:end-1);         % clip time vector size to match road profile vector
    [y t x_closedloop] = lsim(sys,input,t,x0);  % simulate system
    

elseif model == "sine"
    t_end = 20;
    dt = 0.005;
    t = 0:dt:t_end;
    w = 0.1*sin(pi*(2.475*t+0.5).*t);
    spd_kph = 60;              % kilometers per hour test speed
    dist = t_end*spd_kph/3.6;
    dx = dist/size(t,2);
    if cl
        input = [w];
    else
        u = zeros(1,numel(w));  % create empty input vector 
        input = [u; w];
    end
    x0 = [0 0 0 0]; % initial conditions of state variables
    [y t x_closedloop] = lsim(sys,input,t,x0);

elseif model == "bump"
    t_end = 5;
    dt = 0.005;
    t = 0:dt:t_end;
    i = 0;
    w = ones(1,numel(t));
    for idx = 0:0.005:t_end
        i = i+1;
        if idx>1 && idx<1.5
             w(i) = 0.07*(1-cos(4*pi*idx));
        elseif idx>2.5&& idx<3
             w(i) = 0.1*(1-cos(4*pi*idx));
        else
             w(i) = 0;
        end
    end
    spd_kph = 60;              % kilometers per hour test speed
    dist = t_end*spd_kph/3.6;
    dx = dist/size(t,2);
    if cl
        input = [w];
    else
        u = zeros(1,numel(w));  % create empty input vector 
        input = [u; w];
    end
    x0 = [0 0 0 0]; % initial conditions of state variables
    [y t x_closedloop] = lsim(sys,input,t,x0);
end


% post-processing
z2dot = [0 diff(y(:,2))'/dt];  % sprung mass acceleration
data = [t z2dot'];

if cl
    u_cl = zeros(size(x_closedloop,1),1);
    for idx= 1:size(x_closedloop,1)
        u_cl(idx) = K(1)*x_closedloop(idx,1)+K(2)*x_closedloop(idx,2)+...
        K(3)*x_closedloop(idx,3)+K(4)*x_closedloop(idx,4);
    end
    
    figure()
    subplot(2,1,1)
    plot(t,z2dot)
    legend("body accel")
    subplot(2,1,2)
    plot(t,u_cl)
    legend("ctrl force")

    fprintf('Max actuator force: %.2f Watts\n', max(u_cl))
end

%% animate model simulation
addpath("qcar_animation")
saveVid = 1;
filename = 'qcar_anim.mp4';
v = VideoWriter(filename,'MPEG-4');
v.FrameRate = 60;

% z0 = x(1);                   % road elevation
% z1 = x(2);                   % unsprung mass m deviation
% z2 = x(3);                   % sprung mass m deviation
% t = x(4);                    % current time

z0 = w;                        % road elevation
z1 = y(:,3);                   % wheel m position
z2 = y(:,1);                   % body m position
zmf = 1;                       % exaggerate response for better visualization
umf = 1;                       % road scaling factor
road_z = w;
road_x = 0:dx:dist;
road_x = road_x(1:end-1);

if saveVid
    open(v);
end
for i=1:length(t)
    plotsusp([z0(i), z1(i)*zmf, z2(i)*zmf, t(i)],road_x,road_z,road_x(i),umf);
    drawnow
    if saveVid && mod(i,3)==0
        frame = getframe(gcf);
        writeVideo(v,frame);
    end
end
if saveVid
    close(v);
end

%% plot FFT of qcar response to road input

%plot body response
figure()
hold on;
[freq, amp, phase] = simpleFFT(y(:,1),dt);
loglog(freq, amp);
set(gca,'xscale','log','yscale','log');
grid on;
title('FFT frequency profile of sprung mass acceleration')

%% plot FFT of road input

%plot body response
figure()
hold on;
[freq, amp, phase] = simpleFFT(w,dt);
plot(freq, amp);
% set(gca,'xscale','log','yscale','log');
grid on;
title('FFT frequency profile of sprung mass acceleration')


%% launch VDV calculation app
addpath('vibrationdata_9_23_2021\')
% calculate VDV w/ settings:
% Time History w/ Acceleration input, ISO Generic --> ISO 2631
% m/sec^2, no mean removal, Wk z-axis weighting, no segmentation
% start time = 0, end time = t(end)
run('vibrationdata.m')

%% plot displacement of road, wheel, bodyo

figure()
plot(t,w)
hold on
plot(t,y(:,3))
plot(t,y(:,1))
legend('road','wheel','body')
xlabel('time (sec)')
ylabel('displacement (m)');


%% visualize all bump body accelerations
%run after executing open-loop bump simulation:
plot(t,w.*10,'r:','LineWidth',1.2)
hold on;
plot(t,z2dot,'^-','MarkerIndices',5:10:length(z2dot),'LineWidth',0.7)
load('trial_1.mat')
load('trial_2.mat')
load('trial_3.mat')
load('trial_4.mat')
load('trial_5.mat')
load('trial_6.mat')
load('trial_7.mat')
load('trial_8.mat')
load('trial_9.mat')
load('trial_10.mat')
plot(t,bump_trial_1_z2dot,'o--','MarkerIndices',1:10:length(bump_trial_1_z2dot))
plot(t,bump_trial_2_z2dot,'--')
plot(t,bump_trial_3_z2dot)
plot(t,bump_trial_4_z2dot,":",'LineWidth',2)
plot(t,bump_trial_5_z2dot,'r')
plot(t,bump_trial_6_z2dot,'c')
plot(t,bump_trial_7_z2dot,'m')
plot(t,bump_trial_8_z2dot,'g')
plot(t,bump_trial_9_z2dot,'b--')
plot(t,bump_trial_10_z2dot,'k')
legend('Road Profile (scaled x10 for visibility)','open-loop baseline','trial 1','trial 2', ...
    'trial 3','trial 4','trial 5','trial 6','trial 7','trial 8', ...
    'trial 9','trial 10');
ylabel('Sprung Mass Acceleration (m/s^2)')
xlabel('Time (sec)')
title('Effect of LQR weights on Vehicle Body Acceleration Over Speed Bump')

%% visualize all iso body accelerations
%run after executing open-loop bump simulation:
plot(t,w,'r:')
hold on;
plot(t,z2dot)
load('trial_1.mat')
load('trial_2.mat')
load('trial_3.mat')
load('trial_4.mat')
load('trial_5.mat')
load('trial_6.mat')
load('trial_7.mat')
load('trial_8.mat')
load('trial_9.mat')
load('trial_10.mat')
plot(t,iso_trial_1_z2dot,'o')
plot(t,iso_trial_2_z2dot,'--')
plot(t,iso_trial_3_z2dot)
plot(t,iso_trial_4_z2dot)
plot(t,iso_trial_5_z2dot,'r')
plot(t,iso_trial_6_z2dot,'c')
plot(t,iso_trial_7_z2dot,'m')
plot(t,iso_trial_8_z2dot,'g')
plot(t,iso_trial_9_z2dot,'b')
plot(t,iso_trial_10_z2dot,'k')
legend('Road Profile','open-loop baseline','trial 1','trial 2', ...
    'trial 3','trial 4','trial 5','trial 6','trial 7','trial 8', ...
    'trial 9','trial 10');
ylabel('Sprung Mass Acceleration (m/s^2)')
xlabel('Time (sec)')
title('Effect of LQR weights on Vehicle Body Acceleration Over ISO Grade F Road')

%% Analyze trends in Q and R weight effects on VDV reduction

Q_R = [1,1e2, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12, 1e13];
iso_VDV_raw = [7.568,7.568, 7.647, 7.552, 7.012, 6.141, 5.547, 4.643, 3.072, 1.755, 0.9851, 0.4301];
iso_VDV_reduction = (7.658-iso_VDV_raw)./7.658.*100;
bump_VDV_raw = [10,10, 9.98,9.776, 8.449, 5.631, 3.26, 1.507, 0.5182, 0.1787, 0.06362, 0.02212];
bump_VDV_reduction = (10-bump_VDV_raw)./10.*100;

figure()
hold on;
set(gca,'XScale','log');

scatter(Q_R,iso_VDV_reduction,'o','Color','#0072BD')
% Spline interpolation for iso_VDV_reduction
spline_iso_VDV = spline(log10(Q_R), iso_VDV_reduction);
q_interp = logspace(log10(min(Q_R)), log10(max(Q_R)), 100); % Generate finer sampling points for smoother spline plot
iso_VDV_interp = ppval(spline_iso_VDV, log10(q_interp)); % Evaluate the spline at the finer sampling points
plot(q_interp, iso_VDV_interp,'Color','#0072BD'); % Plot the spline curve

scatter(Q_R,bump_VDV_reduction,'g^')
% Spline interpolation for bump_VDV_reduction vs. Q_R (logarithmic Q_R)
spline_bump_VDV = spline(log10(Q_R), bump_VDV_reduction);
q_interp = logspace(log10(min(Q_R)), log10(max(Q_R)), 100); % Generate finer sampling points for smoother spline plot
bump_VDV_interp = ppval(spline_bump_VDV, log10(q_interp)); % Evaluate the spline at the finer sampling points
plot(q_interp, bump_VDV_interp,'g'); % Plot the spline curve

legend('ISO Grade F Road Case', 'ISO Grade F Spline Fit','Speed Bump Case','Speed Bump Spline Fit');
xlabel('Q_R');
ylabel('Percent VDV Reduction');
title('Effect of LQR Weights on VDV Reduction Trend');

%% Analyze trend between VDV reduction and power consumption

iso_VDV_raw = [7.568,7.568, 7.647, 7.552, 7.012, 6.141, 5.547, 4.643, 3.072, 1.755, 0.9851, 0.4301];
iso_VDV_reduction = (7.658-iso_VDV_raw)./7.658.*100;
iso_pwr_kW = [6.7e-7, 6.7e-5, 6.67e-3, 0.064, 0.4726, 1.651, 3.514, 6.224, 9.168, 16.30, 24.14, 27.96];

bump_VDV_raw = [10,10, 9.98,9.776, 8.449, 5.631, 3.26, 1.507, 0.5182, 0.1787, 0.06362, 0.02212];
bump_VDV_reduction = (10-bump_VDV_raw)./10.*100;
bump_pwr_kW = [1.2e-6, 1.2e-4, 1.19e-1, 0.112, 0.761, 2.404, 5.48, 8.441, 10.57, 11.71, 12.20, 12.41];

figure()
hold on;
plot(iso_pwr_kW,iso_VDV_reduction)
plot(bump_pwr_kW,bump_VDV_reduction)

xlabel('Peak Power Required (kW)');
ylabel('Percent VDV Reduction From Passive Baseline')
title('Power Cost of Active Suspension Performance Measured in VDV Reduction')
legend('ISO Grade F Road Case','Speed Bump Case')


