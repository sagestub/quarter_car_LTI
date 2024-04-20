
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
R_k = 0.2
Q_k = diag([1 1 1 1])
[M, K, L] = icare(A,B,Q_k,R_k) 

sys = ss(A-B(:,1)*K(1,:),B,C,D);
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
    
    u = zeros(1,numel(w));  % create empty input vector -K(1,:)*x;%
    input = [u; w];
    
    x0 = [w(1) 0 0 0]; % initial conditions of state variables
    
    spd_mps = spd_kph/3.6;  % convert nominal steady-state speed of simulation to meters/sec
    t_end = dist/spd_mps;   % final time value
    dt = t_end/size(w,2);   % calculate time step increment
    t = 0:dt:dist/spd_mps;  % create time vector
    t = t(1:end-1);         % clip time vector size to match road profile vector
    y = lsim(sys,input,t,x0);  % simulate system
elseif model == "sine"
 t_end = 20;
 dt = 0.005
 t = 0:dt:t_end;
 w = 0.1*sin(pi*(2.475*t+0.5).*t);
 spd_kph = 60;              % kilometers per hour test speed
 dist = t_end*spd_kph/3.6;
 dx = dist/size(t,2);
 u = zeros(1,numel(w));     % create empty input vector
 input = [u; w];
 x0 = [w(1) 0 0 0]; % initial conditions of state variables
 y = lsim(sys,input,t,x0);

elseif model == "bump"
 t_end = 5;
 dt = 0.005
 t = 0:dt:t_end;
 i = 0;
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
 u = zeros(1,numel(w));     % create empty input vector
 input = [u; w];
 x0 = [w(1) 0 0 0]; % initial conditions of state variables
 y = lsim(sys,input,t,x0);
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
z2 = y(:,1);
z2dot = [0 diff(y(:,2))'/dt];  % sprung mass acceleration
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
data = [t' z2dot'] %use this to calculate VDV w/ settings:
% Time History w/ Acceleration input, ISO Generic --> ISO 2631
% m/sec^2, no mean removal, Wk z-axis weighting, no segmentation
% start time = 0, end time = t(end)
run('vibrationdata.m')



