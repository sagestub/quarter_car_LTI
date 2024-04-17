
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
m2 = M_fsus+M_whl  ;     % suspension mass                        (kg)
k1 = Kf;                 % spring constant of suspension system   (N/m)
k2 = Kt;                 % spring constant of wheel and tire      (N/m)
b1 = Bf;                 % damping constant of suspension system  (N.s/m)
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
bodemag(G_u)
bodemag(G_w)
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
   (b1+k1)/m2 0 -(k1+k2)/m2 -(b1+b2)/m2];

B=[0 0 
   1/m1 0
   0 0
   -1/m2 (b2+k2)/m2];

C=[1 0 0 0
   0 0 1 0];

D=[0 0
   0 0];
x0 = [0 0 0 0]; % initial conditions of state variables
sys=ss(A,B,C,D);

% plot open-loop poles
% figure()
% p = eig(A)
% plot(real(p), imag(p), 'r*'); % plot poles of open-loop system
% title('Open-loop System Poles')

dist = 1000;            % meters test distance
dx = 0.1;               % meters distance increment
spd_kph = 60            % kilometers per hour test speed
wlen = dist/dx;         % calculate length of road vector
w = func_roadElevationProfile(6, dist,dx,'figure',false,'fignum',6); %road profile vector
u = zeros(1,numel(w));  % create empty input vector
input = [u; w];

spd_mps = spd_kph/3.6;  % convert nominal steady-state speed of simulation to meters/sec
t_end = dist/spd_mps;   % final time value
dt = t_end/size(w,2);   % calculate time step increment
t = 0:dt:dist/spd_mps;  % create time vector
t = t(1:end-1);         % clip time vector size to match road profile vector
y = lsim(sys,input,t,x0);  % simulate system

%% animate model simulation
addpath("qcar_animation")
z0 = w; % road elevation
z1 = y(:,2);              % wheel cm position
z2 = y(:,1);          % sprung mass position
zmf = 5;              % exaggerate response for better visualization
umf = 1;              % road scaling factor
road_z = w;
road_x = 0:dx:dist;
for i=1:length(t)
    plotsusp([z0(i), z1(i)*zmf, z2(i)*zmf, t(i)],road_x,road_z,road_x(i),umf);
    refresh
end

%% plot FFT of qcar response to road input

figure(7)
hold on;
[freq, amp, phase] = simpleFFT(y,dt);
plot(freq, amp);
title('FFT frequency profile of unsprung mass acceleration')

%% launch VDV calculation app
cd 'vibrationdata_9_23_2021'
data = [t' y] %use this to calculate VDV w/ settings:
% m/sec^2, no mean removal, Wk z-axis weighting, no segmentation
% start time = 0, end time = t(end)
run('vibrationdata.m')



