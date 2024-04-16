clc;
clear all;
close all;
% Constant 'k' Value as per Road roughness ISO 8608
k1=3;
k2=4;
k3=5;
k4=6;
k5=7;
k6=8;
k7=9;
% there is some formula for no. of data point check once
N = 2500; % Number of data points
L = 250; % Length Of Road Profile (m)
B = L/N ; % Sampling Interval (m)
dn = 1/L; % Frequency Band
n0 = 0.1; % Spatial Frequency (cycles/m)
n = dn:dn:dn*N; % Spatial Frequency Band
phi = 2*pi*rand(size(n)); % Random Phase Angle
%here for amp. we can use for loop.
Amp11 = sqrt(dn)*(2^k1)*(1e-3)*(n0./n);% Amplitude for Road Class A-B
Amp12 = sqrt(dn)*(2^k2)*(1e-3)*(n0./n);% Amplitude for Road Class B-C
Amp13 = sqrt(dn)*(2^k3)*(1e-3)*(n0./n);% Amplitude for Road Class C-D
Amp14 = sqrt(dn)*(2^k4)*(1e-3)*(n0./n);
Amp15 = sqrt(dn)*(2^k5)*(1e-3)*(n0./n);
Amp16 = sqrt(dn)*(2^k6)*(1e-3)*(n0./n);
Amp17 = sqrt(dn)*(2^k7)*(1e-3)*(n0./n);
x = 0:B:L-B; % Abscissa Variable from 0 to L
%no need of so many zeros it will takemore space
hx = zeros(size(x));
% hx2 = zeros(size(x));
% hx3 = zeros(size(x));
% hx4 = zeros(size(x));
% hx5 = zeros(size(x));
% hx6 = zeros(size(x));
% hx7 = zeros(size(x));
for i=1:length(x)
    hx1(i) = sum(Amp11.*cos(2*pi*n*x(i)+ phi));
    hx2(i) = sum(Amp12.*cos(2*pi*n*x(i)+ phi));
    hx3(i) = sum(Amp13.*cos(2*pi*n*x(i)+ phi));
    hx4(i) = sum(Amp14.*cos(2*pi*n*x(i)+ phi));
    hx5(i) = sum(Amp15.*cos(2*pi*n*x(i)+ phi));
    hx6(i) = sum(Amp16.*cos(2*pi*n*x(i)+ phi));
    hx7(i) = sum(Amp17.*cos(2*pi*n*x(i)+ phi));
end
%figure is waste command here which is creating an issue i guess so i just
%made it comment from your code
% figure
hold on
%while using hold on command specify all the variable which are to be plotted in between hold on and hold of command as shown
plot(x, hx1);
plot(x, hx2);
plot(x, hx3);
plot(x, hx4);
plot(x, hx5);
plot(x, hx6);
plot(x, hx7);
xlabel('Distance (m)');
ylabel('displacement Of PSD (m)');
grid on
hold off
%no more issues i have found i am also working on PSD if any solution or
%issue do contact me on my university id 1032200009@mitwpu.edu.in