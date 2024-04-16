%
%   Sutherland_air_one_atm.m  ver 1.0  by Tom Irvine
%
%   https://en.wikipedia.org/wiki/Viscosity
%   https://en.wikipedia.org/wiki/Density_of_air

function[dvisc,kvisc]=Sutherland_air_one_atm(tempC)

T0=291.15;      % K
C=120;          % K
mu0=18.27e-06;  % N-sec/m^2 

tempK=tempC+273.15;

dvisc=mu0*((T0+C)/(tempK+C))*((tempK/T0)^1.5);

% table 4 columns
%
%  1  Temperature         T deg C
%  2  Speed of Sound      c (m/sec)
%  3  Air Density         rho (kg/m^3)
%  4  Acoustic Impedance  Z (Pa-sec/m)

table=[35	351.88	1.1455	403.2
30	349.02	1.1644	406.5
25	346.13	1.1839	409.4
20	343.21	1.2041	413.3
15	340.27	1.2250	416.9
10	337.31	1.2466	420.5
5	334.32	1.2690	424.3
0	331.30	1.2922	428.0
-5	328.25	1.3163	432.1
-10	325.18	1.3413	436.1
-15	322.07	1.3673	440.3
-20	318.94	1.3943	444.6
-25	315.77	1.4224	449.1];

[~, index] = min(abs(tempC-table(:,1)));

rho=table(index,3);

kvisc=dvisc/rho;
