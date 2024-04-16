    
%  rainflow_damage_function.m  ver 1.0  by Tom Irvine

%  response = response time history amplitude
%         b = fatigue exponent (typically: 4 <= b <= 8 )
%    damage = relative damage, dimension: (response amplitude unit)^b

%   ASTM E1049-85 method
%
%   https://www.mathworks.com/help/signal/ref/rainflow.html
%

function[damage]=rainflow_damage_function(response,b)

c=rainflow(response);

% The output array c has five columns:   Count, Range, Mean, Start, and End
% Only the first two columns are needed for the damage calculation  
% The relative damage can then be calculated for a fatigue exponent b  
% The amplitude is one-half of range.

cycles=c(:,1);
amp=c(:,2)/2;
damage=sum( cycles.*amp.^b );
