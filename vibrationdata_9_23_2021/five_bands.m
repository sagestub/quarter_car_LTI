

P=[68.3 27.1 4.33 0.26 0.0063]/100;

b=6.4;

A=0;

for i=1:5
    
    A=A+P(i)*i^6.4;
    
end

A^(1/b)