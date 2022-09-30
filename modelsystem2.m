function xp = modelsystem2(t,x,d)
%Representacion de estado
A   = [-0.8 1 1.6; 0 -3 2; 0 0 -6];
Bu  = [0 0 ; 1 0; 0 1];
Bd  = [0.8 0; 0 -1; -0.4 1.2];
L   = [40 0 -20; 0 -20 60];
LBd = [-40 24; 24 -92];
u = [1;1];

xp = A*x + Bu*u + Bd*d; 
end