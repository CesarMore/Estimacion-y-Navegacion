function dxdt = Modelo13_4(t,x,w,wp)
zeta = 0.1;
omegan = 2;
b = -0.4;
dxdt(1,1) = x(2);
dxdt(2,1) = x(3)*x(1) + b*x(2) - x(3)*w;
dxdt(3,1) = wp;