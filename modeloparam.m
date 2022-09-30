%-------Kalman el del los parametros---------

function dxdt = modeloparam(t,x,w,wp)

b = -0.4;

dxdt(1,1) = x(2);
dxdt(2,1) = x(3)*x(1) + b*x(2) - x(3)*w;
dxdt(3,1) = wp;





















