function [dxdt] = ModeloEjem13_2(t,x,w1,w2,w3)
po = 0.0034;
g = 32.2;
k = 22000;

dxdt(1,1) = x(2) + w1;
dxdt(2,1) = (po*(exp(-x(1)/k))*x(2)^22*x(3))/2 - g + w2;
dxdt(3,1) = w3;
end