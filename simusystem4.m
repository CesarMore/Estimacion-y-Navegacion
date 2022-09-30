clc
clear all
close all

ti = 0;
T = 0.01; 
tf = 10;

t = (ti:T:tf)';

u = [0;1];
d = 0;
x = [2 1];
m  = 0.1;
mc = 1;
l  = 0.5;
g  = 9.81;
c = [1 1]';
f = [x(2), (g*sin(x(1))-(m*l*(x(2))^2*cos(x(1))*sin(x(1)))/(mc+m))/(l*(4/3-m*cos(x(1))^2)/(mc+m))];
g1 = [0, (cos(x(1))/(mc+m)) / (l*(4/3-m*cos(x(1))^2)/(mc+m))];
g2 = [0 1];
lx = [0 1]';
px = x(2);
e1 = [x(1) - 0];
e2 = [x(2) - 0];
e  = [e1 e2];
sigma = c*e;
d     = 0.2*sin(x(1)*x(2));
xdp = [0 0];
v   = 0.05;
dh  = -1.0;
alfa = 1;
U = 0;
% if norm(sigma/v) <= 1
%     u = -inv(c*g1)*(c*f-c*xdp+c*g2*dh+alfa*(sigma/v)+(norm(sigma*c*g2))^2);
% else
%     u = -inv(c*g1)*(c*f-c*xdp+c*g2*dh+alfa*sign(sigma/v)+(norm(sigma*c*g2))^2);
% end

for i = 1 : length(t)-1
    [tt, xx] = ode45(@modelsystem3, [t(i) t(i+1)], x(i), [], u(i), d);
    x(i+1)   = xx(end);

    U(i+1)   = inv(-c*g1)*(c*f-c*xdp+c*g2*dh+alfa*sign(sigma/v)+(norm(sigma*c*g2))^2);
end 