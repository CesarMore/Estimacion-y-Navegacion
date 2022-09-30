clear all
close all
clc

ti = 0;
h = 0.01;
tf = 5;

ts = (ti:h:tf)';

x0 = [pi;0;0;0];
xh0 = [-0.8;0.1;0.1;0.1];

config = odeset('RelTol',h,'AbsTol',h,'InitialStep',h,'MaxStep',h);
[t,x] = ode45('robot',ts,[x0;xh0],config);

x1 = x(:,1);
x2 = x(:,2);
x3 = x(:,3);
x4 = x(:,4);
xh1 = x(:,5);
xh2 = x(:,6);
xh3 = x(:,7);
xh4 = x(:,8);

[n,m] = size(t);
uu = zeros(n,1);
for k =1:n
    [Xd,u] = robot(t(k),[x1(k);x2(k);x3(k);x4(k);xh1(k);xh2(k);xh3(k);xh4(k)]);
    uu(k) = u;
end

figure, plot(t,x1,t,xh1);
figure, plot(t,x2,t,xh2);
figure, plot(t,x3,t,xh3);
figure, plot(t,x4,t,xh4);
figure, plot(t,uu);




