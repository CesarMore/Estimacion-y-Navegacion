clc
clear all
close all

ti = 0;
T = 0.01; 
tf = 10;

t = (ti:T:tf)';

u = [1;1];
d = [0;0];
x = [0 0 0];

x1 = [0];
x2 = [0];
x3 = [0];

for i = 1 : length(t)-1
    
    if t(i) > 1 
        d(1) = 2;
    end

    if t(i) > 1
        d(2) = 3;
    end

    [tt, xx] = ode45(@modelsystem2, [t(i) t(i+1)], x(i,:),[],d);
    x(i+1,:)   = xx(end,:);
    
    x1(:,i+1) = x(i+1,1);
    x2(:,i+1) = x(i+1,2);
    x3(:,i+1) = x(i+1,3);
    
end

figure(1)
plot(t,x1,t,x2,t,x3)
grid on
xlabel("Tiempo, [Seg]")
ylabel("Salida, y")


%%
%%%%%%%%%%%%%%%

A   = [-0.8 1 1.6; 0 -3 2; 0 0 -6];
Bu  = [0 0 ; 1 0; 0 1];
Bd  = [0.8 0; 0 -1; -0.4 1.2];
L   = [40 0 -20; 0 -20 60];
LBd = [-40 24; 24 -92];

x = [x1 x2 x2];
z = [0,0];
zp1 = -LBd*(z(:,1) + L*x(:,1)) - L*(A*x(:,1)+Bu*u);



