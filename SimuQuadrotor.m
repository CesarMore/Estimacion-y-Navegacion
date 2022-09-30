clc 
clear all 
close all 

ti = 0;
T = 0.01;
tf = 10;
t = [ti:T:tf];

x = [0 0 0 0 0 0];

Kp = [6 5 7]';
Kd = [4 4.2 5]';

r = 0.2*sin(t);
rp = 0.2*cos(t);
L = 0.5;
tau1 = 0;
tau2 = 0;
tau3 = 0;

ex = 0;
exp = 0;
ey = 0;
eyp = 0;
ez = 0;
ezp = 0;
taux = 0;
tauy = 0;
tauz = 0;
Jx=0.0552;
Jy=0.0552;
Jz=0.1104;
Jm = 0.56;
Kafx = 0.008;
Kafy = 0.008;
Kafz = 0.009;
Omega1 = 0;
Omega2 = 0;
Omega3 = 0;
Omega4 = 0;
Omega = Omega1-Omega2+Omega3+Omega4; 
Kt = 1;
Kl = 1;

for i=1:length(t)-1
    [tt,xx] = ode45(@modeloQuadrotor,[t(i) t(i+1)],x(i,:),[],tau1(i),tau2(i),tau3(i),Omega1(i),Omega2(i),Omega3(i),Omega4(i));
    x(i+1,:) = xx(end,:);
    
    %%error para roll
    ex(i+1) = x(i+1,1) - r(i); 
    exp(i+1) = x(i+1,2) - rp(i);
    %%error para roll
    ey(i+1) = x(i+1,3) - r(i); 
    eyp(i+1) = x(i+1,4) - rp(i);
    %%error para roll
    ez(i+1) = x(i+1,5) - r(i); 
    ezp(i+1) = x(i+1,6) - rp(i);
    
    %%control tau
    taux(i+1) = (Kafx*x(i+1,1) - (Jy-Jz)*x(i+1,4)*x(i+1,6) - Jm*Omega) + Kp(1)*ex(i+1)+Kd(1)*exp(i+1);
    tauy(i+1) = (Kafy*x(i+1,4) - (Jz-Jx)*x(i+1,2)*x(i+1,6) - Jm*Omega) + Kp(2)*ey(i+1)+Kd(2)*eyp(i+1);
    tauz(i+1) = (Kafz*x(i+1,6) - (Jx-Jy)*x(i+1,2)*x(i+1,4)) + Kp(3)*ez(i+1)+Kd(3)*ezp(i+1);
    %%Omegas
    Omega4(i+1) = sqrt(Omega2(i)^2+taux(i+1)/L*Kt); 
    Omega1(i+1) = sqrt(Omega3(i)^2+tauy(i+1)/L*Kt); 
    Omega3(i+1) = sqrt(Omega4(i)^2+Omega2(i)^2-Omega1(i)^2+tauz(i+1)/Kl);
    Omega2(i+1) = sqrt(Omega4(i+1)^2-taux(i+1)/L*Kt);
    
    Omega(i+1) = Omega1(i+1)-Omega2(i+1)+Omega3(i+1)-Omega4(i+1);
end
