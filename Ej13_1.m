%Filtro de Kalman Extendido Continuo

clc
clear all
close all

global R L lambda ua ub q1 q2 q3 F J
R = 0.5;
L = 0.5;
lambda = 0.05;
F = 0.0015;
J = 0.001;

T = 0.01;
t = [0:T:2];

ua = sin(2*pi*t);
ub = cos(2*pi*t);

%----Ruidos--------

q1ds = 0.01;
q2ds = 0.01;
q3ds = 0.5;

VAR = q1ds^2
VAR1 = q2ds^2
VAR2 = q3ds^2
q1 = sqrt(0.25*(VAR))*((2*randn(length(t),1)-1));
q2 = sqrt(0.25*(VAR1))*((2*randn(length(t),1)-1));
q3 = sqrt(0.25*(VAR2))*((2*randn(length(t),1)-1));
% mean(q1)
% mean(q2)
% mean(q3)
% var(q1)
% var(q2)
% var(q3)

v1ds = 0.1;
v2ds = 0.1;

VARv1 = 0.1^2;
VARv2 = 0.1^2;

v1 = sqrt(0.25*(VARv1))*((2*randn(length(t),1)-1));
v2 = sqrt(0.25*(VARv2))*((2*randn(length(t),1)-1));
% 
% mean(v1)
% mean(v2)
% var(v1)
% var(v2)

%--------------------

x = [0 0 0.1 10];


for i=1:length(t)-1
    [tt,xx] = ode45(@modelox,[t(i) t(i+1)],x(i,:),[],ua(i),ub(i),q1(i),q2(i),q3(i));
    x(i+1 ,:) = xx(end,:);
    
end


figure(1)
plot(t,x(:,1));
figure(2)
plot(t,x(:,2));
figure(3)
plot(t,x(:,3));
figure(4)
plot(t,x(:,4));



%----------------------------------------------
%Filtro de Kalman Extendido Continuo
%-----------------------------------------------

%xhat = [mean(x(:,1)) mean(x(:,2)) mean(x(:,3)) mean(x(:,4))];
xhat = [mean(x(1,1)) mean(x(1,2)) mean(x(1,3)) mean(x(1,4))];

%si o si debe ser de A =  4x4 L = 4x4 

A = [-(R/L) 0 (lambda/L)*sin(xhat(4)) (lambda/L)*xhat(3)*cos(xhat(4));
    0 -(R/L) -(lambda/L)*cos(xhat(4)) (lambda/L)*xhat(3)*sin(xhat(4));
    -((3*lambda)/(2*J))*sin(xhat(4)) ((3*lambda)/(2*J))*cos(xhat(4)) -F/J -((3*lambda)/(2*J))*(xhat(1)*cos(xhat(4))+xhat(2)*sin(xhat(4)));
    0 0 1 0];

Lx = [1/L 0 0;
    0 1/L 0;
    0 0 1;
    0 0 0];

C = [1 0 0 0;
    0 1 0 0];

M = [1 0;
    0 1];


Q = [var(q1) 0 0;
    0 var(q2) 0;
    0 0 var(q3)];

Rx = [var(v1) 0;0 var(v2)];


Qtilde = Lx*Q*transpose(Lx);
Rtilde = M*Rx*transpose(M);

ex = x(1,:)'-xhat';
ex1 = (ex)*transpose(ex);
P = [mean(ex1(1,1)) mean(ex1(1,2)) mean(ex1(1,3)) mean(ex1(1,4));
    mean(ex1(2,1)) mean(ex1(2,2)) mean(ex1(2,3)) mean(ex1(2,4));
    mean(ex1(3,1)) mean(ex1(3,2)) mean(ex1(3,3)) mean(ex1(3,4));
    mean(ex1(4,1)) mean(ex1(4,2)) mean(ex1(4,3)) mean(ex1(4,4))];


for i=1:length(t)-1
    
   Ppunto = A*P(:,:,i)+P(:,:,i)*A'+Qtilde-P(:,:,i)*transpose(C)*inv(Rtilde)*C*P(:,:,i);

   P(:,:,i+1) = P(:,:,i)+(Ppunto)*T;

   K = P(:,:,i+1)*C'*inv(Rtilde);

   f = [-(R/L)*xhat(i,1) + (lambda/L)*xhat(i,3)*sin(xhat(i,4)) + (ua(i)+q1(i))/L;
       -(R/L)*xhat(i,2) - (lambda/L)*xhat(i,3)*cos(xhat(i,4)) + (ub(i)+q2(i))/L;
       -((3*lambda)/(2*J))*xhat(i,1)*sin(xhat(i,4)) + ((3*lambda)/(2*J))*xhat(i,2)*cos(xhat(i,4)) - F/J*xhat(i,3)+q3(i);
       xhat(i,3)];

   y = [x(i,1) + v1(i);
       x(i,2) + v2(i)];

   h = [xhat(i,1);
       xhat(i,2)];

   funcion = f + K*(y-h);

   xhat(i+1,:) = xhat(i,:)'+T*funcion;


end



figure(5)
plot(t,xhat(:,1));
figure(6)
plot(t,xhat(:,2));
figure(7)
plot(t,xhat(:,3));
figure(8)
plot(t,xhat(:,4));





figure(9)
plot(t,x(:,1),t,xhat(:,1));
legend('x1 real','x1 filtrada')
title('x1 vs xh 1')
set(gca,'FontSize',13)
xlabel('Tiempo [s]')
ylabel('ia')

figure(10)
plot(t,x(:,2),t,xhat(:,2));
legend('x2 real','x2 filtrada')
title('x2 vs xh 2')
set(gca,'FontSize',13)
xlabel('Tiempo [s]')
ylabel('ib')

figure(11)
plot(t,x(:,3),t,xhat(:,3));
legend('x3 real','x3 filtrada')
title('x3 vs xh 3')
set(gca,'FontSize',13)
xlabel('Tiempo [s]')
ylabel('Omega')

figure(12)
plot(t,x(:,4),t,xhat(:,4));
legend('x4 real','x4 filtrada')
title('x4 vs xh 4')
set(gca,'FontSize',13)
xlabel('Tiempo [s]')
ylabel('Theta')





function dxdt = modelox(t,x,ua,ub,q1,q2,q3)

global R L lambda F J

dxdt(1,1)=-(R/L)*x(1) + (lambda/L)*x(3)*sin(x(4)) + (ua+q1)/L;
dxdt(2,1)=-(R/L)*x(2) - (lambda/L)*x(3)*cos(x(4)) + (ub+q2)/L;
dxdt(3,1)=-((3*lambda)/(2*J))*x(1)*sin(x(4)) + ((3*lambda)/(2*J))*x(2)*cos(x(4)) - F/J*x(3)+q3;
dxdt(4,1)=x(3);
end














