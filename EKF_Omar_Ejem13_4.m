%Filtro de Kalman Extendido Continuo para Estimacion de Parametros

clc
clear all
close all
T = 0.01;
t = [0:T:10];
%----Ruidos--------
q1ds = 0.01;
VARw = q1ds^2;
VARwp = 0.1;
VARv = 0.1;
w = sqrt(0.25*(VARw))*((2*randn(length(t),1)-1));
wp = sqrt(0.25*(VARwp))*((2*randn(length(t),1)-1));
v = sqrt(0.25*(VARv))*((2*randn(length(t),1)-1));
% mean(w)
% mean(wp)
% 
% var(w)
% var(wp)
%--------------------
x = [0.1 2 -8];
y = [0;0];
b = -0.4;
for i=1:length(t)-1
    [tt,xx] = ode45(@modelox,[t(i) t(i+1)],x(i,:),[],w(i),wp(i));
    x(i+1 ,:) = xx(end,:);
    
end
% figure(1)
% plot(t,x(:,1));
% figure(2)
% plot(t,x(:,2));
% figure(3)
% plot(t,x(:,3));
%----------------------------------------------
%Filtro de Kalman Extendido Continuo
%-----------------------------------------------
xhat = [mean(x(1,1)) mean(x(1,2)) mean(x(1,3))];
C = [1 0 0;
    0 1 0];
M = [1;
    1];
Q = [var(w) 0;
    0 var(wp)];
Rx = [var(v)];
Rtilde = M*Rx*transpose(M);
P = [0 0 0;0 0 0;0 0 20];
for i=1:length(t)-1
   A = [0 1 0;
       xhat(i,3) b xhat(i,1);
       0 0 0];
   Lx = [0 0;-xhat(i,3) 0;0 1];
   Qtilde = Lx*Q*transpose(Lx);
   Ppunto = A*P(:,:,i)+P(:,:,i)*A'+Qtilde-P(:,:,i)*transpose(C)*inv(Rtilde)*C*P(:,:,i);
   P(:,:,i+1) = P(:,:,i)+(Ppunto)*T;
   K = P(:,:,i+1)*C'*inv(Rtilde);
   f = [xhat(i,2);
       xhat(i,3)*xhat(i,1) + b*xhat(i,2) - xhat(i,3);
       0];
   y(:,i+1) = [x(i,1) + v(i);
       x(i,2) + v(i)];
   h = [xhat(i,1);
       xhat(i,2)];
   funcion = f + K*(y(:,i)-h);
   xhat(i+1,:) = xhat(i,:)'+T*funcion;
end



function dxdt = modelox(t,x,w,wp)
b = -0.4;
dxdt(1,1)=x(2);
dxdt(2,1)=x(3)*x(1) + b*x(2) - x(3)*w;
dxdt(3,1)=wp;
end