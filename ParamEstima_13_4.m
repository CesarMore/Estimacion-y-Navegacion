%Filtro de Kalman Extendido Continuo para
% Estimación de Parametros
clc
clear all
close all

T = 0.01;
t = [0:T:20];
x = [0 0 0];
VARwp = 0.1;
VARw = 2.5;   %Puede cambiar
w = sqrt(0.25*(VARw))*((2*randn(length(t),1)-1));
wp = sqrt(0.25*(VARwp))*((2*randn(length(t),1)-1));
VARv = 0.009;
v = sqrt(0.25*(VARv))*((2*randn(length(t),1)-1));
x3real = 4;
% mean(w)
% mean(wp)
% var(w)
% var(wp)
for i=1:length(t)-1
    [tt,xx] = ode45(@Modelo13_4,[t(i) t(i+1)],x(i,:),[],w(i),wp(i));
    x(i+1 ,:) = xx(end,:);
    
end
%-----Filtro de Kalman---------
R = [var(v) 0;0 var(v)];
H = [1 0 0;
    0 1 0];
Q = [var(w) 0;0 var(wp)];
xhat = [mean(x(1,1)) mean(x(1,2)) -1];   %Se puede mover
P = [0 0 0;
    0 0 0;
    0 0 20];
% ex = x(1,:)'-xhat';
% ex1 = (ex)*transpose(ex);
% P = [mean(ex1(1,1)) mean(ex1(1,2)) mean(ex1(1,3));
%     mean(ex1(2,1)) mean(ex1(2,2)) mean(ex1(2,3));
%     mean(ex1(3,1)) mean(ex1(3,2)) 20];
b = -0.4;
y = [x(:,1)+v x(:,2)+v];

for i = 1:length(t)-1
F = [0 1 0;
    xhat(i,3) b xhat(i,1);
    0 0 0];
L = [0 0;
    -xhat(i,3) 0;
    0 1];
% L = [0 0;
%     -xhat(i,3);
%     0];
f0 = [xhat(i,2);
    xhat(i,3)*xhat(i,1) + b*xhat(i,2) - xhat(i,3);
    0];
K = P(:,:,i)*H'*inv(R);
Ppunto = F*P(:,:,i) + P(:,:,i)*F' + L*Q*L' - P(:,:,i)*H'*inv(R)*H*P(:,:,i);
P(:,:,i+1) = P(:,:,i) + T*Ppunto;
xhatpunto = f0 + K*(y(i,:)' - H*xhat(i,:)');
xhat(i+1,:) = xhat(i,:) + T*xhatpunto';

end

errorx3real = x3real - xhat(:,3);
errorx1 = x(:,1) - xhat(:,1);
errorx2 = x(:,2) - xhat(:,2);
errorx3 = x(:,3) - xhat(:,3);

figure(1)
plot(t,errorx3real)
title('error')

figure(2)
subplot(2,2,1); hold on; box on;
plot(t,x(:,1)); 
xlabel('Tiempo (Segundos)'); ylabel('');
legend('x1')
subplot(2,2,2); hold on; box on;
plot(t,x(:,2)); 
xlabel('Tiempo (Segundos)'); ylabel('x'); 
legend('x2')
subplot(2,2,3); hold on; box on;
plot(t,x(:,3)); 
xlabel('Tiempo (Segundos)'); ylabel(''); 
legend('x3')

figure(3)
subplot(2,2,1)
plot(t,xhat(:,1));
xlabel('Tiempo (Segundos)'); ylabel('');
legend('xhat1')
subplot(2,2,2)
plot(t,xhat(:,2));
xlabel('Tiempo (Segundos)'); ylabel(''); 
legend('xhat1')
subplot(2,2,3)
plot(t,xhat(:,3));
xlabel('Tiempo (Segundos)'); ylabel('');
legend('xhat1')

figure(4)
subplot(2,2,1)
plot(t,errorx1);
xlabel('Tiempo (Segundos)'); ylabel(''); 
legend('error1')
subplot(2,2,2)
plot(t,errorx2);
xlabel('Tiempo (Segundos)'); ylabel('xhat1');
legend('error2')
subplot(2,2,3)
plot(t,errorx3);
xlabel('Tiempo (Segundos)'); ylabel('');
legend('error3')

% figure(5)
% plot(t,P(:,:,i+1))
% xlabel('Tiempo (Segundos)'); ylabel('Covarianza');
% title('Coravianza')

