% clc 
% clear all 
% close all 
% 
% T = 0.001;
% tf = 2;
% t = [0:T:tf];
% dt = 0.0002;
% x = [0; 0; 0; 0]; % Estado inicial
% xhat = x; % Estado inicial estimada
% %Parametros del sistema 
% L = 0.003;
% J = 0.002;
% Ra = 2;
% lambda = 0.1;
% B = 0.001;
% sigma_i = 0.001;
% sigma_m = 0.1;      
% sigma_T = 0.05;       
% q1 = sqrt(0.25*(sigma_i))*((2*randn(length(t),1)-1));
% q2 = sqrt(0.25*(sigma_m))*((2*randn(length(t),1)-1));
% q3 = sqrt(0.25*(sigma_T))*((2*randn(length(t),1)-1));
% ControlNoise = 0.001; % std dev de incertidumbre en las entradas de control(amps)
% MeasNoise = 0.1; % std dev del ruido de medición (amps)
% R = [MeasNoise^2 0; 0 MeasNoise^2]; % Covarianza del ruido de medición
% xdotNoise = [ControlNoise/L; ControlNoise/L; 0.5; 0];
% % Definición de la covarianza Q del ruido del proceso en tiempo continuo
% Q = [xdotNoise(1)^2 0 0 0;
%      0 xdotNoise(2)^2 0 0;
%      0 0 xdotNoise(3)^2 0;
%      0 0 0 xdotNoise(4)^2];
% % ex = x(:,1)'-xhat';
% % ex1 = (ex)*transpose(ex);
% % P = [mean(ex1(1,1)) mean(ex1(2,1)) mean(ex1(3,1)) mean(ex1(4,1));
% %      mean(ex1(1,2)) mean(ex1(2,2)) mean(ex1(3,2)) mean(ex1(4,2));
% %      mean(ex1(1,3)) mean(ex1(2,3)) mean(ex1(3,3)) mean(ex1(4,3));
% %      mean(ex1(1,4)) mean(ex1(2,4)) mean(ex1(3,4)) mean(ex1(4,4))];
% P = 1*eye(4); % Covarianza de estimación de estado inicial
% ua = sin(2*pi*t);
% ub = cos(2*pi*t);
% VARv1 = MeasNoise^2;
% VARv2 = MeasNoise^2;
% v1 = sqrt(0.25*(VARv1))*((2*randn(length(t),1)-1));
% v2 = sqrt(0.25*(VARv2))*((2*randn(length(t),1)-1));
% 
% for i=1:length(t)-1
%     [tt,xx] = ode45(@ModeloEjem13_1,[t(i) t(i+1)],x(i,:),[],ua(i),ub(i),q1(i),q2(i),q3(i));
%     x(i+1,:) = xx(end,:);
%     x1(:,i+1) = x(i+1,1);
%     x2(:,i+1) = x(i+1,2);
%     x3(:,i+1) = x(i+1,3);
%     x4(:,i+1) = x(i+1,4);
% end
% 
% % Calculo de las matrices
% A = [-(Ra/L) 0 (lambda/L)*sin(xhat(4)) (lambda/L)*xhat(3)*cos(xhat(4));
%     0 -(Ra/L) -(lambda/L)*cos(xhat(4)) (lambda/L)*xhat(3)*sin(xhat(4));
%     -((3*lambda)/(2*J))*sin(xhat(4)) ((3*lambda)/(2*J))*cos(xhat(4)) -B/J -((3*lambda)/(2*J))*(xhat(1)*cos(xhat(4))+xhat(2)*sin(xhat(4)));
%     0 0 1 0];
% C = [1 0 0 0; 0 1 0 0];
% % Calculo de la ganancia de Kalman
% K = P*C'*inv(R);
% % Actualización del estado estimado
% 
% for i=1:length(t)-1
%     Ppunto = A*P(:,:,i)+P(:,:,i)*A'+Q-P(:,:,i)*(C)'*inv(R)*C*P(:,:,i);
%     P(:,:,i+1) = P(:,:,i)+(Ppunto)*T;
% %     Pp = A*P(:,i) + P(:,i)*A' + Q - P(:,i)*C'*inv(R)*C*P(:,i);
% %     P(:,i+1) =  P(:,i) + T*(Pp(:,i));
%     K = P(:,:,i+1)*C'*inv(R);
%     f = [-(Ra/L)*xhat(1,i) + (lambda/L)*xhat(3,i)*sin(xhat(4,i)) + (ua(i)+q1(i))/L;
%        -(Ra/L)*xhat(2,i) - (lambda/L)*xhat(3,i)*cos(xhat(4,i)) + (ub(i)+q2(i))/L;
%        -((3*lambda)/(2*J))*xhat(1,i)*sin(xhat(4,i)) + ((3*lambda)/(2*J))*xhat(2,i)*cos(xhat(4,i)) - B/J*xhat(3,i)+q3(i);
%        xhat(3,i)];    
%     y = [x(1,i) + v1(i);
%          x(2,i) + v2(i)];
%     h = [xhat(1,i);
%          xhat(2,i)];
%     xdot = f + K*(y-h);
%     xhat(:,i+1) = xhat(:,i)'+T*xdot;
%     % Mantener el ángulo estimado entre 0 y 2*pi
% %     xhat(4) = mod(xhat(4), 2*pi);
% 
%     xhat1(:,i+1) = xhat(i+1,1);
%     xhat2(:,i+1) = xhat(i+1,2);
%     xhat3(:,i+1) = xhat(i+1,3);
%     xhat4(:,i+1) = xhat(i+1,4);
% end
% 
% figure(1)
% subplot(2,2,1); hold on; box on;
% plot(t,x1,t,xhat1); 
% xlabel('Tiempo (Segundos)'); ylabel('Corriente A (Amps)');
% legend('x1','x1hat')
% 
% subplot(2,2,2); hold on; box on;
% plot(t,x2,t,xhat2); 
% xlabel('Tiempo (Segundos)'); ylabel('Corriente B (Amps)'); 
% legend('x2','x2hat')
% 
% subplot(2,2,3); hold on; box on;
% plot(t,x3,t,xhat3); 
% xlabel('Tiempo (Segundos)'); ylabel('Aceleración (rad/sec)'); 
% legend('x3','x3hat')
% 
% subplot(2,2,4); hold on; box on;
% plot(t,x4,t,xhat4); 
% xlabel('Tiempo (Segundos)'); ylabel('Posición (rad)'); 
% legend('x4','x4hat')

%%
clc 
clear all 
close all 

T = 0.001;
tf = 2;
t = [0:T:tf];
x = [0.1 0.1 0 0]; % Estado inicial
xhat = x; % Estado inicial estimada
%Parametros del sistema 
L = 0.003;
J = 0.002;
Ra = 2;
lambda = 0.1;
B = 0.001;
sigma_i = 0.001;
sigma_m = 0.1;      
sigma_T = 0.05; 
VAR  = sigma_i^2;
VAR1 = sigma_m^2;
VAR2 = sigma_T^2;
q1 = sqrt(0.25*(VAR))*((2*randn(length(t),1)-1));
q2 = sqrt(0.25*(VAR1))*((2*randn(length(t),1)-1));
q3 = sqrt(0.25*(VAR2))*((2*randn(length(t),1)-1));
ControlNoise = 0.001; % std dev de incertidumbre en las entradas de control(amps)
MeasNoise = 0.1; % std dev del ruido de medición (amps)
R = [MeasNoise^2 0; 0 MeasNoise^2]; % Covarianza del ruido de medición
xdotNoise = [ControlNoise/L; ControlNoise/L; 0.5; 0];
% Definición de la covarianza Q del ruido del proceso en tiempo continuo
Q = [xdotNoise(1)^2 0 0 0;
     0 xdotNoise(2)^2 0 0;
     0 0 xdotNoise(3)^2 0;
     0 0 0 xdotNoise(4)^2];
ex = x(1,:)'-xhat';
ex1 = (ex)*transpose(ex);
P = [mean(ex1(1,1)) mean(ex1(1,2)) mean(ex1(1,3)) mean(ex1(1,4));
     mean(ex1(2,1)) mean(ex1(2,2)) mean(ex1(2,3)) mean(ex1(2,4));
     mean(ex1(3,1)) mean(ex1(3,2)) mean(ex1(3,3)) mean(ex1(3,4));
     mean(ex1(4,1)) mean(ex1(4,2)) mean(ex1(4,3)) mean(ex1(4,4))];
%P = eye(4); % Covarianza de estimación de estado inicial
ua = sin(2*pi*t);
ub = cos(2*pi*t);

VARv1 = MeasNoise^2;
VARv2 = MeasNoise^2;
v1 = sqrt(0.25*(VARv1))*((2*randn(length(t),1)-1));
v2 = sqrt(0.25*(VARv2))*((2*randn(length(t),1)-1));

for i=1:length(t)-1
    [tt,xx] = ode45(@ModeloEjem13_1,[t(i) t(i+1)],x(i,:),[],ua(i),ub(i),q1(i),q2(i),q3(i));
    x(i+1,:) = xx(end,:);
    x1(:,i+1) = x(i+1,1);
    x2(:,i+1) = x(i+1,2);
    x3(:,i+1) = x(i+1,3);
    x4(:,i+1) = x(i+1,4);
end

% Calculo de las matrices
A = [-(Ra/L) 0 (lambda/L)*sin(xhat(4)) (lambda/L)*xhat(3)*cos(xhat(4));
    0 -(Ra/L) -(lambda/L)*cos(xhat(4)) (lambda/L)*xhat(3)*sin(xhat(4));
    -((3*lambda)/(2*J))*sin(xhat(4)) ((3*lambda)/(2*J))*cos(xhat(4)) -B/J -((3*lambda)/(2*J))*(xhat(1)*cos(xhat(4))+xhat(2)*sin(xhat(4)));
    0 0 1 0];
C = [1 0 0 0; 0 1 0 0];
% Calculo de la ganancia de Kalman
%K = P*C'*inv(R);
% Actualización del estado estimado

for i=1:length(t)-1
    Ppunto = A*P(:,:,i)+P(:,:,i)*A'+Q-P(:,:,i)*(C)'*inv(R)*C*P(:,:,i);
    P(:,:,i+1) = P(:,:,i)+(Ppunto)*T;
%     Pp = A*P(:,i) + P(:,i)*A' + Q - P(:,i)*C'*inv(R)*C*P(:,i);
%     P(:,i+1) =  P(:,i) + T*(Pp(:,i));
    K = P(:,:,i+1)*C'*inv(R);
    f = [-(Ra/L)*xhat(i,1) + (lambda/L)*xhat(i,3)*sin(xhat(i,4)) + (ua(i)+q1(i))/L;
       -(Ra/L)*xhat(i,2) - (lambda/L)*xhat(i,3)*cos(xhat(i,4)) + (ub(i)+q2(i))/L;
       -((3*lambda)/(2*J))*xhat(i,1)*sin(xhat(i,4)) + ((3*lambda)/(2*J))*xhat(i,2)*cos(xhat(i,4)) - B/J*xhat(i,3)+q3(i);
       xhat(i,3)];    
    y = [x(i,1) + v1(i);
         x(i,2) + v2(i)];
    h = [xhat(i,1);
         xhat(i,2)];
    xhatdot = f + K*(y-h);
    xhat(i+1,:) = xhat(i,:)'+T*xhatdot;
    % Mantener el ángulo estimado entre 0 y 2*pi
%     xhat(4) = mod(xhat(4), 2*pi);

    xhat1(:,i+1) = xhat(i+1,1);
    xhat2(:,i+1) = xhat(i+1,2);
    xhat3(:,i+1) = xhat(i+1,3);
    xhat4(:,i+1) = xhat(i+1,4);
end

figure(1)
subplot(2,2,1); hold on; box on;
plot(t,x1,t,xhat1); 
xlabel('Tiempo (Segundos)'); ylabel('Corriente A (Amps)');
legend('x1','x1hat')

subplot(2,2,2); hold on; box on;
plot(t,x2,t,xhat2); 
xlabel('Tiempo (Segundos)'); ylabel('Corriente B (Amps)'); 
legend('x2','x2hat')

subplot(2,2,3); hold on; box on;
plot(t,x3,t,xhat3); 
xlabel('Tiempo (Segundos)'); ylabel('Aceleración (rad/sec)'); 
legend('x3','x3hat')

subplot(2,2,4); hold on; box on;
plot(t,x4,t,xhat4); 
xlabel('Tiempo (Segundos)'); ylabel('Posición (rad)'); 
legend('x4','x4hat')
