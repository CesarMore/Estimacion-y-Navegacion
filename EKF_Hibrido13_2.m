clc
clear all
close all

T = 0.004;
t = [0:T:30];
x = [100000 -6100 1/2500];
%x = [100000 -6100 1/2500];
xhatminus = [100000; -6100; 1/2500];
VAR  = 0.01^2;
VAR1 = 0.01^2;
VAR2 = 0.01^2;
w1 = sqrt(0.25*(VAR))*((2*randn(length(t),1)-1));
w2 = sqrt(0.25*(VAR1))*((2*randn(length(t),1)-1));
w3 = sqrt(0.25*(VAR2))*((2*randn(length(t),1)-1));
% mean(w1)
% mean(w2)
% mean(w3)
VARv =0.01^2;
v = sqrt(0.5*(VARv))*((2*randn(length(t),1)-1));

% mean(v)
po = 0.0034;
g = 9.82;
k = 6705.4;
for i=1:length(t)-1
    [tt,xx] = ode45(@ModeloEjem13_2,[t(i) t(i+1)],x(i,:),[],w1(i),w2(i),w3(i));
    x(i+1,:) = xx(end,:);
    x1(:,i+1) = x(i+1,1);
    x2(:,i+1) = x(i+1,2);
    x3(:,i+1) = x(i+1,3);
  
end

figure(1)
subplot(2,2,1); hold on; box on;
plot(t,x1); 
xlabel('Tiempo (Segundos)'); ylabel('p33');

subplot(2,2,2); hold on; box on;
plot(t,x2); 
xlabel('Tiempo (Segundos)'); ylabel('xhat3'); 
legend('wn^2')

subplot(2,2,3); hold on; box on;
plot(t,x3); 
xlabel('Tiempo (Segundos)'); ylabel('xhat1'); 
legend('xhat1')


% C = [1 0 0];
% Pmas = [500 0 0; 0 2000 0; 0 0 1/250000];
% xhatmas = [100010; -6100; 1/2500];
% po = 0.0034;
% g = 32.2;
% k = 2200;
% Q = [VAR 0 0;
%       0     VAR 0;
%       0      0    VAR];
% L = [1 0 0; 0 1 0; 0 0 1];
% H = [1 0 0];
% R = [cov(v)];
% for i=1:length(t)-1
% %     Ppunto(:,:,i) = A*Pmas(:,:,i)+Pmas(:,:,i)*A'+L*Q*L';
% %     Pmas(:,:,i+1) = Pmas(:,:,i)+(Ppunto(:,:,i))*T;
% %     K = Pmas(:,:,i+1)*H'*inv(R);
% %     f = [xhat(2,i) + w1(i); po*(exp(-xhat(1,i)/k))*((xhat(2,i)^2)/2*xhat(3,i)) - g + w2(i); w3(i)];
% %     xhatdot = f;
% %     xhat(:,i+1) = xhat(:,i) + T*xhatdot;
%     
%     A = [0 1 0;
%     -po*exp(-xhatminus(i,1)/k)*(xhatminus(i,2)^2)*(1/2*k*xhatminus(i,3)) po*exp(-xhatminus(i,1)/k)*(xhatminus(i,2))*(1/xhatminus(i,3)) -po*exp(-xhatminus(i,1)/k)*(xhatminus(i,2)^2)*(1/2*xhatminus(i,3)^2);
%     0 0 0];
%     f = [xhatminus(i,2);
%     po*exp(-xhatminus(i,1)/k)*(xhatminus(i,2)^2)*(1/2*xhatminus(i,3)) - g;
%     0];
%     xhatpunto = f;
%     %Integrar xhpunto
%     xhatminus(i+1,:) = xhatminus(i,:)'+T*xhatpunto;
%     Ppunto = A*Pmas+Pmas*(A)' + L*Q*(L)';
%     %Integrar Ppunto
%     Pminus(:,:,i+1) = Pminus(:,:,i)+(Ppunto)*T;
%     xhaxhat1(i+1,:) = xhat(1,i+1);
%     xhat2(i+1,:) = xhat(2,i+1);
%     xhat3(i+1,:) = xhat(3,i+1);
%         xhat2(i+1,:) = xhatminus(2,i+1);
%     xhat3(i+1,:) = xhatminus(3,i+1);
% end 
% 
% F = [0 1 0; 
%         po*(exp(-xhat(1,i)/k))*((xhat(2,i)^2)/(2*k*xhat(3,i))), po*(exp(-xhat(1,i)/k))*((xhat(2,i))/xhat(3,i)), po*(exp(-xhat(1,i)/k))*((xhat(2,i)^2)/(2*xhat(3,i)^2));
%          0 0 0];
% T2 = 0.5;
% td = [0:T2:10];
% xmas = x;
% 
% for i=2:length(td)
%     x(i,:) = F*x(i-1,:)';
%     y(i) = H*x(i,:)' + v(i);
%     %Prediccion
%     xmenos(i,:) = F*xmas(i-1,:)';
%     Pmenos = F*Pmas*F';
%     %Pminus(:,:,i) = Pmenos;
%     %Actualicasion
%     K = Pmenos*H'*inv(H*Pmenos*H' +R);
%     Pmas = (I-K*H)*Pmenos*(I-K*H)' + K*M*R*M'*K';
%     %Pplus(:,:,i) = Pmas;
%     xmas(i,:) = xmenos(i,:) + K'*(y(i) - H*xmenos(i,:)');
% end