
clc
clear all
close all

T = 5;
tf = 60;
t = [0:T:tf];
% tiempo = []
% TIEMPO = []
% for i=1:tf
%     tiempo(i) = t(i) + 0.01;
% end
% for j=1:tf
%     TIEMPO(end+1) = t(j);
%     TIEMPO(end+1) = tiempo(j);
% end

x = [35;7;5];
%x = 200*x;
xhat = [40; 2; 2.5];
%xhat = 200*xhat;
P = eye(3);
P = 500*P;
R = 30^2;
%R = 0.02;
xhat_1 = [];
xhat_2 = [];
P_p = [];
P_p2 = [];
m = 0;
n = 1;
error_med= [];
error_estim = [];

for k=1:length(t)-1
H = [1 0 0];
y = H*x;
y = awgn(y,0.5,'measured');
F = [1 T (T^2)/2;0 1 T;0 0 1];
P = F*P*transpose(F);
P_p(k+m) = P(1,1);
m = m+1;
xhat = F*xhat;
xhat_1(k) = xhat(1);
error_estim(k) = x(1)-xhat(1);
K = P*transpose(H)*inv(H*P*transpose(H) + R);
xhat = xhat + K*(y - H*xhat);
xhat_2(k) = xhat(1);
error_med(k) = x(1)-xhat(1);
P = P - K*H*P;
P_p(k+n) = P(1,1);
n = n+1;
%xhat = xhatpriori;
end

figure(1)
plot(t,P_p,'LineWidth',2,'Color',[0 1 0])
grid on
figure(2)
grid on
plot(error_med)
hold on
plot(error_estim)
figure(3)
grid on
% plot(xhat_1)
% % hold on
plot(xhat_2)
figure(4)
plot(t,y,t,xhat_1)









