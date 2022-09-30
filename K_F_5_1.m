clc
clear all
close all

T = 1;
t = [0:T:60];
x = [22; 19; 18];
xhat = [19; 17; 15];
a = 0.1;
v = a*randn(length(t),1);
P = eye(3);
R = 30^2;
%R = 0.01;
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
y = H*x + v(k);
% y = awgn(y,50,'measured');
F = [1 T (T^2)/2;0 1 T;0 0 1];
P = F*P*transpose(F);
P_p(k+m) = round(P(1,1),1);
m = m+1;
xhat = F*xhat;
xhat_1(k) = xhat(1);
error_med(k) = x(1)-xhat(1);
K = P*transpose(H)*inv(H*P*transpose(H) + R);
xhat = xhat + K*(y - H*xhat);
xhat_2(k) = xhat(1);
error_estim(k) = x(1)-xhat(1);
P = P - K*H*P;
P_p(k+n) = round(P(1,1),2);
n = n+1;
%xhat = xhatpriori;
end

figure(1)
plot(P_p,'LineWidth',2,'Color',[0 1 0])
grid on
figure(2)
grid on
plot(error_med)
hold on
plot(error_estim)
figure(3)
grid on
plot(xhat_1)
hold on
plot(xhat_2)