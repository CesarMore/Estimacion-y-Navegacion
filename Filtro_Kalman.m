clc
clear all
close all

T = 1;
t = [0:T:60];
x = [20; 19; 18];
xhat = [19; 17; 15];
P = eye(3);
R = 0.01;
%R = 30^2;
%v = sqrt(R)*randn;
P_step= [];
xhat_1 = [];
xhat_2 = [];
% e_1 = [];
I = eye(3);
a = 0.01;
v = a*randn(length(t),1);
n = 0;
m = 0;
error_med= [];
error_estim = [];
%
tiempo = [];
TIEMPO = [];
%
for i=1:60
    tiempo(i) = t(i)+0.01;
end
for j=1:60
    TIEMPO(end+1) = t(j);
    TIEMPO(end+1) = tiempo(j);
end
%
for k=1:length(t)-1
    H = [1 0 0];
    y = H*x + v(k);
    F = [1 T (T^2)/2; 0 1 T; 0 0 1];
    P = F*P*F';
    P_p(k+m) = round(P(1,1),1);
    m = m+1;
    xhat = F*xhat;
    xhat_1(k) = xhat(1);
    error_med(k) = x(1)-xhat(1);
    K = P*(H)'*inv(H*P*(H)' + R);
    xhat = xhat + K*(y - H*xhat);
    error_estim(k) = x(1)-xhat(1);
    xhat_2(k) = xhat(1);
    e_1(k) = x(1) - xhat(1);
    %e_1(k) = x(1) - xhat(1);
    P = (I - K*H)*P;
%     P_step(k) = (P(1,1));
    P_p(k+n) = round(P(1,1),2);
    n = n+1;
end

figure(1)
plot(P_p)
hold on 
legend('P_step')
ylabel('covarianza')
xlabel('Tiempo [seg]')
figure(2)
plot(xhat_1,xhat_2)
% ylim([-80 80])
hold on 
legend('xhat1')
ylabel('Posición')
xlabel('Tiempo [seg]')
figure(3)
plot(e_1)
% ylim([-80 80])
hold on 
legend('e_1')
ylabel('Error')
xlabel('Tiempo [seg]')


