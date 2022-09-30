clc
clear all
close all

T = 1;
t = [0:T:100];
x = [10; 5];
xhat = [8; 7];
P = eye(2);
R = 0.01;
%v = sqrt(R)*randn;
P_paso1 = [];
P_paso2 = [];
xhat_v1 = [];
xhat_v2 = [];
I = eye(2);
a = 0.1;
v = a*randn(length(t),1);

for k=1:length(t)-1
    H = [1 0.99^(k-1)];
    y = H*x + v(k);
    K = P*(H)'*inv(H*P*(H)' + R);
    xhat = xhat + K*(y - H*xhat);
    xhat_v1(k) = xhat(1);
    xhat_v2(k) = xhat(2);
    P = (I - K*H)*P;
    P_paso1(k) = P(1,1);
    P_paso2(k) = P(2,2);
end

figure(1)
plot(P_paso1)
ylim([0 1])
hold on 
grid on 
plot(P_paso2)
legend('P(1,1)', 'P(2,2)')
ylabel('Varianza')
xlabel('Tiempo [seg]')
figure(2)
plot(xhat_v1)
hold on 
grid on
plot(xhat_v2)
legend('xhat1', 'xhat2')
ylabel('Valores estimados de X')
xlabel('Tiempo [seg]')

