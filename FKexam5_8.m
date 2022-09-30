clc
clear all
close all

T = 5;
t = [0:T:150];
% t =[0:150];
VAR = 1000000;
VAR1 = 0.000289;
VAR2 = 1178.78;
VAR3 = 0.000000013;
v1=sqrt(0.25*(VAR))*((2*randn(length(t),1)-100));
mean(v1)
var(v1)
v2=sqrt(0.25*(VAR1))*((2*randn(length(t),1)-1));
mean(v2)
var(v2)
w1=sqrt(0.25*(VAR2))*((2*randn(length(t),1)-0.0099999));
mean(w1)
var(w1)
w2=sqrt(0.25*(VAR3))*((2*randn(length(t),1)-0.0099999));
mean(w2)
var(w2)
v = [v1;v2];
sigmarsqrd = 1000000;
sigmathetasqrd = 0.017^2;
sigma1sqrd = (103/3)^2;
sigma2sqrd = 0.000000013;
rho = 0.5;
F = [1 T 0 0 0 0;
    0 1 1 0 0 0;
    0 0 rho 0 0 0;
    0 0 0 1 T 0;
    0 0 0 0 1 1;
    0 0 0 0 0 rho];
Q = [0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 sigma1sqrd 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 sigma2sqrd];
x = [0 0 0 0 0 0];
H = [1 0 0 0 0 0;
    0 0 0 1 0 0];
w = [0;0;w1;0;0;w2];
%y = H*x';
R=VAR;
xmas = x;
% pmas = 100*eye(6);
Pmas = [sigmarsqrd sigmarsqrd/T 0 0 0 0;
        sigmarsqrd/T 2*sigmarsqrd/T^2+sigma1sqrd 0 0 0 0;
        0 0 sigma1sqrd 0 0 0;
        0 0 0 sigmathetasqrd sigmathetasqrd/T 0;
        0 0 0 sigmathetasqrd/T 2*sigmathetasqrd/T^2+sigma2sqrd 0;
        0 0 0 0 0 sigma2sqrd];

for i=2:length (t)
    x(i,:) = F*x(i-1,:)' + w(i);
    y(i,:) = H*x(i,:)' + v(i);
end

for i=2:length(t)
    %Prediccion
    xmenos(i,:) = F*xmas(i-1,:)';
    Pmenos = F*Pmas*F'+ Q;
    Pminus(:,:,i)=Pmenos;
    %Actualizacion
    K = Pmenos*H'*inv(H*Pmenos*H' +R);
    Pmas = Pmenos - K*H*Pmenos;
    Pplus(:,:,i) = Pmas;
    xmas(i,:) = xmenos(i,:)'+ K*(y(i)-H*xmenos(i,:)');
end

figure(1)
plot(t,y(:,1),t,xmas(:,1),t,x(:,1))
legend('y1','xposteriori','xactual')
% ylim ([-60000 500])

figure(2)
subplot(211)
plot(t,xmas(:,2),t,x(:,2))
legend('x2posteriori','x2actual')
subplot(212)
plot(t,xmas(:,3),t,x(:,3))
legend('x3posteriori','x3actual')

figure(3)
plot(t,x(:,1)-y',t,x(:,1)-xmas(:,1))
legend('E_{medicion}','E_{estimado}')
%plot(t,x(:,1)-y',t,x(:,1)-xmas(:,1))

figure(4)
hold on
for i=1:50
    plot([i i+1],[Pplus(1,1,i) Pminus(1,1,i+1)],'b','LineWidth',2)
    plot([i+1 i+1],[Pminus(1,1,i+1) Pplus(1,1,i+1)],'r','LineWidth',2)
end

