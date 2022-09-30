clc
clear all 
close all

T =0.01;
t = [0:T:11];


% f1 = 5*t + sin(t) + 0.01*cos(10*t) + 0.01*rand(1,length(t));
f1 = 0.5*sin(0.5*t) + 0.5*cos(t);% + 0.01*rand(1,length(t));
f1p = 0.25*cos(0.5*t) - 0.5*sin(t);
f2p = -0.125*sin(0.5*t) - 0.5*cos(t);

x1 = 0;
x2 = 0;
x3 = 0;
% alpha = 40;  %k2
% lamda = 20;  %k1
L = 10;

for i=1:length(t)-1

    if i==1
    f1p_euler = 0;
    else
    f1p_euler(i+1) = (f1(i) - f1(i-1))/T;
    end

    v0 = (-3*L^(1/3)*abs(x1(i) - f1(i))^(2/3)*sign(x1(i) - f1(i))+x2(i));

    x1(i+1) = x1(i) + T*v0;

    v1 = (-1.5*L^(1/2)*abs(x2(i) - v0)^(1/2)*sign(x2(i) - v0)+x3(i));

    x2(i+1) = x2(i) + T*v1;

    x3(i+1) = x3(i) + T*(-1.1*L*sign(x3(i) - v1));

    f1p_robusto(i+1) = v0;
    f2p_robusto(i+1) = v1;

end

%Filtro de primer orden 
tf1 = tf([1 0], conv([0.05 1], [0.05 1]));
f1p_tf = lsim(tf1,f1,t);

tf2 = tf([1], [0.05 1]);
f1p_rfilt = lsim(tf2,f1p_robusto,t);

tf3 = tf([1], [0.05 1]);
f2p_rfilt = lsim(tf3,f2p_robusto,t);

% plot(t,f1p,t,f1p_euler);
% legend('f1p', 'f1p_euler')
plot(t,f2p,t,f2p_robusto,t,f2p_rfilt);
legend('f2p', 'f2p_robusto','f2p_rfilt');
