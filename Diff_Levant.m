
clc
clear all 
close all

T =0.01;
t = [0:T:11];

% f1 = 5*t + sin(t)  + 0.01*rand(1,length(t));
% f1p = 5 + cos(t);

f1 = 5*t + sin(t) + 0.01*cos(10*t) + 0.01*rand(1,length(t));
f1p = 5 + cos(t) -0.1*sin(10*t);
x1 = 0;
x2 = 0;
alpha = 40;  %k2
lamda = 20;  %k1

for i=1:length(t)-1
    if i==1
    f1p_euler = 0;
    else
    f1p_euler(i+1) = (f1(i) - f1(i-1))/T;
    end

    u = (-lamda*abs(x1(i) - f1(i))^(1/2)*sign(x1(i) - f1(i)) + x2(i));
    
    x1(i+1) = x1(i) + T*u;
    x2(i+1) = x2(i) + T*(-alpha*sign(x1(i) - f1(i)));

    f1p_robusto(i+1) = u;
end

%Filtro de primer orden 
tf1 = tf([1 0], conv([0.05 1], [0.05 1]));
f1p_tf = lsim(tf1,f1,t);

tf2 = tf([1], [0.05 1]);
f1p_rfilt = lsim(tf2,f1p_robusto,t);

% plot(t,f1p,t,f1p_euler);
% legend('f1p', 'f1p_euler')
plot(t,f1p, t, f1p_tf,t, f1p_rfilt);
legend('f1p', 'f1p_tf', 'f1p_rfilt')
