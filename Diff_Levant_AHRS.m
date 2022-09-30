clc
clear all 
close all

datos = readmatrix("Datos.txt");

T =0.01;
t = [T:T:length(datos)*T];
f1 = datos(:,2)';
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

tf1 = tf([1 0], conv([0.05 1], [0.05 1]));
f1p_tf = lsim(tf1,f1,t);

tf2 = tf([1], [0.05 1]);
f1p_rfilt = lsim(tf2,f1p_robusto,t);

plot( t, f1p_tf, t, f1p_rfilt, t, f1);
legend('levant', 'levant filtrado', 'real')