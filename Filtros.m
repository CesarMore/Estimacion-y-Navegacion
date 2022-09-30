clc
clear all
close all

datos1 = csvread('datos4.csv');

T = 0.05;
tau = 1;
tauL = 5;
x = datos1(:,1);
xref = datos1(:,2);
t = [0:0.05:119.85];

alpha = .005;  %k2
lambda = 1.9;  %k1
xf(1) = x(1);
x1 = 0;
x2 = 0;
for i=1:length(t)-1
    xf(i+1) = xf(i) + (T/tau)*(x(i) - xf(i));
    %Derivada de Levant
    u = (-lambda*abs(x1(i) - x(i))^(1/2)*sign(x1(i) - x(i)) + x2(i));
    x1(i+1) = x1(i) + T*u;
    x2(i+1) = x2(i) + T*(-alpha*sign(x1(i) - x(i)));
    x_levant(i+1) = u;
   
end

x_levantf = x_levant;
for i=1:length(t)-1
    x_levantf(i+1) = x_levantf(i) + (T/tau)*(x_levant(i) - x_levantf(i));
end


figure(1)
% plot(t,x,t,xf,t,x1_levant) legend('x','xfil','x1 levant')
plot(t,x,t,xf,'LineWidth',1.7)
legend('x','xfil','x1 levant')
set(gca,'FontSize',13)
xlabel('X-axis [m]')
ylabel('Y-axis [m]')


%Filtro alpha-beta
%-------------------Ruido de medición ----------------------
VAR = 5;         %Varianza para el ruido
MED =0;             %Media del ruido
v=sqrt(3*(VAR))*((2*rand(length(t),1)-1)); %Generar el ruido
%hist(v)
media= mean(v);
varianza=var(v);
desviacion=std(v); %sqrt(varianza)
%-------------------Ruido del sistema ----------------------
VAR1 = 0.5;          %Varianza para el ruido del sistema
MED2 =0;             %Media del ruido del sistema
Q=[1/4*(T^4) 1/2*(T^3);
   1/2*(T^3)       T^2]*VAR1^2;
varianza_Q=var(Q);

 w1=sqrt(3*varianza_Q(1,1))*((2*rand(length(t),1)-1)); %Generar el ruido
 w2=sqrt(3*varianza_Q(1,2))*((2*rand(length(t),1)-1)); %Generar el ruido

w=[w1.'; w2.'];

F = [1 T ; 0 1 ];     %Matriz F
x = [-2 0.04];        %Valores iniciales de x
H = [1 1];            %Matriz H
y = H*x';
R=VAR;

xmas = x;

%Calcular K
lambda=(VAR1*T^2)/VAR;
k1=-1/8*(lambda^2+8*lambda-(lambda+4)*sqrt(lambda^2+8*lambda));
k2=(1/(4*T))*(lambda^2+4*lambda-lambda*sqrt(lambda^2+8*lambda));
K=[k1;k2];
p11_min=(k1*VAR1)/(1-k1);
p12_min=(k2*VAR1)/(1-k1);
p22_min=((k1/T)+(k2/2))*p12_min;
p_min=[ p11_min  p12_min; p12_min  p22_min];
%--------------------------------   
p11_mas=k1*VAR;
p12_mas=k2*VAR;
p22_mas=((k1/T)-(k2/2))*p12_min;
p_mas=[p11_mas  p12_mas; p12_mas  p22_mas];

for i=2:length(t)
    %Prediccion 
    xmenos(i,:) = F*xmas(i-1,:)';
    %Actualizacion
    xmas(i,:) = xmenos(i,:) + K'*(y - H*xmenos(i,:)');
end

x_ABp = xmas(2);
for i=1:length(t)-1
    x_ABp(i+1) = (xmas(i+1)-xmas(i))/T;
end

% figure(1)
% % plot(t,x,t,xf,t,x1_levant) legend('x','xfil','x1 levant')
% plot(t,x,t,xf,'LineWidth',1.7)
% legend('x','xfil','x1 levant')
% set(gca,'FontSize',13)
% xlabel('X-axis [m]')
% ylabel('Y-axis [m]')

figure(2)
plot(t,xf,t,xref)
title('X filtrada')

xp = 0;
for i=1:length(t)-1
    xp(i+1) = (xf(i+1)-xf(i))/T;
end
xp2 = diff(x)/T;

figure(3)
plot(t,xp)
legend('xp_E')

xpf(1) = xp(1);
for i=1:length(t)-1
    xpf(i+1) = xpf(i) + (T/tau)*(xp(i) - xpf(i));    
end
xrefp = diff(xref)/T;

figure(4)
plot(t,xpf,t(1:end-1),xrefp,t,x_levantf)
legend('xpf Euler','xp Ref')
title('X punto Filtrada')

figure(5)
plot(t,xp,t,x_levantf,t, x_ABp)
legend('xp, x levantf','x_ABF')
