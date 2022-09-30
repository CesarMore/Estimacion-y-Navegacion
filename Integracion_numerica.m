clc 
clear all 
close all 

T1 = 0.1;
T2 = 0.05;
T3 = 0.025;
t  = [0:T1:5];
x  = 0.5;
x1 = 0.5;
x2 = 0.5;
%Integración rectangular
for i=1:length(t)-1
    x(i+1) = x(i) + T1*(cos(i)); 
end
%Integración trapecio
for i=1:length(t)-1
    x1(i+1) = x1(i) + T2*(cos(i) + cos(i+1)); 
end
%Integración Runge-Kutta
for i=1:length(t)-1
    x2(i+1) = x2(i) + T3*(cos(i) + cos(i+1) + cos(i+2) +cos(i+3)); 
end
%Gráficas
figure
subplot(2,2,1); hold on; box on;
plot(t,x); title('Integración rectangular', 'FontSize', 12);
xlabel('Tiempo (Segundos)'); ylabel('Área [m^2]');

subplot(2,2,2); hold on; box on;
plot(t,x1); title('Integración trapecio', 'FontSize', 12);
xlabel('Tiempo (Segundos)'); ylabel('Área [m^2]'); 

subplot(2,2,3); hold on; box on;
plot(t,x2); title('Integración Runge-Kutta', 'FontSize', 12);
xlabel('Tiempo (Segundos)'); ylabel('Área [m^2]'); 
