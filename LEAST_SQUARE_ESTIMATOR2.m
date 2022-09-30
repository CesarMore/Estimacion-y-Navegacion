%Ejemplo 3.5
%--------> CON RUIDO GAUSSIANO CONSTANTE

clc
clear all
close all

T = 1;
t = [0:T:100];               % Intervalo de tiempo
x = [10; 5];                 % valores reales de x
xhat = [8; 7];               % x_gorrito en cero
P = eye(2);                  % Matriz de covarianza P_0
R = 0.01;
%Proponer el ruido Gausiano variante en el tiempo
% a = 1;
% b = 0;
% v = a*randn(1000,1) + b;
%Proponer el ruido Gausiano constante
v = sqrt(R)*randn;

for k=1:length(t)-1
%La matriz Hk es variable    
H = [1 0.99^(k-1)];                %Matriz Hk
% Calcular la medición con ruido
 y = H*x + v;                      %Valores de salida con ruido NO  CTE
% y = H*x + v(k);                  %Valores de salida con ruido CTE
 % 1.- Calcular Kk
Kk = P*(H)'*inv(H*P*(H)'+ R);      %Matriz Kk de ganancias
 % 2.- Calcular xgorrito
xhat = xhat + Kk*(y - H*xhat);     %Valor real estimado
 % 3.- Actualizar la matriz de covarianza
 P = P - Kk*H*P;                   %Matriz de covarianza
 %Separar las variables para poder graficar
xhat_1(k) = xhat(1);
xhat_2(k) = xhat(2);
P_1(k) = P(1,1);
P_2(k) = P(2,2);
end

%Grafica de la covarianza
figure(1)
plot(P_1)
ylim([0 1])
hold on 
plot(P_2)
legend('P(1,1)','P(2,2)') % poner leyendas en la gráfica
ylabel('Varianza')
xlabel('Tiempo [seg]')
%Grafica  de x estimada
figure(2)
plot(xhat_1)
hold on 
plot(xhat_2)
legend('x_1 estimada','x_2 estimada') % poner leyendas en la gráfica
ylabel('Valores estimados de X')
xlabel('Tiempo [seg]')
