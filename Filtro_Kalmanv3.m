% Filtro Kalman: Ejemplo 5.1
clc
clear all
close all

T = 1;
t =[0:T:60];                
x = [35; 7; 5];
xhat = [40; 2; 2.5];
P = eye(3); 
R = 30^2;
a = 0.05;
b = 0;
%v = a*randn(length(t),1) + b; 

for k=1:length(t)-1
    H = [1  0  0];   
    %y = H*x + v(k);                  
    y = H*x;
    y = awgn(y,0.5,'measured');
    F=[1 T (T^2)/2; 0 1 T; 0 0 1]; 
    xhat=F*xhat; 
    xg_1(k)=xhat(1);
    error_estimacion(k)= x(1)-xhat(1);
    P=F*P*(F)';
    PP(k)=P(1,1);
    Kk = P*(H)'*inv(H*P*(H)'+ R);    
    xhat = xhat + Kk*(y - H*xhat);    
    xhat_1(k) = xhat(1); 
    error_medicion(k)= x(1)-xhat(1);
    P = P - Kk*H*P;                  
    PP(k)=P(1,1);

end

figure(1)
plot(PP)
ylabel('Varianza')
xlabel('Tiempo [seg]')
figure(2)
plot(xhat_1)
hold on 
plot(xg_1)
legend('x a posteriori','x a priori') 
ylabel('Valores estimados de X')
xlabel('Tiempo [seg]')
figure(3)
plot(error_medicion)
hold on 
plot(error_estimacion)
legend('Error de medición','Error de estimacion') 
ylabel('Errores')
xlabel('Tiempo [seg]')
figure(5)
plot(xhat_1)
legend('x a posteriori') 
ylabel('Valores estimados de X')
xlabel('Tiempo [seg]')