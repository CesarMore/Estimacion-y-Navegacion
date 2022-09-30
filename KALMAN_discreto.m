clc
clear all
close all

T = 5;              %Intervalo de muestreo
t = [0:0.1:150];    %Vector de tiempo

%----------------------- Ruido de medición -------------------------
sigma_r=1000^2;      %Varianza para el ruido v1
sigma_t= 0.017^2;    %Varianza para el ruido v2
MED =0;              %Media del ruido
v1=sqrt(0.25*(sigma_r))*((2*randn(length(t),1)-100));
v2=sqrt(0.5*(sigma_t))*((2*randn(length(t),1)-1));
media_v=mean(v2)
varianza_v=var(v2)
% figure(20)
% hist(v)
v=[v1;
   v2]; 
%-------------------- Ruido del sistema ------------------------------
sigma_1=(103/3)^2;    %Varianza para el ruido w1
sigma_2=0.000000013;  %Varianza para el ruido w2  
MED =0;               %Media del ruido
w1=sqrt(0.25*(sigma_1))*((2*randn(length(t),1)-0.009999)); %Generar el ruido w1
w2=sqrt(0.25*(sigma_2))*((2*randn(length(t),1)-1)); %Generar el ruido w2
%hist(v)
media_w= mean(w2)
varianza_w=var(w2)
%--------------------------- Sistema --------------------------------

x = [0 0 0 0 0 0];                %Valores iniciales de x

H = [1 0 0 0 0 0;
     0 0 0 1 0 0];                %Matriz H

                           
 %y = H*x';                        %Valores medidos con ruido

rho=0.5;

F=[1 T  0  0 0  0;
   0 1  1  0 0  0;
   0 0 rho 0 0  0;
   0 0  0  1 T  0;
   0 0  0  0 1  1;
   0 0  0  0 0 rho];  


%----------------- Parámetros e inicializacion -----------------------

Pmas=[sigma_r          (sigma_r)/T              0        0                          0                   0;
     (sigma_r)/T  ((2*sigma_r)/(T^2))+sigma_1   0        0                          0                   0;
     0                         0                  sigma_1    0                          0                   0;
     0                         0                      0     sigma_t               (sigma_t)/T             0; 
     0                         0                      0    (sigma_t)/T   ((2*sigma_t)/(T^2))+sigma_2    0;
     0                         0                      0        0                          0               sigma_2];


Q=[ 0 0    0       0 0     0;
    0 0    0       0 0     0;
    0 0 sigma_1  0 0     0;
    0 0    0       0 0     0;
    0 0    0       0 0     0;
    0 0    0       0 0  sigma_2];

R=[sigma_r   0;
     0      sigma_t];

xmas = x;

w=[0;
   0;
   w1;
   0;
   0;
   w2];
%Calcular todos los datos del sistema  
for i=2:length (t)
    x(i,:) = F*x(i-1,:)'+w(i); % (real)
    y(i,:) = H*x(i,:)'+v(i); % (Ruido)
end

%Calcular todos los datos (estimados)
for i=2:length(t)
    
    %Prediccion
    
    xmenos(i,:)=F*xmas(i-1,:)';
    Pmenos = F*Pmas*F'+Q;
    pminus(:,:,i)=Pmenos;
    
    %Actualizacion
    
    K = Pmenos*H'*inv(H*Pmenos*H'+R);
    Pmas = Pmenos - K*H*Pmenos;
    Pplus(:,:,i)=Pmas;
    xmas(i,:) = xmenos(i,:)'+K*(y(i)-H*xmenos(i,:)');
end

%Grafica de la salida (posición) con ruido y la estimada
figure(1)
plot(t,y(:,1),t,xmas(:,1),t,x(:,1))
ylabel(' Range of the vehicle  ')
xlabel('Tiempo[segundos]')
legend('Salida y','r_{estimado}','r_{real}')

%Grafica de la velocidad  con ruido y la estimada
figure(2)
plot(t,xmas(:,2),t,x(:,2))
ylabel(' Range rate of the vehicle ')
xlabel('Tiempo[segundos]')
legend('rp_{estimado}','rp_{real}')

%Grafica del error de posición
figure(3)
e_m = x(:,1)-y(:,1);
e_e = x(:,1)-xmas(:,1);
plot(t,e_m,t,e_e)
legend('E_{medicion}','E_{estimado}')

%Grafica de la varianza
figure(4)
hold on
for i=1:30
    plot([i i+1],[Pplus(1,1,i) pminus(1,1,i+1)],'r','LineWidth',2)
    plot([i+1 i+1],[pminus(1,1,i+1) Pplus(1,1,i+1)],'r','LineWidth',2)
end