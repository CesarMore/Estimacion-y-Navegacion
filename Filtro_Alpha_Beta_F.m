clc
clear all
close all

T = 0.05;              %Intervalo de muestreo
t = [0:T:119.85];    %Vector de tiempo
tam_t=length(t)

%------------------- MEDICIONES CON RUIDO -----------------
datos1 = csvread('datos4.csv');
%PONER PERIODO DE MUESTREO
y_sf = datos1(:,1);          %Posición con ruido
xref = datos1(:,2);       %Posición real
%hist(y)
media_y= mean(y_sf)
varianza_y=var(y_sf)
desviacion_y=std(y_sf)
%--------------- Filtrar la señal de salida y --------------
filter1.gain = 0.0011381669272451226;
filter1.numerator 	= [1, 2, 1] * filter1.gain;
filter1.denominator = [1, -1.902329, 0.906882];

y=filter(filter1.numerator,filter1.denominator,y_sf);

%-------------------Ruido de medición ----------------------
VAR = 5;         %Varianza para el ruido
MED =0;             %Media del ruido
v=sqrt(3*(VAR))*((2*rand(length(t),1)-1)); %Generar el ruido
%hist(v)
media= mean(v)
varianza=var(v)
desviacion=std(v) %sqrt(varianza)
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
x = [-2 0.1 ];        %Valores iniciales de x
H = [1 0];            %Matriz H

R=VAR;

xmas = x;
% p_mas = 1*eye(2);

%Calcular todos los datos del sistema con ruido
for i=2:length (t)
    x(i,:)=F*x(i-1,:)'+w(:,i-1);
    %y no se debe calcular ya se conoce
end

%Calcular K
    lambda=(VAR1*T^2)/VAR;
    
    k1=-1/8*(lambda^2+8*lambda-(lambda+4)*sqrt(lambda^2+8*lambda));
    k2=(1/(4*T))*(lambda^2+4*lambda-lambda*sqrt(lambda^2+8*lambda));
    K=[k1;k2];
    %------------------------------
    p11_min=(k1*VAR1)/(1-k1);
    p12_min=(k2*VAR1)/(1-k1);
    p22_min=((k1/T)+(k2/2))*p12_min;
    
    p_min=[ p11_min  p12_min;
            p12_min  p22_min];
   %--------------------------------   
    p11_mas=k1*VAR;
    p12_mas=k2*VAR;
    p22_mas=((k1/T)-(k2/2))*p12_min;
    
    p_mas=[p11_mas  p12_mas;
           p12_mas  p22_mas];
    

for i=2:length(t)
    
    %Prediccion 
    xmenos(i,:)=F*xmas(i-1,:)';
   
    %Actualizacion
    
    xmas(i,:)=xmenos(i,:)+K'*(y(i)-H*xmenos(i,:)');
end

%Meter filtro para la posición con alpha-beta
filter1.gain = 0.04187682823477176;
filter1.numerator 	= [1, 3, 3, 1] * filter1.gain;
filter1.denominator = [1, -1.280319, 0.775103, -0.159769];
xf_alpha=filter(filter1.numerator,filter1.denominator,xmas(:,1));

%Graficar el ruido 
% figure(4) 
% plot(t,w2)
% title ('Ruido uniforme')
% %gráfica de la posición
% figure(1)
% plot(t,y,t,xref)
% legend('Posición(ruido)','Posición(real)')

%gráfica de la posición
figure(1)
plot(t,y,t,xmas(:,1),t,xref)
legend('Posición(ruido)','Posición(estimada)','Posición(referencia)')

figure(9)
plot(t,y,t,y_sf,t,xref)
legend('y=Posición(ruido) filtro','y=Posición(ruido)sin filtro','x_{ref}')

% figure(8)
% plot(t,y,t,xmas(:,1),t,xref,t,xf_alpha(:,1),t,y_sf)
% legend('y=Posición(ruido)filtro','xg=Posición(estimada)','x_r=Posición(referencia)','xgf=Posición(estimada)filtro','y=Posición(ruido)sin filtro')

% %Gráfica de la velocidad
% figure(2)
% subplot(211)
% plot(t,xmas(:,2),t,x(:,2))

% figure(3)
% plot(t,x(:,1)-y',t,x(:,1)-xmas(:,1))
% legend('E_{medicion}','E_{estimado}')
% plot(t,x(:,1)-y',t,x(:,1)-xmas(:,1))

% figure(4)
% hold on
% for i=1:4
%     plot([i i+1],[pplus(1,1,i) pminus(1,1,i+1)],'r','LineWidth',2)
%     plot([i+1 i+1],[pminus(1,1,i+1) pplus(1,1,i+1)],'r','LineWidth',2)
% end