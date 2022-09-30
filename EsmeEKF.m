% ------------------- FILTRO DE KALMAN CONTINUO ----------------------
% Ejemplo 13.1  pag.410
clc
clear all
close all

T = 0.01;            %Intervalo de muestreo
t = [0:T:10];    %Vector de tiempo
 %------------------ Parámetros ----------------------
 Resis=2;       % Se cambio R por Resis donde R= resistencia
 Induc=0.03;       % Se cambio L por Induc donde L= inductancia
 lambda=0.1;
 F=0.0015;
 J=0.002;
 %----------------------- Ruido en la salida ----------------------
 media1=0;
 desviacion1=0.1;
 VAR=0.1^2;
 v1=sqrt(0.28*(VAR))*((2*randn(length(t),1)-1));
 v2=sqrt(0.28*(VAR))*((2*randn(length(t),1)-1));
 media_v1=mean(v1);
 varianza_v1=var(v1);
 media_v2=mean(v2);
 varianza_v2=var(v2);
 % hist(v)
 %-------------------- Ruidos del sistema --------------------------
 % Ruido de la entrada de control 
 media2=0;
 desviacion2=0.01;
 VAR1=0.01^2;
 q1=sqrt(0.28*(VAR1))*((2*randn(length(t),1)-0.2)); %Generar el ruido
 q2=sqrt(0.28*(VAR1))*((2*randn(length(t),1)-0.2)); %Generar el ruido
 media_q1= mean(q1);
 varianza_q1=var(q1);
 media_q2= mean(q2);
 varianza_q2=var(q2);
 % Ruido q3
 media3=0;
 desviacion3=0.5;
 VAR2=0.5^2;
 q3=sqrt(0.28*(VAR2))*((2*randn(length(t),1)-0.2)); %Generar el ruido
 media_q3= mean(q3);
 varianza_q3=var(q3);
%------------------------Señales con ruido ---------------------
ua=sin(2*pi*t);
ub=cos(2*pi*t);
%---------------------- Matrices de covarianza -----------------
%La varianza = (desviación estandar)^2
% Q=[0.01^2  0    0   0;
%     0   0.01^2  0   0;
%     0    0   0.5^2  0;
%     0    0    0     0];
  Q=[0.01^2  0    0 ;
     0   0.01^2  0  ;
     0    0   0.5^2 ];

R= [0.1^2  0;
     0  0.1^2];

%x=[ia ib omega theta]'   Vector de estados
%--------------------------------------------------------------------------
%Valores iniciales para la integración
x  =[0 0 0 0];   %Valores iniciales de x (estados)
x_e=[0;0;0;0];
ia_n=0;
ib_n=0;
omega_n=0;
theta_n=0;
xp_e=[0;0;0;0];
xp=[0;0;0;0];
f=[0;0;0;0];     %Valores iniciales de f
y=[0;0];         %Valores iniciales de y
y0=[0;0];        %Valores iniciales de y0

%-------------------------------------
% L=[1/Induc   0     0   0;
%      0    1/Induc  0   0;
%      0       0    1/J  0;
%      0       0     0   0];
 L=[1/Induc   0     0  ;
     0    1/Induc  0   ;
     0       0    1/J  ;
     0       0      0];
 
C=[1 0 0 0;
   0 1 0 0];

M=[1 0;
   0 1];

P=eye(4);        %Valores iniciales de la Matriz P


Qtilde=L*Q*L';
Rtilde=M*R*M';

v=[v1;
   v2];

%INTEGRAR POR ODE45 PARA OBTENER x_real (Modelo original)
for i=1:length(t)-1            
    [tt2, xx2] = ode45(@modelo_x,[t(i) t(i+1)],x(i,:),[],ua(i),ub(i),q1(i),q2(i),q3(i));
    x(i+1,:) = xx2(end,:);
end


% for i=i:length(t)-1
% %Hacer la integración rectangular para el sistema original
%     x(:,i+1)= x(:,i)+T*xp(:,i);
%     ia(i+1)=x(1);
%     ib(i+1)=x(2);
%     omega(i+1)=x(3);
%     theta(i+1)=x(4);
%     
%     y1(i+1)=ia(i)+v1(i);
%     y2(i+1)=ib(i)+v2(i);
%     
%     y=[y1;
%        y2];
%  % ------------------- Sistema original ---------------------------
%  % f(x,u,w,t)
%  iap(i+1)=-(Resis/Induc)*ia(i) +(omega(i)*lambda/Induc)*sin(theta(i))+((ua(i)+q1(i))/Induc);
%  ibp(i+1)=-(Resis/Induc)*ib(i) -(omega(i)*lambda/Induc)*cos(theta(i))+((ub(i)+q2(i))/Induc);
%  omegap(i+1)=-(3*lambda/2*J)*ia(i)*sin(theta(i)) +(3*lambda/2*J)*ib(i)*cos(theta(i))-(F*omega(i)/J)+q3(i);
%  thetap(i+1)=omega(i);
%  
%  xp=[iap;
%      ibp;
%      omegap;
%      thetap]; 
% end

%--------------------------------------------------------------------------
%Calcular la estimación por filtro de Kalman extendido
for i=1:length(t)-1
    
    y(:,i+1)=[x(i,1)+v1(i);
              x(i,2)+v2(i)];
    %Calculo de la matriz A
    
    A=[-Resis/Induc                                   0                    lambda*sin(theta_n(i))/Induc             omega_n(i)*lambda*cos(theta_n(i))/Induc;
          0                                      -Resis/Induc             -lambda*cos(theta_n(i))/Induc             omega_n(i)*lambda*sin(theta_n(i))/Induc;
    -3*lambda*sin(theta_n(i))/(2/J)   3*lambda*cos(theta_n(i))/(2/J)           -F/J                     -3*lambda*(ia_n(i)*cos(theta_n(i))+ib_n(i)*sin(theta_n(i)))/(2/J);
          0                                        0                             1                                                   0] ;



Pp= A*P+P*A'+Qtilde-P*C'*inv(Rtilde)*C*P;

P = P + T*(Pp);

K=P*C'*inv(Rtilde);
xp_e=f(:,i)+K*(y(:,i)-y0(:,i));


%funcion = f + K*(y(:,i)-h);

%Hacer la integración rectangular para Pp
    
%Hacer la integración rectangular para xp_e
    x_e(:,i+1)= x_e(:,i)+T*xp_e;
    ia_n(i+1)=x_e(1,i+1);
    ib_n(i+1)=x_e(2,i+1);
    omega_n(i+1)=x_e(3,i+1);
    theta_n(i+1)=x_e(4,i+1);
    %Calculo de la salida nominal    
    y0(:,i+1)=[x_e(1,i+1);
        x_e(2,i+1)];
 % ------------------- Sistema nominal ---------------------------
 % f(xe,u,w0,t)
 iap_n(i+1)=-(Resis/Induc)*ia_n(i) +((omega_n(i)*lambda)/Induc)*sin(theta_n(i))+(ua(i)/Induc);
 ibp_n(i+1)=-(Resis/Induc)*ib_n(i) -((omega_n(i)*lambda)/Induc)*cos(theta_n(i))+(ub(i)/Induc);
 omegap_n(i+1)=-(3*lambda/(2*J))*ia_n(i)*sin(theta_n(i)) +(3*lambda/(2*J))*ib_n(i)*cos(theta_n(i))-(F*omega_n(i)/J);
 thetap_n(i+1)=omega_n(i);
 
 f(:,i+1)=[iap_n(i+1);
    ibp_n(i+1);
    omegap_n(i+1);
    thetap_n(i+1)];

end

%Grafica de la corriente ia real y estimada
figure(1)
plot(t,y(1,:),t,x_e(1,:),t,x(:,1))
ylabel(' Corriente ')
xlabel('Tiempo[segundos]')
legend('Salida y_1','i_a(Estimada)','i_a(Real)')
%Grafica de la corriente ib real y estimada
figure(2)
plot(t,y(2,:),t,x_e(2,:),t,x(:,2))
ylabel(' Corriente ')
xlabel('Tiempo[segundos]')
legend('Salida y_2','i_b(Estimada)','i_b(Real)')

% %Grafica de la velocidad  con ruido y la estimada
% figure(2)
% plot(t,xmas(:,2),t,x(:,2))
% ylabel(' Range rate of the vehicle ')
% xlabel('Tiempo[segundos]')
% legend('rp_{estimado}','rp_{real}')

% %Grafica del error de posición
% figure(3)
% e_m=x(:,1)-y(:,1)
% e_e=x(:,1)-xmas(:,1)
% plot(t,e_m,t,e_e)
% legend('E_{medicion}','E_{estimado}')




% ------------------- Sistema original ---------------------------
function dxdt = modelo_x(t,x,ua,ub,q1,q2,q3)
%------------------ Parámetros ----------------------
 Resis=2;       % Se cambio R por Resis donde R= resistencia
 Induc=3;       % Se cambio L por Induc donde L= inductancia
 lambda=0.1;
 F=0.0015;
 J=0.002;

 % f(x,u,w,t)
 dxdt(1,1)=-(Resis/Induc)*x(1) +(x(3)*lambda/Induc)*sin(x(4))+((ua+q1)/Induc);
 dxdt(2,1)=-(Resis/Induc)*x(2) -(x(3)*lambda/Induc)*cos(x(4))+((ub+q2)/Induc);
 dxdt(3,1)=-(3*lambda/2*J)*x(1)*sin(x(4)) +(3*lambda/2*J)*x(2)*cos(x(4))-(F*x(3)/J)+q3;
 dxdt(4,1)=x(3);

end