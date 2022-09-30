
% CONTROL POR MODOS DESLIZANTES BASADO EN UN OBSERVADOR DE PERTURBACIONES
clc
clear 
close all
%Declarar  las variables y sus valores
T=0.01; %Intervalo de muestreo
t=[0:T:2] % Tiempo inicial y final  con 0.1 entre tiempos
x=[2 1]; %valores iniciales X0

%Parámetros
mc=1;    % Masa del péndulo
m=0.1;   % Mada del vehiculo
l=0.5;   % Longitud
g=9.81;  % Gravedad

%---------------------------------------------------------------
tau=0;      % valor inicial de tau (u) 
dg= -1;     % Vector inicial de la perturbación estimada
t1d=0;      % Valores deseados de theta 
t2d=0;      % Valores deseados de theta punto
xdp=[0;
     0];    % Vector de theta deseada punto
K= 1 ;      % En el artículo no viene que valor se uso, se va a proponer
d=0;        % Valor de perturbacion real inicial

for i=1:length(t)-1
[tt,xx]= ode45(@modeloMDBEODP,[t(i) t(i+1)],x(i,:),[],tau(i),d(i));
%                           [x:y:z]      [a,b],[],[tau]
% ODE45 = integrar la ecuacion diferencial
%[tt,xx]  tt= vector de tiempo xx= vector o matriz que devuelve el estado
%[x:y:z] x= inicio del tiempo y=intervalo de tiempo z= tiempo final
%,[a,b] = condiciones iniciales de x1=1 y x2=0.1
x(i+1,:)= xx(end,:); % Se va a guardar solamente el último renglon en el tiempo de muestreo T=0.1
%x(i+1,:)=[x] - El vector esta formado por las componentes

%---------------------------------------------------------------
%----------------- Superficie de deslizamiento -----------------

e1=x(i+1,1)-t1d;  % error 1
e2=x(i+1,2)-t2d;  % error 2
C=[1 ,1];
S= C*[e1;
      e2];

% ---------------------- PERTURBACION ---------------------------
d(i+1)=0.2*sin(x(i+1,1)*x(i+1,2)); % Perturbación REAL

%------------------------- OBSERVADOR ---------------------------
lx=[0 1];
px=x(i+1,2);
Bx=1;

a=(m*l*x(i+1,2)*x(i+1,2)*cos(x(i+1,1))*sin(x(i+1,1)));
b= (l*(4/3-(m*cos(x(i+1,1))*cos(x(i+1,1))/(mc+m))));
fx=[x(i+1,2);
   ((g*sin(x(i+1,1)))- a)/(l*b)];

g1=[0;
   ((cos(x(i+1,1))*sin(x(i+1,1)))/(mc+m))/(l*(4/3-((m*cos(x(i+1,1))*cos(x(i+1,1)))/(mc+m))))];

g2=[0;
    1];

z=dg-px;                            %ecuación (12)
zp=-lx*fx-lx*g1-lx*g2*z-lx*g2*px;   %ecuación (11)

%¿PARA OBTENER dgorrito (dg) SE VA A INTEGRAR ZP?

%----------------------------------------------------------------
%Ley de control  tau (u) 
j=inv(C*g1);
tau(i+1)=-j*(C*fx-C*xdp+C*g2*dg+K*sign(S)+ (abs(S*C*g2)^2)/Bx);

%----------------------------------------------------------------
% Calcular la salida (y)
%y=x1=theta
y(i+1)=x(i+1,1);  
end


%Grafica de los estados
figure(1)
plot(t,x(:,1)) % generar la gráfica
hold on
plot(t,x(:,2)) % generar la gráfica
hold on
legend('x1[Grados]','x2[Grados/s]') % poner leyendas en la gráfica
ylabel('Estados')
xlabel('Tiempo [seg]')
grid
hold off

%Grafica de la perturbación real (d)
% d=0.2*sin(x(1)*x(2))
figure(2)
plot(t,d) % generar la gráfica
legend('d') % poner leyendas en la gráfica
ylabel('Perturbación real (d)')
xlabel('Tiempo [seg]')
grid

%Grafica de la ley de contro de tau (u)
figure(3)
plot(t,tau) % generar la gráfica
legend('u [N]') % poner leyendas en la gráfica
ylabel('Entrada de control (u)')
xlabel('Tiempo [seg]')
grid




