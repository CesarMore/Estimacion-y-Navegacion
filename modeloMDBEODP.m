function dxdt = modeloMDBEODP(t,x,tau,d)
%Parámetros
mc=1;
m=0.1;
l=0.5;
g=9.81;

% Definir las matrices f(x),g1(x),g2(x) (VECTORES COLUMNA)

a=(m*l*x(2)*x(2)*cos(x(1))*sin(x(1)));
b= (l*(4/3-(m*cos(x(1))*cos(x(1))/(mc+m))));
fx=[x(2);
   ((g*sin(x(1)))- a)/(l*b)];

g1=[0;
   ((cos(x(1))*sin(x(1)))/(mc+m))/(l*(4/3-((m*cos(x(1))*cos(x(1)))/(mc+m))))];

g2=[0;
    1];

%x1= theta
%x2= thetapunto 
dxdt(1:2,1)=fx+g1*tau+g2*d; % se va a obtener theta y theta punto
                                             %x(1)  y  x(2)

%¿DEBE ENTRAR d AL MODELO?
%d=0.2*sin(x(i+1,1)*x(i+1,2));
% o debe entrar dgorrito =dg




