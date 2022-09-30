%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Control de altidud de un cuadrirotor...
% basado en el metodo LADRC 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all 
close all
ti = 0;
tf = 60;
T = 0.01;
t = [ti:T:tf];
%---------Parametros
b1=0.426;
b2=0.0327;
b3=0.426;
wc1=10;
wc2=30;
wc3=15;
kp1=wc1^2;
kd1=2*0.5*wc1;
kp2=wc2^2;
kd2=2*0.5*wc2;
kp3=wc3^2;
kd3=2*0.7*wc3;
r02p=0;
%
wo2=6; %pitch
wo1=5; %roll
wo3=5; %yaw
%Valores de beta para pitch
beta1=3*wo2;
beta2=3*wo2^2;
beta3=wo2^3;
%Valores de beta para roll
beta4=3*wo1;
beta5=3*wo1^2;
beta6=wo1^3;
%Valores de beta para yaw
beta7=3*wo3;
beta8=3*wo3^2;
beta9=wo3^3;
%
l=0.1;
kax=0.008;
kay=0.008;
Jz=0.1; 
kaz=0.009;
Jx=0.0552;
Jy=0.0552;
ku=0.1188;
ks=0.0036;
%%%%%%%%%%%%%%%%%%%%%%%%
x = [0 0 0 0 0 0 0 0 0];
z1p = 0;
z2p = 0;
z3p = 0;
z4p = 0;
z5p = 0;
z6p = 0;
z7p = 0;
z8p = 0;
z9p = 0;
U1 = 0;
U2 = 0;
U3 = 0;
uo1 = 0;
uo2 = 0;
uo3 = 0;
d1 = 5;
d2 = 5;
d3 = 5;
A = -4;
w = 0.251;
r0 = A*square(w*t);
e1 = 0;
e2 = 0;
z1 = 0;
z2 = 1;
z3 = 0;
z4 = 0;
z5 = 0;
z6 = 0;
z7 = 0;
z8 = 0;
z9 = 0;

for i=1:length(t)-1
    [tt,xx,zz] = ode45(@modeloLADRC,[t(i) t(i+1)],x(i,:),[],U1(i),U2(i),U3(i));
    x(i+1,:) = xx(end,:);
    %z(i+1,:) = zz(end,:);

    %%Observador de estado extendido
    z1p(i+1) = z2(i)-beta1*(z1(i)-x(i+1,1));
    z2p(i+1) = -(kay/Jy)*z2(i)+z3(i)+l*(ku/Jy)*U2(i)-beta2*(z1(i)-x(i+1,1));
    z3p(i+1) = -beta3*(z1(i)-x(i+1,1));
    z4p(i+1) = z5(i)-beta4*(z4(i)-x(i+1,4));
    z5p(i+1) = -(kax/Jx)*z5(i)+z6(i)+l*(ku/Jx)*U1(i)-beta5*(z4(i)-x(i+1,4));
    z6p(i+1) = -beta6*(z4(i)-x(i+1,4));
    z7p(i+1) = z8(i)-beta7*(z7(i)-x(i+1,7));
    z8p(i+1) = -(kaz/Jz)*z8(i)+z9(i)+(ks/Jz)*U3(i)-beta8*(z7(i)-x(i+1,7));
    z9p(i+1) = -beta9*(z7(i)-x(i+1,7));
    %Integrar zp
    z1(i+1) = T*(z1p(i+1))+z1(i);
    z2(i+1) = T*(z2p(i+1))+z2(i);
    z3(i+1) = T*(z3p(i+1))+z3(i);
    z4(i+1) = T*(z4p(i+1))+z4(i);
    z5(i+1) = T*(z5p(i+1))+z5(i);
    z6(i+1) = T*(z6p(i+1))+z6(i);
    z7(i+1) = T*(z7p(i+1))+z7(i);
    z8(i+1) = T*(z8p(i+1))+z8(i);
    z9(i+1) = T*(z9p(i+1))+z9(i);
    %Referencia 
    %%Control PD uoi
    r0p = 0;
    % Para pitch
    uo2(i+1)=kp2*(r0(i)-z1(i+1)) + kd2*(r0p-z2(i+1)) + r02p;
    % Para roll
    uo1(i+1)=kp1*(r0(i)-z4(i+1)) + kd1*(r0p-z5(i+1)) + r02p;
    % Para yaw
    uo3(i+1)=kp3*(r0(i)-z7(i+1)) + kd3*(r0p-z8(i+1)) + r02p;
    %%Control Ui
    % Para roll
    U1(i+1)=(uo1(i)-z6(i))/b1;
    %Para pitch
    U2(i+1)=(uo2(i)-z3(i))/b2;
    %para yaw
    U3(i+1)=(uo3(i)-z9(i))/b3;
    %%%Error de seguimiento 
    %e1(i+1) = r0(i) - z1(i);
    %e2(i+1) = r0p(i) - z2(i);
    %%Salida
    y(i+1) = x(i+1,1); 
end
%Grafica de los estados
figure(1)
plot(t,x(:,1),'LineWidth',1) 
hold on
plot(t,z1,'LineWidth',1)
hold on
legend('x1[Grados]','z1[Grados/s]')
ylabel('Estados')
xlabel('Tiempo [seg]')
grid
hold off
figure(2)
plot(t,x(:,1),'LineWidth',1) 
hold on
plot(t,r0,'LineWidth',1)
hold on
legend('x1(roll)','ro(referencia)')
ylabel('Estados')
xlabel('Tiempo [seg]')
grid
hold off
figure(3)
plot(t,uo1)
hold on
legend('uo1')
% ylabel('Estados')
% xlabel('Tiempo [seg]')
grid
hold off