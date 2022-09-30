
function dxdt = modeloQuad2(t,x,tau1,tau2,tau3,OMEGA)
%function dxdt = modeloQuad2(t,x,tau1,tau2,tau3)

l=0.1969;
Kax=0.008;
Kay=0.008;
Jz=0.1104; 
Kaz=0.009;
Jx=0.0552;
Jy=0.0552;
ku=0.1188; %ku=kt
ks=0.0036;
Jm=0.1;
% JM=Jr= Inercia del motor
% omega= velocidad angular
%omega=omega1+omega2+omega3+omega4;
% p=phip
% q=thetap
% r=psip

%MODELO DINÁMICO

%Sean los cambios de variable
% x1=theta  ---> x1p=x2
% x2=thetap  ---> x2p=theta2p
% x3=phi     ---> x3p=x4
% x4=phip   ----> x4p=phi2p
% x5=psi     ---> x5p=x6
% x6=psip   ---->x6p=psi2p

% Modelo para pitch
%theta2p=(Jz-Jx)/Jy*phip*psip+Jm*omega/Jy+tau2/Jy-Kay*thetap;
dxdt(1,1) = x(2);
dxdt(2,1) =(1/Jy)*((Jz-Jx)*x(4)*x(6)+Jm*OMEGA+tau2-Kay*x(2));
% Modelo para roll
% phi2p=(Jy-Jz)/Jx*thetap*psip+Jm*omega/Jx+tau1/Jx-Kax*phip;
dxdt(3,1) = x(4);
dxdt(4,1) =(1/Jx)*((Jy-Jz)*x(2)*x(6)+Jm*OMEGA+tau1-Kax*x(4));
% Modelo para yaw
% yaw2p=(Jx-Jy)/Jz*phip*thetap+tau3/Jz-Kaz*psip;
dxdt(5,1) = x(6);
dxdt(6,1) = (1/Jz)*((Jx-Jy)*x(4)*x(2)+tau3-Kaz*x(6));




