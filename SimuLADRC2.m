clc
clear all
close all

ti = 0;
tf = 60;
T = 0.01;

t = [ti:T:tf];

x = [1 0 0 0 0 0];

U1 = 0;
U2 = 0;
U3 = 0;

z1p = 0;
z2p = 0;
z3p = 0;
z4p = 0;
z5p = 0;
z6p = 0;
z7p = 0;
z8p = 0;
z9p = 0;

z1 = 0;
z2 = 0;
z3 = 0;
z4 = 0;
z5 = 0;
z6 = 0;
z7 = 0;
z8 = 0;
z9 = 0;


wo2=60; %pitch
wo1=50; %roll
wo3=50; %yaw

beta1=3*wo2;       %pitch
beta2=3*wo2^2;
beta3=wo2^3;

beta4=3*wo1;        %roll
beta5=3*wo1^2;
beta6=wo1^3;

beta7=3*wo3;         % yaw
beta8=3*wo3^2;
beta9=wo3^3;

kax=0.008;
kay=0.008;
Jz=0.1; 
kaz=0.009;
l=0.1;
Jx=0.0552;
Jy=0.0552;
ku=0.1188;
ks=0.0036;

wc1=40;
wc2=25;
wc3=15;
kp1=wc1^2;
kd1=2*0.5*wc1;
kp2=wc2^2;
kd2=2*1*wc2;
kp3=wc3^2;
kd3=2*0.7*wc3;
r02p=0;

% r1 = sin(t);
% r2 = cos(t);

A = -4;
w = 0.251;
r1 = A*square(w*t);
r2 = 0;

% b1=0.426;
% b2=0.0327;
% b3=0.426;
%varias las b
b1=1;
b2=1;
b3=0.1;

for i=1:length(t)-1
    [tt,xx]=ode45(@modeloLADRC,[t(i) t(i+1)],x(i,:),[],U1(i),U2(i),U3(i));
    x(i+1,:) = xx(end,:);
    %Observador de estado extendido
    z1p(i+1) = z2(i)-beta1*(z1(i)-x(i+1,1));
    z2p(i+1) = (-kay*z2(i))/Jy+z3(i)+(l*ku*U2(i))/Jy-beta2*(z1(i)-x(i+1,1));
    z3p(i+1) = -beta3*(z1(i)-x(i+1,1));
    z4p(i+1) = z5(i)-beta4*(z4(i)-x(i+1,3));
    z5p(i+1) = (-kax*z5(i))/Jx+z6(i)+(l*ku*U1(i))/Jx-beta5*(z4(i)-x(i+1,3));
    z6p(i+1) = -beta6*(z4(i)-x(i+1,3));
    z7p(i+1) = z8(i)-beta7*(z7(i)-x(i+1,5));
    z8p(i+1) = (-kaz*z8(i))/Jz+z9(i)+(ks*U3(i))/Jz-beta8*(z7(i)-x(i+1,5));
    z9p(i+1) = -beta9*(z7(i)-x(i+1,5));
    %Derivada del observador extendido
    z1(i+1) = T*z1p(i+1)+z1(i);
    z2(i+1) = T*z2p(i+1)+z2(i);
    z3(i+1) = T*z3p(i+1)+z3(i);
    z4(i+1) = T*z4p(i+1)+z4(i);
    z5(i+1) = T*z5p(i+1)+z5(i);
    z6(i+1) = T*z6p(i+1)+z6(i);
    z7(i+1) = T*z7p(i+1)+z7(i);
    z8(i+1) = T*z8p(i+1)+z8(i);
    z9(i+1) = T*z9p(i+1)+z9(i);
    
%     u01(i+1) = kp1*(r1(i)-z4(i))+kd1*(r2(i)-z5(i))+r02p;   %roll
%     u02(i+1) = kp2*(r1(i)-z1(i))+kd2*(r2(i)-z2(i))+r02p;   %pitch
%     u03(i+1) = kp3*(r1(i)-z7(i))+kd1*(r2(i)-z8(i))+r02p;   %yaw

    u01(i+1) = kp1*(r1(i)-z4(i))+kd1*(r2-z5(i))+r02p;   %roll
    u02(i+1) = kp2*(r1(i)-z1(i))+kd2*(r2-z2(i))+r02p;   %pitch
    u03(i+1) = kp3*(r1(i)-z7(i))+kd1*(r2-z8(i))+r02p;   %yaw

    U1(i+1) = (u01(i)-z6(i))/b1;
    U2(i+1) = (u02(i)-z3(i))/b2;
    U3(i+1) = (u03(i)-z9(i))/b3;

end

figure(1), plot(t,x(:,1),t,r1)
hold on
legend('roll','r0')
ylabel('Roll [rad/s]')
xlabel('Tiempo [seg]')
grid
hold off
figure(2), plot(t,x(:,3),t,r1)
hold on
legend('pitch','r0')
ylabel('Pitch [rad/s]')
xlabel('Tiempo [seg]')
grid
hold off
figure(3), plot(t,x(:,5),t,r1)
hold on
legend('yaw','r0')
ylabel('Yaw [rad/s]')
xlabel('Tiempo [seg]')
grid
hold off
