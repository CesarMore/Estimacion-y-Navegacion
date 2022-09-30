clc
clear all
close all

global K C B A u zeta omega k1 k2


T = 0.01;
t = [0:T:10];
%t1 =[0:100];

%------------------------
VAR = 4.47;
VAR1 = 0.00001;
w=sqrt(0.25*(VAR))*((2*randn(length(t),1)-1));
v=sqrt(0.25*(VAR1))*((2*randn(length(t),1)-1));
% w=0*w;
% v=0*v;
%------------------------

omega = 5;
zeta = 0.2;
u = 0;


x = [0 0];
xhat = [0 0];
P = [2 0;0 2];

A = [0 1;-omega^2 -2*zeta*omega];
B = [0;0];
C = [1 0];
Q = var(w);
R = var(v);
Qc = Q/T;
Rc = R/T;

%---------P

P=[2 0 0 2];

[tt, PP] = ode45(@modeloPP,[0 T],P,[]);
z(:) = PP(end,:);

P=[z(1) z(2);
   z(3) z(4)]

figure(3)
plot(tt,PP(:,1),tt,PP(:,2),tt,PP(:,4))
%-----------------------------

%-------- K
K = P*transpose(C)*inv(Rc);
k1=K(1);
k2=K(2);
%-------------


for i=1:length(t)-1
    [tt,xx] = ode45(@modeloKalmanContinuo,[t(i) t(i+1)],x(i,:),[],w(i));
    x(i+1 ,:) = xx(end,:);
    
end


figure(1)
plot(t,x(:,1));
figure(2)
plot(t,x(:,2));


y = x(:,1)+v;

for i=1:length(t)-1
    [ttt,xxxx] = ode45(@modeloxhat,[t(i) t(i+1)],xhat(i,:),[],y(i));
    xhat(i+1 ,:) = xxxx(end,:);
    

end

figure(4)
plot(t,xhat(:,1));
figure(5)
plot(t,xhat(:,2));

figure(6)
plot(t,x(:,1),t,xhat(:,1));
figure(7)
plot(t,x(:,2),t,xhat(:,2));



%%%  Modelos ------------------------


function dxdt = modeloPP(t,P)
omega=5;
zeta=0.2;

A=[0 1; -omega^2 -2*zeta*omega];
A11=0;
A12=1;
A21=-omega^2;
A22=-2*zeta*omega;

Rc=0.01;
Qc=4.47;

P11=P(1);
P12=P(2);
%P21=P(3);
P22=P(4);

dxdt(1,1) = 2*A11*P11+2*A12*P12-((P11^2)/Rc)+Qc;
dxdt(2,1) = P12*(A11+A22)+A12*P22+P11*A21-((P11*P12)/Rc)+Qc;
dxdt(3,1) = P12*(A11+A22)+A12*P22+P11*A21-((P11*P12)/Rc)+Qc;
dxdt(4,1) = 2*A21*P12+2*A22*P22-((P12^2)/Rc)+Qc;
end

function dxdt = modeloKalmanContinuo(t,x,w)

omega = 5;
zeta = 0.2;



dxdt(1,1)=x(2);
dxdt(2,1)=-omega^2*x(1)-2*zeta*omega*x(2)+w+12;

end

function dxdt = modeloxhat(t,xhat,y)

global k1 k2 zeta omega

dxdt(1,1)=xhat(2)+k1*(y-xhat(1));
dxdt(2,1)=-omega^2*xhat(1)-2*zeta*omega*xhat(2)+k2*(y-xhat(1));

end