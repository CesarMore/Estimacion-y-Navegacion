clc 
clear all 
close all 

global K C B A u zeta omega k1 k2
T = 0.01;
t = [0:T:10];
x = [0 0];
xhat = [0 0];
VAR = 4.47;
VAR1 = 0.0001;
w=sqrt(0.25*(VAR))*((2*randn(length(t),1)-1));
v=sqrt(0.25*(VAR1))*((2*randn(length(t),1)-1));
omega = 5;
zeta  = 0.2;
A = [0 1;-omega^2 -2*zeta*omega];
B = [0;0];
C = [1 0];
Q = var(w);
R = var(v);
Qc = Q/T;
Rc = R/T;
Pp = [0 0; 0 0];
P  = [2 0; 0 2];
%-------- K
K = P*transpose(C)*inv(Rc);
k1=K(1);
k2=K(2);
u = 0;

for i=1:length(t)-1
    [tt,xx] = ode45(@modeloFK,[t(i) t(i+1)],x(i,:),[],w(i));
    x(i+1,:) = xx(end,:);   

    x1(:,i+1) = x(i+1,1);
    x2(:,i+1) = x(i+1,2); 
end
y = x(:,1)+v;

for i=1:length(i)-1
    [tt, PP] = ode45(@modeloPp,[0 T],P,[]);
    z(i,:) = PP(end,:);

    P=[z(1,1) z(1,2); z(1,3) z(1,4)];
end

for i=1:length(t)-1
    [tt,xx] = ode45(@modeloxhat,[t(i) t(i+1)],xhat(i,:),[],y(i));
    xhat(i+1,:) = xx(end,:);   
    xhat1(:,i+1) = xhat(i+1,1);
    xhat2(:,i+1) = xhat(i+1,2);   
end
%
figure(1)
plot(t,x1,t,xhat1);
legend('x1 (real)','xhat1 (estimada)')
figure(2)
plot(t,x2,t,xhat1);
legend('x2 (real)','xhat2 (estimada)')
% figure(3)
% plot(t,P(:,1),t,P(:,2),t,P(:,4))
% legend('P11','P12','P22')
% figure(4)
% plot(t,xhat(:,1))
%
function dxdt = modeloPp(t,P)
omega=5; zeta=0.2; A=[0 1; -omega^2 -2*zeta*omega]; A11=0; A12=1;
A21=-omega^2; A22=-2*zeta*omega;
P11=P(1); P12=P(2); P21=P(3); P22=P(4);
Rc=0.01;
Qc=4.47;
dxdt(1,1) = 2*A11*P11+2*A12*P12-((P11^2)/Rc)+Qc;
dxdt(2,1) = P12*(A11+A22)+A12*P22+P11*A21-((P11*P12)/Rc)+Qc;
dxdt(3,1) = P12*(A11+A22)+A12*P22+P11*A21-((P11*P12)/Rc)+Qc;
dxdt(4,1) = 2*A21*P12+2*A22*P22-((P12^2)/Rc)+Qc;
end

function dxdt = modeloxhat(t,xhat,y)
global k1 k2 zeta omega
dxdt(1,1)=xhat(2)+k1*(y-xhat(1));
dxdt(2,1)=-omega^2*xhat(1)-2*zeta*omega*xhat(2)+k2*(y-xhat(1));
end

