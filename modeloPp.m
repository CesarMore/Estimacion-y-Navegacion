function dxdt = modeloPp(t,P)
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
P21=P(3);
P22=P(4);

dxdt(1,1) = 2*A11*P11+2*A12*P12-((P11^2)/Rc)+Qc;
dxdt(2,1) = P12*(A11+A22)+A12*P22+P11*A21-((P11*P12)/Rc)+Qc;
dxdt(3,1) = P12*(A11+A22)+A12*P22+P11*A21-((P11*P12)/Rc)+Qc;
dxdt(4,1) = 2*A21*P12+2*A22*P22-((P12^2)/Rc)+Qc;