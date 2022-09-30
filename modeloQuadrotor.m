function[dxdt] = modeloQuadrotor(t,x,taux,tauy,tauz,Omega1,Omega2,Omega3,Omega4)
Jx=0.0552;
Jy=0.0552;
Jz=0.1104;
Jm = 0.56;
Kafx = 0.008;
Kafy = 0.008;
Kafz = 0.009;
Omega1 = 0;
Omega2 = 0;
Omega3 = 0;
Omega4 = 0;
Omega = Omega1-Omega2+Omega3+Omega4;

dxdt(1,1) = x(2);
dxdt(2,1) = (Jy-Jz)*x(4)*x(6)/Jx + Jm*Omega/Jx + taux/Jx - Kafx*x(2)/Jx;
dxdt(3,1) = x(4);
dxdt(4,1) = (Jz-Jx)*x(2)*x(6)/Jy + Jm*Omega/Jy + tauy/Jy - Kafy*x(4)/Jy;
dxdt(5,1) = x(6); 
dxdt(6,1) = (Jx-Jy)*x(2)*x(4)/Jz + tauz/Jx - Kafz*x(6)/Jz;
end