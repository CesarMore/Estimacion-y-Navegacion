clc 
clear 
close all
T = 0.01;             %Intervalo de muestreo 
t = [0:T:30];         % Tiempo inicial y final  con 0.1 entre tiempos       
x = [0 0 0 0 0 0 ];   % vector variables de estado(valores iniciales)                             

%Parámetros
l=0.1969;
Kax=0.008;
Kay=0.008;
Jz=0.1104; 
Kaz=0.009;
Jx=0.0552;
Jy=0.0552;
ku=0.1188; 
ks=0.0036;
Jm=0.1;
%-------------------
tau1 =0;               %Entrada de control  para roll
tau2 =0;               %Entrada de control  para pitch
tau3 =0;               %Entrada de control  para yaw

%--------------------------------Ganancias---------------------------------
% se pueden modificar
%theta
kd2 = 50;
kp2 = 80;
%-----------
%roll
kd1 = 50;
kp1 = 100;
%-----------
%yaw
kd3 = 50;
kp3 = 100;
%Valores iniciales para omega
omega1=0;
omega2=0;
omega3=0;
omega4=0;
omega=0;
OMEGA=0;
              
%--------------------------------------------------------------------------
%------------------------- REFERENCIA  deseada ---------------------------
r0  = 0.5*cos(t);           
r0p = -0.5*sin(t);
%--------------------------------------------------------------------------
for i=1:length(t)-1            
    [tt, xx] = ode45(@modeloQuad2,[t(i) t(i+1)],x(i,:),[],tau1(i),tau2(i),tau3(i),OMEGA(i));
    %[tt, xx] = ode45(@modeloQuad2,[t(i) t(i+1)],x(i,:),[],tau1(i),tau2(i),tau3(i));
    x(i+1,:) = xx(end,:);
    
    %----------------------------------------------------------------------
    theta = x(i+1,1);
    thetap = x(i+1,2);
    phi = x(i+1,3);
    phip = x(i+1,4);
    psi = x(i+1,5);
    psip = x(i+1,6);   
    %----------------------- ley de control para pitch --------------------
    et  = r0(i+1)  -x(i+1,1);
    etp = r0p(i+1) -x(i+1,2);
    ephi = r0(i+1) -x(i+1,5);
    ephip = r0p(i+1) -x(i+1,6);
    epsi = r0(i+1)  -x(i+1,3);
    epsip= r0p(i+1) -x(i+1,4);


    tau2(i+1) = Jy*(-(Jz-Jx)*psi*phi/Jy-Jm*OMEGA(i)+Kay*theta/Jy+kp2*et+kd2*etp);
    %tau2(i+1) = Jy*((Kay*theta - (Jz-Jx)*phi*psi - Jm*OMEGA(i)) + kp2*et+kd2*etp);
    %----------------------- ley de control para roll --------------------
    tau1(i+1) = Jx*(-(Jy-Jz)*theta*psi/Jx-Jm*OMEGA(i)+Kax*phi/Jx+kp1*epsi+kd1*epsip);
    %tau1(i+1) = Jx*(Jx*(Kax*phi - (Jy-Jz)*theta*psi - Jm*OMEGA(i)) + kp1*ephi+kd1*ephip);
    %----------------------- ley de control para yaw --------------------
    tau3(i+1) = Jz*(-(Jx-Jy)*phi*theta/Jz-Jm*OMEGA(i)+Kaz*psi/Jz+kp3*ephi+kd3*ephip);
    %tau3(i+1) = Jz*((Kaz*psi - (Jx-Jy)*phi*psi) + kp3*epsi+kd3*epsip);

    %Calcular omega
    tau=[tau1(i+1);
         tau2(i+1);
         tau3(i+1)];
    L=0.197; %rotor arm length
    kt=0.1188; %lift coefficient
    kl=0.01; %torque coefficient
    
    U=[ 0    -L*kt   0     L*kt;
       L*kt   0    -L*kt   0;
       kl   -kl     kl    -kl];
   
    uinv=pinv(U); %Valor constante

   omegacuad=uinv*tau; %Genera imaginarios
   omega=real(sqrt(omegacuad));
   omega1=omega(1);
   omega2=omega(2);
   omega3=omega(3);
   omega4=omega(4);
   OMEGA(i+1)=omega1-omega2+omega3-omega4;

end

% graficas de orientacion
figure(1)
subplot(311)
plot(t,r0,t,x(:,1))
legend('referencia','\phi')
ylabel('Orientacion roll[rad]')
subplot(312)
plot(t,r0,t,x(:,3))
legend('referencia','\theta')
ylabel('Orientacion pitch[rad]')
subplot(313)
plot(t,r0,t,x(:,5))
legend('referencia','\psi')
ylabel('Orientacion yaw[rad]')
xlabel('Tiempo[segundos]')




