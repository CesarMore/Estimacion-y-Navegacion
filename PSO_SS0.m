
clc
clear all 

Niter = 100;   
Npar = 5*20;   
c1 = 2;        
c2 = 2;      

w = 1;
wd = 10;
wt = 10;
Xpart = rand(Npar,1)-1;
Ypart = rand(Npar,1)-1;
Vxpart = zeros(Npar,1);
Vypart = zeros(Npar,1);
Pbest = 100*ones(Npar,1);
BestX = zeros(Npar,1);
BestY = zeros(Npar,1);
X = Xpart;
Y = Ypart;
Gx=23;
Gy=30;
Ox=0;
Oy=0;
Xops = 1.5*rand(1,50)*23;
Yops = 1.5*rand(1,50)*30;
% r1 = 1.1*randi([0 1]);
% r2 = 1.1*randi([0 1]);
r1 = 0.1;
r2 = 0.9;

for iter = 1:Niter
    
    % Actualizar la mejor solución de cada partícula
    for i = 2:Npar                                          
        %Funcion a minimizar
        f1 = sqrt(abs((Xpart(i)-Gx)^2 + (Ypart(i)-Gy)^2));
        d = sqrt((Xpart(i)-Xpart(i-1))^2 + (Ypart(i)-Ypart(i-1))^2);
        v = Vxpart(i) + Vypart(i);
        t = d/v;
        f2 = 1/(wd*d + wt*t);
  
        func_val= 0.9*f1 + 0.1*f2;
        
       if func_val < Pbest(i)
            for j=1:50
                if (Xpart(i)~=Xops(j)) && (Ypart(i)~=Yops(j))
                    BestX(i) = Xpart(i);
                    BestY(i) = Ypart(i);
                    Pbest(i) = func_val;
                end
            end
        end
    end

    %Buscar la mejor solución de todas las partículas
    [gbest,posbest] = min(Pbest);
    Xb(iter) = Xpart(posbest);
    Yb(iter) = Ypart(posbest);
    P(iter) = gbest;
         
    % Actualizar los vectores de posición y velocidad
    for i = 1:Npar  
    %+ c1*r1*(BestX(i) - Xpart(i))
    %+ c1*r1*(BestY(i) - Ypart(i))
        Vxpart(i) = rand*w*Vxpart(i) + c1*rand*(BestX(i) - Xpart(i)) + c2*rand*(Xb(iter) - Xpart(i));
        Vypart(i) = rand*w*Vypart(i) + c1*rand*(BestY(i) - Ypart(i)) + c2*rand*(Yb(iter) - Ypart(i));
        
        Xpart(i) = Xpart(i) + Vxpart(i);
        Ypart(i) = Ypart(i) + Vypart(i);        
    end   

    X(:,iter) = Xpart;
    Y(:,iter) = Ypart;
end
%Grafica de la mejor particula
% figure(1)
% plot(P)

f2=figure(2)
clf(f2)
hold on
plot(0,0,'o','MarkerFaceColor','b')
plot(23,30,'o','MarkerFaceColor','b')
hold on
plot(Xb,Yb,'-s','MarkerFaceColor','r')
%hold on
%plot(Xb(end),Yb(end),'-s','MarkerFaceColor','r')
hold on
plot(Xops(:),Yops(:),'o','MarkerFaceColor','k')
grid minor
ylim([0 35])
xlim([0 25])
%