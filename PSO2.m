clc
clear all

Niter = 50;
Npar = 5;
c1 = 2;
c2 = 2;
w = 1;
Xpart = 2*rand(Npar,1)*10-5;
Ypart = 1*rand(Npar,1)*10-5;
Vxpart = zeros(Npar,1);
Vypart = zeros(Npar,1);
Pbest = 1000*ones(Npar,1);
BestX = zeros(Npar,1);
BestY = zeros(Npar,1);
X = Xpart;
Y = Ypart;

for iter = 1:Niter
    
    % Actualizar la mejor solución de cada partícula
    for i = 1:Npar                                          
        
        func_val = Xpart(i)^2 + Ypart(i)^2;
        %func_val = (Xpart(i) - 3.14)^2 + (Ypart(i) - 2.72).^2 + sin(3*Xpart(i) + 1.41) + sin(4*Ypart(i) - 1.73);
        
        if func_val < Pbest(i)
            BestX(i) = Xpart(i);
            BestY(i) = Ypart(i);
            Pbest(i) = func_val;
        end
    end

    %Buscar la mejor solución de todas las partículas
    [gbest posbest] = min(Pbest);
    Xb(iter) = Xpart(posbest);
    Yb(iter) = Ypart(posbest);
    P(iter) = gbest;
         
    % Actualizar los vectores de posición y velocidad
    for i = 1:Npar                
        Vxpart(i) = rand*w*Vxpart(i) + c1*rand*(BestX(i) - Xpart(i)) +c2*rand*( Xb(iter) - Xpart(i));
        Vypart(i) = rand*w*Vypart(i) + c1*rand*(BestY(i) - Ypart(i)) +c2*rand*( Yb(iter) - Ypart(i));
        
        Xpart(i) = Xpart(i) + Vxpart(i);
        Ypart(i) = Ypart(i) + Vypart(i);        
    end   

    X(:,iter) = Xpart;
    Y(:,iter) = Ypart;
end

figure(1)
plot(P)

f2=figure(2)
clf(f2)
x = linspace(-5,5);
y = linspace(-5,5);
[X1,Y1] = meshgrid(x,y);
Z = X1.^2 + Y1.^2;
contour3(X1,Y1,Z,100)
hold on
plot(Xb,Yb,'r--o')
%

f3=figure(3)
clf(f3)
contour3(X1,Y1,Z,20)
hold on
plot(X(1,:),Y(1,:),'r--o',X(2,:),Y(2,:),'g--o')

