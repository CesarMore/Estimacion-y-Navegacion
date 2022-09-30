
for iter = 1:Niter
    
    % Actualizar la mejor solución de cada partícula
    for i = 1:Npar                                          
                
        func_val = (Xpart(i) - 3.14)^2 + (Ypart(i) - 2.72).^2 + sin(3*Xpart(i) + 1.41) + sin(4*Ypart(i) - 1.73);
        
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
x = linspace(0,5,200);
y = linspace(0,5,200);
[X1,Y1] = meshgrid(x,y);
Z = (X1 - 3.14).^2 + (Y1 - 2.72).^2 + sin(3*X1 + 1.41) + sin(4*Y1 - 1.73);
contour3(X1,Y1,Z,100)
hold on
plot(Xb,Yb,'r--o')
plot(Xb(end),Yb(end),'r--o')

f3=figure(3)
clf(f3)
contour3(X1,Y1,Z,20)
hold on
plot(X(1,:),Y(1,:),'r--o',X(2,:),Y(2,:),'g--o')
