clc
clear all


Niter = 50;    
Npar = 200;     
c1 = 2;       
c2 = 2;        
wmin = 0.8;    
wmax= 1.2;
Xpart = 1*rand(Npar,1)*30;
Ypart = 1*rand(Npar,1)*30;
Vxpart = zeros(Npar,1);
Vypart = zeros(Npar,1);
Pbest = 1000*ones(Npar,1);
BestX = zeros(Npar,1);
BestY = zeros(Npar,1);
X = Xpart;
Y = Ypart;
Gx = 170;
Gy = 170;
Ox = 100;
Oy = 100;
Xops = rand(1,8)*170;
Yops = rand(1,8)*170;

for iter = 1:Niter
    
    % Actualizar la mejor solución de cada partícula
    for i = 1:Npar                                          
        %Funcion a minimizar
        f1= sqrt(abs((Xpart(i)-Gx)^2 + (Ypart(i)-Gy)^2));
        %Funcion para evadir objeto
        f2 = sqrt(abs((Xpart(i)-Ox)^2 + (Ypart(i)-Oy)^2));
%         f3 = sqrt(abs((Ox1)^2 + (Oy1)^2));
%         f4 = sqrt(abs((Ox2)^2 + (Oy2)^2));
%         f5 = sqrt(abs((Ox3)^2 + (Oy3)^2));
        func_val = f1 + f2;

        if func_val < Pbest(i)
            BestX(i) = Xpart(i);
            BestY(i) = Ypart(i);
            Pbest(i) = func_val;
        end
    end

    %Buscar la mejor solución de todas las partículas
    [gbest,posbest] = min(Pbest);
    Xb(iter) = Xpart(posbest);
    Yb(iter) = Ypart(posbest);
    P(iter) = gbest;
         
    % Actualizar los vectores de posición y velocidad
    for i = 1:Npar  
        w=wmax-((wmax-wmin)*i)/Npar;
        Vxpart(i) = rand*w*Vxpart(i) + c1*rand*(BestX(i) - Xpart(i)) +c2*rand*( Xb(iter) - Xpart(i));
        Vypart(i) = rand*w*Vypart(i) + c1*rand*(BestY(i) - Ypart(i)) +c2*rand*( Yb(iter) - Ypart(i));
        
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
plot(25,25,'o','MarkerFaceColor','b')
plot(170,170,'o','MarkerFaceColor','b')
hold on
plot(Xb,Yb,'-s','MarkerFaceColor','r')
%hold on 
%plot(Xb(end),Yb(end),'-s','MarkerFaceColor','r')
hold on
plot(Xops(:),Yops(:),'d','MarkerFaceColor','k')
grid minor
ylim([0 200])
xlim([0 200])





% rectangle('Position',position,'Curvature',[1 1],'FaceColor','blue','EdgeColor','red','LineWidth',1,'LineStyle','-.')
% hold on
% plot(30,30,'g--o')
% hold on
% position = [90 90 3 3]; 
% rectangle('Position',position,'Curvature',[1 1],'FaceColor','blue','EdgeColor','red','LineWidth',1,'LineStyle','-.')
% hold on
% position = [20 20 3 3]; 
% rectangle('Position',position,'Curvature',[1 1],'FaceColor','blue','EdgeColor','red','LineWidth',1,'LineStyle','-.')
% hold on
% position = [150 150 3 3]; 
% rectangle('Position',position,'Curvature',[1 1],'FaceColor','blue','EdgeColor','red','LineWidth',1,'LineStyle','-.')
% hold on
% position = [120 160 3 3]; 
% rectangle('Position',position,'Curvature',[1 1],'FaceColor','blue','EdgeColor','red','LineWidth',1,'LineStyle','-.')
% hold on
