function dxdt = modelsystem1(t,x,u,d)
    a = 0.7;
    %Representacion de estado
    dxdt = -a*x + u + d; 
end