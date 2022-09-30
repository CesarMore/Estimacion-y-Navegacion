clc
clear all
close all

a       = 0.7;
yr      = 2;
t       = [0: 0.01 : 50];

x(1, 1) = 0;
y(1)    = 0;
u(1)    = 0;

k0 = 6;
k1 = 15;

k       = [k0 k1];
d       = 0;
ey(1)   = yr - y(1);
intey(1)  = 0;

for i = 1 : length(t)-1
    if( t(i) > 10) 
        d = 2;
    end

    if( t(i) > 20 )
        d = 2 + 2*sin(t(i));
    end

    [tt, xx] = ode45(@modelsystem1, [t(i) t(i+1)], x(i), [], u(i), d);
    y(i+1)   = xx(end);
    
    x(i+1)   = y(i+1);
    ey(i+1)  = yr - y(i+1);

    %Integración de error de seguiento
    intey(i+1) = intey(i)+(t(i+1)-t(i))*ey(i+1);

    u(i+1)   = k(2)*ey(i+1) + k(1)*intey(i+1) + a*yr;
end 

figure(1)
plot(t,x)
grid on
xlabel("Tiempo, [Seg]")
ylabel("Salida, y")

figure(2)
plot(t,u)
grid on
xlabel('Tiempo [seg]')
ylabel('Entrada de control, u')


