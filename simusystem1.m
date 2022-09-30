clc
clear all
close all

a       = 0.7;
yr      = 2;
t       = [0:0.01:2];

x(1,1)  = 0;
y(1)    = 0;
u(1)    = 0;
k       = [6 20 200];

for j = 1 : length(k)
    ey(j,1)   = yr - y(1);
    u(j,1)    = k(j)*ey(1) + a*yr;
    for i = 1 : length(t)-1
        if( t(i) > 1 )
            d = 2;
        else
            d = 0;
        end
        [tt, xx]  = ode45(@modelsystem1, [t(i) t(i+1)], x(j, i), [], u(j,i), d);
        y(i+1)    = xx(end);
        
        x(j,i+1)  = y(i+1);
        ey(j,i+1) = yr - y(i+1);
        u(j,i+1)  = k(j)*ey(j,i+1) + a*yr;
    end 
    x(j+1,1) = 0;
end

figure(1)
plot(t,x(1,:))
hold on
plot(t,x(2,:))
hold on
plot(t,x(3,:))
hold on
grid on
legend('k=6','k=20', 'k=200')
xlabel('Tiempo [seg]')
ylabel('Salida, y')

figure(2)
plot(t,u(1,:))
hold on
plot(t,u(2,:))
hold on
plot(t,u(3,:))
hold on
grid on
legend('k=6','k=20', 'k=200')
xlabel('Tiempo [seg]')
ylabel('Entrada de control, u')

figure(3)
plot(t,ey(1,:))
hold on
plot(t,ey(2,:))
hold on
plot(t,ey(3,:))
hold on
grid on
legend('k=6','k=20', 'k=200')
xlabel('Tiempo [seg]')
ylabel('Error de seguimiento, e_y')
