clc
clear all
close all

T = 0.01;
t = [0:T:20];

VARwp = 0.1;
VARw = 2.5;   %Puede cambiar

w = sqrt(0.25*(VARw))*((2*randn(length(t),1)-1));
wp = sqrt(0.25*(VARwp))*((2*randn(length(t),1)-1));

VARv = 0.009;
v = sqrt(0.25*(VARv))*((2*randn(length(t),1)-1));

x3real = -4;

% mean(w)
% mean(wp)
% var(w)
% var(wp)



x = [0 0 0];

for i=1:length(t)-1
    [tt,xx] = ode45(@Modelo13_4,[t(i) t(i+1)],x(i,:),[],w(i),wp(i));
    x(i+1 ,:) = xx(end,:);
    
end

figure(1)
plot(t,x(:,1))
title('x1')

figure(2)
plot(t,x(:,2))
title('x2')

figure(3)
plot(t,x(:,3))
title('x3')


%-----Filtro de Kalman param---------

R = [var(v) 0;0 var(v)];

H = [1 0 0;
    0 1 0];

Q = [var(w) 0;0 var(wp)];

xhat = [mean(x(1,1)) mean(x(1,2)) -1];   %Se puede mover

P = [0 0 0;
    0 0 0;
    0 0 20];


% ex = x(1,:)'-xhat';
% ex1 = (ex)*transpose(ex);
% P = [mean(ex1(1,1)) mean(ex1(1,2)) mean(ex1(1,3));
%     mean(ex1(2,1)) mean(ex1(2,2)) mean(ex1(2,3));
%     mean(ex1(3,1)) mean(ex1(3,2)) 20];


b = -0.4;

y = [x(:,1)+v x(:,2)+v];


for i = 1:length(t)-1

F = [0 1 0;
    xhat(i,3) b xhat(i,1);
    0 0 0];

L = [0 0;
    -xhat(i,3) 0;
    0 1];

% L = [0 0;
%     -xhat(i,3);
%     0];

f0 = [xhat(i,2);
    xhat(i,3)*xhat(i,1) + b*xhat(i,2) - xhat(i,3);
    0];

K = P(:,:,i)*H'*inv(R);

Ppunto = F*P(:,:,i) + P(:,:,i)*F' + L*Q*L' - P(:,:,i)*H'*inv(R)*H*P(:,:,i);

P(:,:,i+1) = P(:,:,i) + T*Ppunto;

xhatpunto = f0 + K*(y(i,:)' - H*xhat(i,:)');

xhat(i+1,:) = xhat(i,:) + T*xhatpunto';




end


figure(4)
plot(t,x(:,1),t,xhat(:,1))
title('x1 xhat 1')

figure(5)
plot(t,x(:,2),t,xhat(:,2))
title('x2  xhat 2')

figure(6)
plot(t,x(:,3),t,xhat(:,3))
title('x3  xhat 3')


errorx1 = x(:,1) - xhat(:,1);
errorx2 = x(:,2) - xhat(:,2);
errorx3 = x(:,3) - xhat(:,3);


figure(7)
plot(t,errorx1)
title('ex1')

figure(8)
plot(t,errorx2)
title('ex2')

figure(9)
plot(t,errorx3)
title('ex3')


errorx3real = x3real - xhat(:,3);

figure(10)
plot(t,errorx3real)
title('ex3real')




