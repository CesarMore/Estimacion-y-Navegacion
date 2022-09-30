% clc 
% clear all
% close all 
% 
% T = 0.01;
% t = [0:T:10];
% 
% xk = [8 7];
% 
% yk = 0;
% vk = 0;
% Pk = 0;
% K  = 0;
% x = sin(2*pi*t);
% 
% vk = wgn(1,1,-2);
% 
% for k=2:length(t)
%     I = eye(2);
%     P0 = eye(2); 
%     Rk = 0.01;
%     Hk = [1 0.99^(k-1)];
%     K = P0*Hk'*inv(Hk*P0*Hk' + Rk);
%     yk = xk(1) + (0.99^(k-1))*xk(2) + vk;
%     xk = xk(k-1) + K*(yk - Hk*xk(k-1));
%     Pk = (I - K*Hk)*Pk(k-1);
% 
% end

clear all
close all
clc


Rk = [0.01];


xh = [8;7];
Pk = [1 0;0 1];

t = [];
m = [-2];
n = [-1];

for k=2:100

    Hk = [1 0.99^(k-1)];

    m = m+1;
    n = n+1;

    
    Kk = Pk(1:2,k+m:k+n)*transpose(Hk)*inv(Hk*Pk(1:2,k+m:k+n)*transpose(Hk)+Rk);
end

