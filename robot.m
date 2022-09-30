function [Xd,u]=robot(t,X)

x1 = X(1);
x2 = X(2);
x3 = X(3);
x4 = X(4);

x = [x1;x2;x3;x4];

xh1 = X(5);
xh2 = X(6);
xh3 = X(7);
xh4 = X(8);

xh = [xh1;xh2;xh3;xh4];

mgl = 10;
I = 1;
J = 1;
k = 100;

a = mgl/I;
b = k/I;
c = k/J;
d = 1/J;

z1 = x1;
z2 = x2;
z3 = -a*sin(x1)-b*(x1-x3);
z4 = -a*x2*cos(x1)-b*(x2-x4);

z = [z1;z2;z3;z4];

k1 = 4788;
k2 = 2318.5;
k3 = 419;
k4 = 33.5;

K = [k1 k2 k3 k4];
ref = 2*pi;

%v = -K*z;
v = -k1*(z1-ref)-k2*z2-k3*z3-k4*z4;

u = (I*J)/(k)*(v-a*sin(x1)*(x2^2+a*cos(x1)+b)-b*(x1-x3)*(b+c+a*cos(x1)));



A = [0 1 0 0;-b 0 b 0;0 0 0 1;c 0 -c 0];

psi = [0;-a*sin(x1);0;d*u ];

C = [1 0 0 0];

xd = A*x+psi;

y = C*x;


%place 5 veces mas rapido con C

H = [46;591;14.3;-419];

Xhd = A*xh +psi+H*(y-C*xh);

Xd = [xd;Xhd];



end
