
ti = 0:0.01:60;

A = -4;
w = 0.251;
r0 = A*square(w*ti);
function y = r0(x)

rop = diff(y,1);
% plot(ti,r0) 