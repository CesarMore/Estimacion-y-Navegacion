clc
clear all 
close all
syms P11 P12 P22 K1 K2 T sigma R

eq1 = K1*P11        - 2*T*P12       + T^2*P22 == (T^4*sigma^2)/4
eq2 = 0*P11         + K1*P12        - T*P22   == -(T^3*sigma^2)/2
eq3 = 0*P11         + K2*P12        - 0*P22   == T^2*sigma^2
%eq4 = -P11/(P11+R)  + 0*P12         + 0*P22   + K1   + 0*K2 == 0
%eq5 = 0*P11         - P12/(P11+R)   + 0*P22   + 0*K1 + K2   == 0
[A B] = equationsToMatrix([eq1 eq2 eq3],[P11 P12 P22])

R = linsolve(A,B)
% 
% eq4 = K1   + 0*K2  == P11/(P11+R)
% eq5 = 0*K1 + K2    == P12/(P11+R)
% 
% [C D] = equationsToMatrix([eq4 eq5],[K1 K2])
% Res = linsolve(C,D)


