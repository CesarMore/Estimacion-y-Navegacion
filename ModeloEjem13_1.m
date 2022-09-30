function [dxdt] = ModeloEjem13_1(t,x,ua,ub,q1,q2,q3)
%Parametros del sistema 
L = 0.003;
J = 0.002;
R = 2;
lambda = 0.1;
F = 0.001;
% q1 = 0;
% q2 = 0;
% q3 = 0;
%Representación en espacio de estado
%ia    = x1
%ib    = x2
%omega = x3
%theta = x4
dxdt(1,1) = -(R/L)*x(1) + ((x(3)*lambda)/L)*sin(x(4)) + (ua+q1)/L;
dxdt(2,1) = -(R/L)*x(2) - ((x(3)*lambda)/L)*cos(x(4)) + (ub+q2)/L;
dxdt(3,1) = -((3*lambda)/(2*L))*x(1)*sin(x(4)) + ((3*lambda)/(2*J))*x(2)*cos(x(4)) - ((F*x(3))/J) + q3;
dxdt(4,1) = x(3);
end


%%
% function [dxdt] = ModeloEjem13_1(t,x,ua,ub,q1,q2,q3)
% %Parametros del sistema 
% L = 0.003;
% J = 0.002;
% R = 2;
% lambda = 0.1;
% F = 0.001;
% % q1 = 0;
% % q2 = 0;
% % q3 = 0;
% %Representación en espacio de estado
% %ia    = x1
% %ib    = x2
% %omega = x3
% %theta = x4
% dxdt(1,1) = -(R/L)*x(1) + ((x(3)*lambda)/L)*sin(x(4)) + (ua+q1)/L;
% dxdt(2,1) = -(R/L)*x(2) - ((x(3)*lambda)/L)*cos(x(4)) + (ub+q2)/L;
% dxdt(3,1) = -((3*lambda)/(2*L))*x(1)*sin(x(4)) + ((3*lambda)/(2*J))*x(2)*cos(x(4)) - ((F*x(3))/J) + q3;
% dxdt(4,1) = x(3);
% end