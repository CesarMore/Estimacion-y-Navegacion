function dxdt = modeloxhat(t,xhat,y,P11) 
omega=5; 
zeta=0.2; 

dxdt(1,1) = 0*xhat(1) + xhat(2) + (P11/0.01)*(y - xhat(1));
dxdt(2,1) = -omega^2*xhat(1) -2*zeta*omega*xhat(2);
end
