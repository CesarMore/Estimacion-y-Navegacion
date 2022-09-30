function [dxdt] = modeloFK(t,x,w)
    omega = 5;
    zeta  = 0.2;
    
    dxdt(1,1) = x(2);
    dxdt(2,1) = -omega^2*x(1) - 2*zeta*omega*x(2) + w + 12;
end

