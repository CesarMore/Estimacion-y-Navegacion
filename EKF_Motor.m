function EKF_Motor

% Discrete-time extended Kalman filter simulation for two-phase
% step motor. Estimate the stator currents, and the rotor position
% and velocity, on the basis of noisy measurements of the stator
% currents.
Ra = 2; % Winding resistance
L = 0.003; % Winding inductance
lambda = 0.1; % Motor constant
J = 0.002; % Moment of inertia
F = 0.001; % Coefficient of viscous friction
ControlNoise = 0.001; % std dev of uncertainty in control inputs (amps)
MeasNoise = 0.1; % std dev of measurement noise (amps)
R = [MeasNoise^2 0; 0 MeasNoise^2]; % Measurement noise covariance
xdotNoise = [ControlNoise/L; ControlNoise/L; 0.5; 0];
% Define the continuous-time process noise covariance Q
Q = [xdotNoise(1)^2 0 0 0;
 0 xdotNoise(2)^2 0 0;
 0 0 xdotNoise(3)^2 0;
 0 0 0 xdotNoise(4)^2];

P = 1*eye(4); % Initial state estimation covariance
dt = 0.0002; % Simulation step size (seconds)
T = 0.001; % how often measurements are obtained
tf = 1.5; % Simulation length
t = [0:T:2];
x = [0; 0; 0; 0]; % Initial state
xhat = x; % Initial state estimate
w = 2 * pi; % Control input frequency (rad/sec)
Q = Q * T; % discrete-time process noise covariance
% Initialize arrays for plotting at the end of the program
NumTimeSteps = round(tf / T); % number of time steps
xArray = zeros(4, NumTimeSteps); % true state
xhatArray = zeros(4, NumTimeSteps); % estimated state
trPArray = zeros(1, NumTimeSteps); % trace of estimation error covariance
tArray = zeros(1, NumTimeSteps); % time array
i = 0; % loop index
% Begin simulation loop
for t = 0 : T : tf-T
 y = [x(1); x(2)] + MeasNoise*randn(2,1); % noisy measurement
 % Save data for plotting

 i = i + 1;
 xArray(:, i) = x;
 xhatArray(:, i) = xhat;
 trPArray(i) = trace(P);
 tArray(i) = t;

 % System simulation
 for tau = 0:dt:T-dt
 time = t + tau;
 ua = sin(w*time);
 ub = cos(w*time);
 xdot = [-Ra/L*x(1) + x(3)*lambda/L*sin(x(4)) + ua/L;
 -Ra/L*x(2) - x(3)*lambda/L*cos(x(4)) + ub/L;
 -3/2*lambda/J*x(1)*sin(x(4)) + ...
 3/2*lambda/J*x(2)*cos(x(4)) - F/J*x(3);
 x(3)];
 xdot = xdot + xdotNoise.*randn(4,1);
 x = x + xdot*dt; % rectangular integration
 x(4) = mod(x(4),2*pi); % keep the angle between 0 and 2*pi
 end

 % Compute the partial derivative matrices
 A = eye(4) + T * [-Ra/L, 0, lambda/L*sin(xhat(4)), xhat(3)*lambda/L*cos(xhat(4));
 0, -Ra/L, -lambda/L*cos(xhat(4)), xhat(3)*lambda/L*sin(xhat(4));
 -3/2*lambda/J*sin(xhat(4)), 3/2*lambda/J*cos(xhat(4)), -F/J, ...
 -3/2*lambda/J*(xhat(1)*cos(xhat(4))+xhat(2)*sin(xhat(4)));
 0 0 1 0];
 C = [1 0 0 0; 0 1 0 0];
 % Compute the Kalman gain
 K = P*C'*inv(C*P*C'+R);
 %K = P*C'*inv(R);
 % Update the state estimate
 ua = sin(w*t);
 ub = cos(w*t);
 deltax = [-Ra/L*xhat(1) + xhat(3)*lambda/L*sin(xhat(4)) + ua/L;
 -Ra/L*xhat(2) - xhat(3)*lambda/L*cos(xhat(4)) + ub/L;
 -3/2*lambda/J*xhat(1)*sin(xhat(4)) + ...
 3/2*lambda/J*xhat(2)*cos(xhat(4)) - F/J*xhat(3);
 xhat(3)] * T;
 xhat = xhat + deltax + K * (y - [xhat(1); xhat(2)]);
 % keep the angle estimate between 0 and 2*pi
 xhat(4) = mod(xhat(4), 2*pi);
 % Update the estimation error covariance.
 P = A * ((eye(4) - K * C) * P) * A' + Q;
end

% Gráfica de resultados 
close all;
figure; 

subplot(2,2,1); hold on; box on;
plot(tArray, xArray(1,:), 'b-', 'LineWidth', 2);
plot(tArray, xhatArray(1,:), 'r:', 'LineWidth', 2)
set(gca,'FontSize',12); ylabel('Current A (Amps)');

subplot(2,2,2); hold on; box on;
plot(tArray, xArray(2,:), 'b-', 'LineWidth', 2);
plot(tArray, xhatArray(2,:), 'r:', 'LineWidth', 2)
set(gca,'FontSize',12); ylabel('Current B (Amps)');

subplot(2,2,3); hold on; box on;
plot(tArray, xArray(3,:), 'b-', 'LineWidth', 2);
plot(tArray, xhatArray(3,:), 'r:', 'LineWidth', 2)
set(gca,'FontSize',12);
xlabel('Time (Seconds)'); ylabel('Speed (Rad/Sec)');

subplot(2,2,4); hold on; box on;
plot(tArray, xArray(4,:), 'b-', 'LineWidth', 2);
plot(tArray,xhatArray(4,:), 'r:', 'LineWidth', 2)
set(gca,'FontSize',12);
xlabel('Time (Seconds)'); ylabel('Position (Rad)');
legend('True', 'Estimated');



