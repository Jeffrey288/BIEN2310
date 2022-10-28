function Nf = lotkaVolterraLogistic_HartantoKwee_Jeffrey(param, tf)
% LOTKAVOLTERRA. Solve the Lotka-Volterra predator-prey population dynamics 
% problem.
% Usage: Nf = lotkaVolterra(param, tf), where
% - ttest is the time at which the population estimates are desired 
% - Nf is a row vector of [N_prey N_predator] at time = tf
% - param is a column vector of parameters [a; b; c; d; N0_prey; N0_predator]
%       a = birth rate of prey
%       b = death rate of prey per predator
%       c = birth rate of predator per prey
%       d = death rate of predator
% Unpack the parameters
% lotkaVolterra([0.5, 0.01, 0.005, 0.4, 50, 80], 200)
% lotkaVolterraLogistic_HartantoKwee_Jeffrey([0.5, 0.01, 0.005, 0.4, 5000, 50, 80], 2000)

preyB = param(1);
preyD = param(2);
predatorB = param(3);
predatorD = param(4);
preyK = param(5);
prey0 = param(6);
predator0 = param(7);

N0 = [prey0; predator0];

t0 = 0;
t1 = tf; 

opt = odeset('RelTol', 1e-90, 'Events', @reachSS);
[tout, Nout] = ode45(@f, [t0; t1], N0, opt);

if (tout(end) < tf)
    Nf = Nout(end)
else
    Nf = interp1(tout, Nout, tf); % use interpolation to find N's at tf
end

% NOTE: Comment out the following plotting code if using curve fitting

% function plot
figure(1);
plot(tout, Nout(:,1), 'r-', tout, Nout(:,2), 'b-');
title('Hartanto Kwee Jeffrey, function plot');
xlabel('Time');
ylabel('Population');
legend('Prey', 'Predator');
set(gca, 'FontSize', 16);


% 3-D function plot
figure(2);
plot3(Nout(:, 1), Nout(:, 2), tout, 'r-');
title("Hartanto Kwee Jeffrey, 3D function plot");
xlabel('Prey population');
ylabel('Predator population');
zlabel('Time');
grid on;
set(gca, 'FontSize', 16);


% phase space plot
figure(3);
plot(Nout(:,1), Nout(:,2), 'k-');
title("Hartanto Kwee Jeffrey, phase space plot");
xlabel('Prey population');
ylabel('Predator population');
grid on;
set(gca, 'FontSize', 16);


% Nested function to calculate dN/dt
function dNdt = f(t, N)
  prey = N(1, 1);
  predator = N(2, 1);
  dNdt(1, 1) = preyB .* prey .* (1 - prey ./ preyK) - preyD .* prey .* predator;
  dNdt(2, 1) = predatorB .* prey .* predator - predatorD .* predator; 
end

function [position, isterminal, direction] = reachSS(t, N)
dNdt = f(t, N);
R = N(1, 1);
F = N(2, 1);
dRdt = dNdt(1, 1);
dFdt = dNdt(2, 1);

% checks relative magnitude
condition = (abs(dRdt ./ R) < 0.00001 && abs(dFdt ./ F) < 0.00001);

position = condition - 0.5;
isterminal = 1;
direction = 1;
end

end