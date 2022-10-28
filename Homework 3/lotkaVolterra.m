function Nf = lotkaVolterra(param, tf)
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

preyB = param(1);
preyD = param(2);
predatorB = param(3);
predatorD = param(4);
prey0 = param(5);
predator0 = param(6);

N0 = [prey0; predator0];

t0 = 0;
t1 = tf; 

opt = odeset('RelTol', 1e-7);
[tout, Nout] = ode45(@f, [t0; t1], N0, opt);

Nf = interp1(tout, Nout, tf); % use interpolation to find N's at tf

% NOTE: Comment out the following plotting code if using curve fitting

% function plot
figure(1);
plot(tout, Nout(:,1), 'r-', tout, Nout(:,2), 'b-');
xlabel('Time');
ylabel('Population');
legend('Prey', 'Predator');
set(gca, 'FontSize', 16);

% 3-D function plot
figure(2);
plot3(Nout(:, 1), Nout(:, 2), tout, 'r-');
xlabel('Prey population');
ylabel('Predator population');
zlabel('Time');
grid on;
set(gca, 'FontSize', 16);

% phase space plot
figure(3);
plot(Nout(:,1), Nout(:,2), 'k-');
xlabel('Prey population');
ylabel('Predator population');
grid on;
set(gca, 'FontSize', 16);

% Nested function to calculate dN/dt
function dNdt = f(t, N)
  prey = N(1, 1);
  predator = N(2, 1);
  dNdt(1, 1) = preyB .* prey - preyD .* prey .* predator;
  dNdt(2, 1) = predatorB .* prey .* predator - predatorD .* predator; 
end
end