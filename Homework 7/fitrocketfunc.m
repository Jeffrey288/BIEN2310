function [] = fitrocketfunc(tdata, ydata1, ydata2)
% FITROCKET - Nonlinear regression to find the parameters of the rocket
% problem
% Usage: [] = fitrocket(tdata, ydata1, ydata2), where 
% - (tdata, ydata1) is the first set of data points, generated from one
%   flight of the rocket at one initial velocity, and (tdata, ydata2) is
%   the second set of data points, for another initial velocity.
%
% fitrocketfunc(hw7_q2_data.t, hw7_q2_data.y1, hw7_q2_data.y2)

% first, we try to fit only the first data-set 

beta0 = [0.01; 100]; % initial guess

% if we want to fit g too, use the following instead
% beta0 = [0.01; 100; 10]; % initial guess

modelfunc = @(params, t) rocketfunc(t, params, false);

[lsbeta, residual, jacobian, CovB, mse] = nlinfit(tdata, ydata1, modelfunc, beta0);
% lsbeta
% residual

ci = nlparci(lsbeta, residual, 'Jacobian', jacobian);

figure(1);
hold off;
rocketfunc(max(tdata), lsbeta, 1);
hold on;
plot(tdata, ydata1, 'rx');


% The following will plot the prediction intervals
% [Ypred delta] = nlpredci(modelfunc, tdata, lsbeta, residual, 'Jacobian', jacobian, 'PredOpt', 'observation');
% plot(tdata, Ypred + delta, 'r:', tdata, Ypred - delta, 'r:');


ti = 'Rocket altitudes';
Dstr1 = sprintf('$$D=%.4g\\ (CI:\\ [%.4g, %.4g])$$', lsbeta(1), ci(1, 1), ci(1, 2));
v0str = sprintf('$$v_0=%.4g\\ (CI:\\ [%.4g, %.4g])$$', lsbeta(2), ci(2, 1), ci(2, 2));
% gstr = sprintf('$$g=%.4g\\ (CI:\\ [%.4g, %.4g])$$', lsbeta(3), ci(3, 1), ci(3, 2));
title({ti, Dstr1, v0str}, 'Interpreter', 'latex');
% title({ti, Dstr1, v0str, gstr}, 'Interpreter', 'latex');


hold off;

% next, we try to fit two 2 data-sets with different initial velocities,
% but presumably the same D
 
beta0 = [0.01; 100; 200];
modelfunc = @(params, t) rocketfunc2(t, params, false);
 
[lsbeta2, resnorm, residual, exitflag, output, lambda, jacobian] = lsqcurvefit(modelfunc, beta0, tdata, [ydata1, ydata2]);

ci2 = nlparci(lsbeta2, residual, 'Jacobian', jacobian);

figure(2);
hold off;
rocketfunc2(max(tdata), lsbeta2, 1);
hold on;
plot(tdata, ydata1, 'rx', tdata, ydata2, 'bx');

ti = 'Rocket altitudes';
Dstr2 = sprintf('$$D=%.4g\\ (CI:\\ [%.4g, %.4g])$$', lsbeta2(1), ci2(1, 1), ci2(1, 2));
v0_1str = sprintf('$$v_{0,1}=%.4g\\ (CI:\\ [%.4g, %.4g])$$', lsbeta2(2), ci2(2, 1), ci2(2, 2));
v0_2str = sprintf('$$v_{0,2}=%.4g\\ (CI:\\ [%.4g, %.4g])$$', lsbeta2(3), ci2(3, 1), ci2(3, 2));
title({ti, Dstr2, v0_1str, v0_2str}, 'Interpreter', 'latex');

% The following will plot the predictor intervals
% [Ypred2 delta2] = nlpredci(modelfunc, tdata, lsbeta2, residual, 'Jacobian', jacobian, 'PredOpt', 'observation')
% tlen = length(tdata);
% plot(tdata, Ypred2(1:tlen) + delta2(1:tlen), 'r:', tdata, Ypred2(1:tlen) - delta2(1:tlen), 'r:')
% plot(tdata, Ypred2(tlen+1:end) + delta2(tlen+1:end), 'b:', tdata, Ypred2(tlen+1:end) - delta2(tlen+1:end), 'b:')

hold off;


end