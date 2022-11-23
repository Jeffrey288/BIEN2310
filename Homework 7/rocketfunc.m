function y = rocketfunc(t, params, showplot)
% ROCKETFUNC - Solve the 2nd-order ODE: d2y/dt2 = -g - D * sign(ydot) * ydot ^ 2
% Usage: ys = rocketfunc(t, params, showplot), where 
% - y is a column vector of the y positions at the time points specified by t. 
% - t is the time points (a column vector) at which the y positions are needed 
% - params is a vector [D; v0], where D is the drag coefficient, and v0 is
%   the initial velocity
% - showplot is a boolean, whether or not to show the plots
g = 9.8;
D = params(1);
v0 = params(2);
% g = params(3); % if we want to try to fit g too, uncomment
t0 = 0;
tf = max(t) .* 1.01;
[tt, YY] = ode45(@f, [t0, tf], [0; v0]);
yy = YY(:, 1); % positions
if (showplot)
    plot(tt, yy, 'r-');
    xlabel('Time / s');
    ylabel('Altitude / m');
    grid on;
    titlestr = sprintf('D=%.3g, v_0=%.3g', D, v0);
    title(titlestr);
end
y = interp1(tt, yy, t);
    function dYdt = f(t, Y)
        y = Y(1, 1);
        ydot = Y(2, 1);
        dydt = ydot;
        dydotdt = -g - D .* sign(ydot) .* (ydot .^ 2);
        dYdt = [dydt; dydotdt];
    end
end