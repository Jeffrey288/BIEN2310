function ys = rocketfunc2(t, params, showplot)
% ROCKETFUNC2 - Solve the 2nd-order ODE: d2y/dt2 = -g - D * sign(ydot) * ydot ^ 2
% for two different initial velocities
% Usage: ys = rocketfunc2(t, params, showplot), where 
% - ys are the y positions at the time points specified by t. It is a
%   matrix with two columns (y_1 and y_2), y_1 and y_2 are the y positions
%   reached using the initial velocities v0_1, v0_2 respectively, 
%   and the rows corresponding to each time point
% - t is the time points (a column vector) at which the y positions are needed 
% - params is a vector [D; v0_1; v0_2], where D is the drag coefficient,
%   v0_1 and v0_2 are the two initial velocities to simulate
% - showplot is a boolean, whether or not to show the plots

g = 9.8;

D = params(1);
v0_1 = params(2);
v0_2 = params(3);

t0 = 0;
tf = max(t) .* 1.01;

[tt1, YY1] = ode45(@f, [t0, tf], [0; v0_1]);
[tt2, YY2] = ode45(@f, [t0, tf], [0; v0_2]);

yy1 = YY1(:, 1); % positions
yy2 = YY2(:, 1);

if (showplot)
    plot(tt1, yy1, 'r-', tt2, yy2, 'b-');
    xlabel('Time / s');
    ylabel('Altitude / m');
    grid on;
    titlestr = sprintf('D=%.3g, v_0,1=%.3g, v_0,2=%.3g', D, v0_1, v0_2);
    title(titlestr);
end

y1 = interp1(tt1, yy1, t);
y2 = interp1(tt2, yy2, t);

ys = [y1, y2];

    function dYdt = f(t, Y)

        y = Y(1, 1);
        ydot = Y(2, 1);

        dydt = ydot;
        dydotdt = -g - D .* sign(ydot) .* (ydot .^ 2);

        dYdt = [dydt; dydotdt];
    end

end