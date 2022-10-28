function [tground, maxy] = rocket(v0, D, showplot)
% rocket(v0, D, showplot) - Solve the 2nd-order ODE: d2y/dt2 = -g  * R^2 / (y^2 + R^2) - D * sign(ydot) * ydot ^ 2
% Usage: rocket(v0, D), where 
% - v0 is the initial velocity 
% - D is the drag coefficient 
% - showplot is a boolean, whether or not to show the plots

R = 6e6; 
g = 9.8;

t0 = 0;
t1 = 100;

opt = odeset('MaxStep', 1);

[tt, YY] = ode45(@f, [t0, t1], [0.6; -0.0077], opt);

yy = YY(:, 1); % positions
vv = YY(:, 2); % velocities = dy/dt

if (showplot)
    figure(1);
    plot(tt, yy, 'r-');
    xlabel('Time / s');
    ylabel('Altitude / m');
    grid on;
    titlestr = sprintf('D=%.3g, v_0=%.3g', D, v0);
    title(titlestr);
 %   set(gca, 'FontSize', 14);
    
%    figure(2);
%    plot(tt, vv, 'g-');
%    xlabel('Time / s');
%    ylabel('Velocity / (m/s)');
%    grid on;
 %   set(gca, 'FontSize', 14);

end

tground = tt(end);
maxy = max(yy);

    function dYdt = f(t, Y)

        y = Y(1, 1);
        ydot = Y(2, 1);

        dydt = ydot;
        dydotdt = (1/0.3)*(0.1*ydot+0.01*y);

        dYdt = [dydt; dydotdt];
    end

    function [position, isterminal, direction] = reachGround(t, Y)
       y = Y(1, 1); 
       % stop the simulation when y reaches 0 again (from positive)
       position =  y;
       isterminal = 1; 
       direction = -1; 
    end

end