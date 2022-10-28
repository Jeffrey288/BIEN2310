function [] = twoObjects(xinit, yinit)

% xinit = [x1,0 ; x2,0]
% yinit = [y1,0 ; y2,0]

m = 1;
q = 0.2;
alpha = 0.1;
delta = 0.5;

% ZZ = [x1; x2; y1; y2; xdot1; xdot2; ydot1; ydot2]

opt = odeset('RelTol', 1e-3, 'Events', @stop);
[tt, ZZ] = ode45(@f, [0, 100], [xinit; yinit; 0; 0; 0; 0], opt);

XX = ZZ(:, 1:2);
YY = ZZ(:, 3:4);

figure(1);
plot3(XX, YY, tt);
xlabel('x');
ylabel('y');
zlabel('Time');
grid on;
title('Trajectory of objects on xy plane');


xmin = min(min(XX));
xmax = max(max(XX)) + 0.01;
ymin = min(min(YY));
ymax = max(max(YY)) + 0.01;

for (i = 1:1:length(tt))
    figure(2);
    plot(XX(i, :), YY(i, :), 'o');
    axis([xmin, xmax, ymin, ymax]);
    axis equal;
    grid on;
    xlabel('x');
    ylabel('y');
    title('Force directed graph');
    F(i) = getframe();
end
movie(F, 1)

    function dZdt = f(t, Z)
        x1 = Z(1, 1);
        x2 = Z(2, 1);
        y1 = Z(3, 1);
        y2 = Z(4, 1);
        xdot1 = Z(5, 1);
        xdot2 = Z(6, 1);
        ydot1 = Z(7, 1);
        ydot2 = Z(8, 1);

        dx1dt = xdot1;
        dx2dt = xdot2;

        r = sqrt((x2 - x1) .^ 2 + (y2 - y1) .^ 2) + 1e-4;

        F1x = -q .* (x2 - x1) ./ (r .^ 3) + alpha .* (x2 - x1) - delta .* xdot1;
        F2x = -q .* (x1 - x2) ./ (r .^ 3) + alpha .* (x1 - x2) - delta .* xdot2;
 
        dxdot1dt = F1x ./ m;
        dxdot2dt = F2x ./ m;

        dy1dt = ydot1;
        dy2dt = ydot2;

        F1y = -q .* (y2 - y1) ./ (r .^ 3) + alpha .* (y2 - y1) - delta .* ydot1;
        F2y = -q .* (y1 - y2) ./ (r .^ 3) + alpha .* (y1 - y2) - delta .* ydot2;

        dydot1dt = F1y ./ m;
        dydot2dt = F2y ./ m;

        dZdt = [dx1dt; dx2dt; dy1dt; dy2dt; dxdot1dt; dxdot2dt; dydot1dt; dydot2dt];

    end

    function [position, isterminal, direction] = stop(t, Z)
        xdot1 = Z(5, 1);
        xdot2 = Z(6, 1);
        ydot1 = Z(7, 1);
        ydot2 = Z(8, 1);

        condition = (xdot1 > 0.001 || xdot2 > 0.001 || ydot1 > 0.001 || ydot2 > 0.001);

        position = condition - 0.5;
        isterminal = 1;
        direction = -1;

    end
end