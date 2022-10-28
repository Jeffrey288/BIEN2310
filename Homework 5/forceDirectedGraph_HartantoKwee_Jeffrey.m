function [] = forceDirectedGraph_HartantoKwee_Jeffrey(xinit, yinit)

assert(all(size(xinit) == size(yinit)));
n = size(xinit);
n = n(1, 1);

m = 1;
q = 0.2;
alpha = 0.1;
delta = 0.5;

% [dxdt dydt x y]
XX0 = [zeros(n,2) xinit yinit];
XX0 = XX0(:);
opt = odeset('RelTol', 1e-10, 'Events', @stop);
[tt, XX] = ode45(@f, [0 1e2], XX0, opt);

% XX(1:10, :)
POS = XX(:, 2*n+1:4*n);
XPOS = POS(:, 1:n);
YPOS = POS(:, n+1:2*n);
plot3(XPOS, YPOS, tt);
figure(1);
xlabel('x');
ylabel('y');
zlabel('Time');
grid on;
title('Trajectory of objects on xy plane');

xmin = min(min(XPOS));
xmax = max(max(XPOS)) + 0.01;
ymin = min(min(YPOS));
ymax = max(max(YPOS)) + 0.01;

for (i = 1:1:length(tt))
    figure(2);
    plot(XPOS(i, :), YPOS(i, :), 'o');
    axis([xmin, xmax, ymin, ymax]);
    axis equal;
    grid on;
    xlabel('x');
    ylabel('y');
    title('Force directed graph');
    F(i) = getframe();
end



    function dXXdt = f(t, XX)
        Xdot = [XX(1:n, 1) XX(n+1:2*n, 1)];
        X = [XX(2*n+1:3*n, 1) XX(3*n+1:4*n, 1)];

%         [Xdot X]

        dXdotdt = zeros(n, 2);
        dXdt = Xdot;
        for k = 1:n
            for i = 1:n
                if (i == k) 
                    continue
                end
            r = sqrt((X(i, 1) - X(k, 1)) .^ 2 + (X(i, 2) - X(k, 2)) .^ 2) + 1e-15;
            dXdotdt(k, :) = dXdotdt(k, :) ...
                + 1 / m * ((-q./r.^3+alpha) .* (X(i, :) - X(k, :)) ...
                - delta .* Xdot(k, :));
%             X(i, :) - X(k, :)
%             [i, X(i, :)]
%             [k, X(k, :)]
%             dXdotdt(k, :) = dXdotdt(k, :) ...
%                 + 1 / m * ((-q./r.^3+alpha) .* (X(i, :) - X(k, :)));
            end
        end
        dXXdt = [dXdotdt dXdt];
        dXXdt = dXXdt(:);

%         input("press enter to continue...");
    end

    function [position, isterminal, direction] = stop(t, XX)
        condition = ~all(abs(XX(1:2*n, :)) < 1e-3); % only test velocity
        
        position = condition - 0.5;
        isterminal = 1;
        direction = -1;
    end

end