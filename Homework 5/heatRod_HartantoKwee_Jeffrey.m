function [] = heatRod_HartantoKwee_Jeffrey(alpha, L, Tinit)

% heatRod_HartantoKwee_Jeffrey(2, 2*pi, sin(linspace(0, 2*pi, 30)).')
% heatRod_HartantoKwee_Jeffrey(2, 3*pi, sin(linspace(0, 3*pi, 30)).')
% heatRod_HartantoKwee_Jeffrey(2, 5*pi, sin(linspace(0, 5*pi, 30)).')
% heatRod_HartantoKwee_Jeffrey(2, 1, linspace(0, 1, 30).')
% heatRod_HartantoKwee_Jeffrey(2, 3*pi, 100*sin(linspace(0, 3*pi, 30)).')

    n = size(Tinit);
    n = n(1, 1);
    l = L / n;

    opt = odeset('RelTol', 1e-5, 'Events', @stopCondition);
    [ttsol, TTsol] = ode45(@f, [0 1e5], Tinit, opt);

    h = surf(linspace(0, L, n), ttsol, TTsol);
    set(h,'LineStyle','none');
    xlabel("position")
    ylabel("time")
    zlabel("temperature")
%     plot3(ttsol, l.*((1:n) - 1), TTsol);

    ymax = max(max(TTsol));
    ymin = min(min(TTsol));
    for (i = 1:1:length(ttsol))
        figure(2);
        plot(linspace(0, L, n), TTsol(i, :), 'r-');
        axis([0 L ymin ymax]);
        grid on;
        xlabel('position');
        ylabel('temperature');
        F(i) = getframe();
    end

    function dTTdt = f(~, TT)
        dTTdt = zeros(n, 1);
        dTTdt(1) = alpha .* (-TT(1) + TT(2)) ./ l^2;
        for i = 2 : n-1
            dTTdt(i) = alpha .* (TT(i-1) - 2*TT(i) + TT(i+1)) ./ l^2;
        end
        dTTdt(n) = alpha .* (TT(n-1) - TT(n)) ./ l^2;
        
    end

    function [position, terminal, direction] = stopCondition(tt, TT)
        dTTdt = f(tt, TT);
        condition = ~all(abs(dTTdt ./ TT) < 0.01);
        position = condition - 0.5;
        terminal = 1;
        direction = -1;
    end


end
