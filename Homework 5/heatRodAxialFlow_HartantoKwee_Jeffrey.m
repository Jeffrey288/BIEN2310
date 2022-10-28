function [] = heatRodAxialFlow_HartantoKwee_Jeffrey(alpha, L, Tinit, phi)

%  heatRodAxialFlow_HartantoKwee_Jeffrey(2, 1, linspace(0, 1, 30).', 1.9)

    n = size(Tinit);
    n = n(1, 1);
    l = L / n;

    opt = odeset('RelTol', 1e-5, 'Events', @stopCondition);
    [ttsol, TTsol] = ode45(@f, [0 1e5], Tinit, opt);

    figure(1);
    h = surf(linspace(0, L, n), ttsol, TTsol);
    set(h,'LineStyle','none');
    xlabel("position")
    ylabel("time")
    zlabel("temperature")
    ymax = max(max(TTsol));
    ymin = min(min(TTsol));

    for (i = 1:3:length(ttsol))
        figure(2);
        plot(linspace(0, L, n), TTsol(i, :), 'r-');
        axis([0 L ymin ymax]);
        grid on;
        xlabel('position');

        ylabel('temperature');
        F(i) = getframe();
    end


    function dTTdt = f(~, TT)
%         TT
        % E = F*Area*dt = rho(l*Area)*C*dT
        % dT/dt = F/rho/C/l = phi / l
        dTTdt = zeros(n, 1);
        dTTdt(1) = alpha .* (-TT(1) + TT(2)) ./ l^2 + phi / l;
        for i = 2 : n-1
            dTTdt(i) = alpha .* (TT(i-1) - 2*TT(i) + TT(i+1)) ./ l^2;
        end
        dTTdt(n) = alpha .* (TT(n-1) - TT(n)) ./ l^2 - phi / l;
        
    end

    function [position, terminal, direction] = stopCondition(tt, TT)
        dTTdt = f(tt, TT);
        condition = ~all(abs(dTTdt ./ TT) < 0.05);
        position = condition - 0.5;
        terminal = 1;
        direction = -1;
    end


end
