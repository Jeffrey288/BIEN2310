function [] = heatRodRadialFlow_HartantoKwee_Jeffrey(alpha, L, ...
    Tinit, phi, Tenv, beta)

%  heatRodRadialFlow_HartantoKwee_Jeffrey(2, 1, 100*linspace(0, 1, 30).', 19, 30, 0.5)

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

    for (i = 1:40:length(ttsol))
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
        % For the qenv
        % qenv*Area of SURFACE*dt = rho(l*Area)*C*dT
        % qenv*2piR*l*dt = rho(l*piR^2)*C*dT
        % dT/dt = 2h*(T-Tenv)/(kR)/rhoC
        % = 2h/kR * k/rhoC * (T-Tenv)
        % = beta * alpha * (T-Tenv)
        dTTdt = zeros(n, 1);
        dTTdt(1) = alpha .* (-TT(1) + TT(2)) ./ l^2 + phi / l - ...
            alpha * beta * (TT(1) - Tenv);
        for i = 2 : n-1
            dTTdt(i) = alpha .* (TT(i-1) - 2*TT(i) + TT(i+1)) ./ l^2 - ...
                alpha * beta * (TT(i) - Tenv);
        end
        dTTdt(n) = alpha .* (TT(n-1) - TT(n)) ./ l^2 - phi / l - ...
             alpha * beta * (TT(n) - Tenv);
    end

    function [position, terminal, direction] = stopCondition(tt, TT)
        dTTdt = f(tt, TT);
        condition = ~all(abs(dTTdt ./ TT) < 0.01);
        position = condition - 0.5;
        terminal = 1;
        direction = -1;
    end


end
