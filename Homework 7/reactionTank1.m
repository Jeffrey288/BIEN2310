function Bf = reactionTank1(tfs, FdivV, Ain, k1, showplot)

    tf = max(tfs);
    if (tf == 0)
        Bf = 0;
        return;
    end

    C0 = [0; 0];
%     opt = odeset('RelTol', 1e-6, 'Events', @term_cond);
    [tt, CC] = ode45(@f, [0, tf], C0);
    Cf = interp1(tt, CC, tfs);
    Af = Cf(:, 1);
    Bf = Cf(:, 2);
    
    if (showplot)
%         plot(tt, CC(:, 1), 'b-', tt, CC(:, 2), 'r-')
        plot(tt, CC(:, 2), 'r-')
        legend('A', 'B');
        xlabel("time");
        ylabel("conc")
    end

    function dCdt = f(t, C)
        A = C(1);
        B = C(2);
        dCdt = [FdivV*Ain-k1.*A-FdivV.*A ; ... %dAdt
                k1.*A-FdivV*B]; %dBdt
    end

    function [position, direction, isterminal] = term_cond(t, C)
        dCdt = f(t, C);
        condition = all(abs(dCdt ./ C) < 1e-3);
        position = ~condition - 0.5;
        direction = 1;
        isterminal = -1;
    end

end