function [] = van_der_Pol_oscillator(mu, y0, tf, F)

    opt = odeset('RelTol', 1e-4);
    [tt, YY] = ode45(@f, [0 tf], [0; y0], opt);
    [peaks, locs] = findpeaks(YY(:, 2));
    diffs = diff(tt(locs));
    period = mean(diffs)
    figure(1);
    plot(tt, YY(:, 2));
    xlabel("time");
    ylabel("y");
    figure(2);
    plot(YY(:, 2), YY(:, 1));
    xlabel("y");
    ylabel("dy/dt");
    

    function dYYdt = f(tt, YY)
        Y = YY(2, 1);
        Ydot = YY(1, 1);
        dYdt = Ydot;
        dYdotdt = mu .* (1 - Y.^2) * dYdt - Y + F(tt);
        dYYdt = [dYdotdt; dYdt];
    end

end