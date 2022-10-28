function [Af,Bf] = bacteriaCompetition_HartantoKwee_Jeffrey(param)
    rA = param(1);
    rB = param(2);
    alpha = param(3);
    beta = param(4);
    K = param(5);
    A0 = param(6);
    B0 = param(7);

    t0 = 0;
    tf = 1e4;

    opt = odeset('RelTol', 1e-5, 'Events', @termCond);
    [tt, XX] = ode45(@model, [t0; tf], [A0; B0], opt);
%     Xf = interp1(tt, XX, tf);
%     Af = Xf(1, 1);
%     Bf = Xf(1, 2);
    Af = XX(end, 1);
    Bf = XX(end, 2);

    figure(1);
    plot(tt, XX(:, 1), "b-", tt, XX(:, 2), 'r-');
    legend("A", "B");
    title("Hartanto Kwee Jeffrey!");
    subtitle("my params are: [0.5, 0.2, 0.2, 0.6, 1000, 20, 20]");
    grid on;
    
    function dXdt = model(t, X)
        A = X(1, 1);
        B = X(2, 1);
        dAdt = rA .* (1 - (A + beta .* B) ./ K) .* A;
        dBdt = rB .* (1 - (alpha .* A + B) ./ K) .* B;
        dXdt = [dAdt; dBdt];
    end

    function [position, isterminal, direction] = termCond(t, X)
        dXdt = model(t, X);
        A = X(1, 1);
        B = X(2, 1);
        dAdt = dXdt(1, 1);
        dBdt = dXdt(2, 1);
        
        % checks relative magnitude
        condition = (abs(dAdt ./ A) < 0.00001 && abs(dBdt ./ B) < 0.00001);
        
        position = condition - 0.5;
        isterminal = 1;
        direction = 1;
    end
end