function [Abrain, Bmax] = brainDrug_HartantoKwee_Jeffrey(A0, kb, kr, eA, eB)

    Abrain_time = 0; 
    Bmax = 0;

    C0 = [A0 ; 0];
    opt = odeset('RelTol', 1e-5, 'Events', @reachSS);
    [tt, CC] = ode45(@f, [0, 1e4], C0, opt);
    AA = CC(:, 1);
    BB = CC(:, 2);

    Abrain = Abrain_time * kb;

    figure;
    plot(tt, AA, "-", tt, BB, "-");
    xlabel('t / s');
    ylabel('Concentration / (mol / L)')
    legend('A', 'B');

    function dCdt = f(t, C)
        A = C(1, 1);
        B = C(2, 1);
        if (B > Bmax)
            Bmax = B;
        end
        if (A <= 0 && t > 0)
            if (Abrain_time == 0)
                Abrain_time = t;
            end
            dAdt = 0;
        else
            dAdt = -kb-kr*A-eA*A;
        end
        if (B <= 0 && t > 0)
            dBdt = 0;
        else
            dBdt = kr*A-eB*B;
        end
        dCdt = [dAdt; dBdt];
    end

    function [position, isterminal, direction] = reachSS(t, C)
        
        A = C(1, 1);
        B = C(2, 1);
        dCdt = f(t, C);
        dAdt = dCdt(1, 1);
        dBdt = dCdt(2, 1);
        condition = (abs(dAdt ./ A) > 0.0001 || abs(dBdt ./ B) > 0.0001);
        position = condition - 0.5;
        isterminal = 1;
        direction = -1;
    end

end