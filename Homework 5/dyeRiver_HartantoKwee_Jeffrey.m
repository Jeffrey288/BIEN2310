function [] = dyeRiver_HartantoKwee_Jeffrey(D, v, k, C0) 
   
    L = 10*max(sqrt(D/k), v/k)


%     Using shooting method
%     Works for > dyeRiver_HartantoKwee_Jeffrey(0.1, 0.05, 0.2, 1)
%     Doesn't work for > dyeRiver_HartantoKwee_Jeffrey(0.3, 0.1, 0.01, 0.6)

    top = 0;
    bottom = -1;
    mid = (top + bottom) / 2;
    dCdtf = simulate(mid);
    iterations = 0;
    mid = fsolve(@simulate, C0) % USE FSOLVE INSTEAD OF BISECTION SEARCH
    
%     while abs(dCdtf) > 1e-3 && iterations < 1e3
%         if (dCdtf > 0)
%             top = mid;
%         else
%             bottom = mid;
%         end
%         mid = (top + bottom) / 2;
%         dCdtf = simulate(mid);
%         iterations = iterations + 1;
%         [mid, dCdtf, iterations]
%     end
    simulate(mid, true)
   
    function res = simulate(dCdx0, showplot)
        if nargin == 1
            showplot = false;
        end

        opt = odeset('RelTol', 1e-11, 'MaxStep', 1, 'AbsTol', 1e-11);
        [XX, CC] = ode89(@g, [0 L], [dCdx0; C0], opt);

        function dCCdx = g(xx, CC)

            Cdot = CC(1, 1);
            C = CC(2, 1);
            dCdx = Cdot;
            dCdotdx = (1 / D) .* (v .* Cdot + k .* C);
            dCCdx = [dCdotdx; dCdx];
%             input("press enter to contninue");
        end
%         res = CC(end, 1); % THIS IS WRONG, USE INTERP
        res = interp1(XX, CC(:, 1), L);
        if showplot
            plot(XX, CC(:, 2), 'b-');
        end

    end

%     Using finite difference method

    n = 100;
    h = L/n;

    CC0 = C0/2 * ones(n, 1); % column vector of initial concentrations

    opt = optimoptions('fsolve', 'MaxIterations', 4e4, 'MaxFunctionEvaluations', 4e4);
    CC = fsolve(@f, CC0, opt);
    CC = [C0; CC];
    XX = linspace(0, L, n+1);

    (-1.5 .* CC(1,1) + 2 .* CC(2,1) - 0.5 .* CC(3,1)) ./ h

    hold on;
    figure(1);
    plot(XX, CC, 'r-');
    title("Hartanto Kwee Jeffrey");
    legend("Shooting Method", "Finite Difference");
    hold off;

    function zz = f(yy)
        zz = ones(n-1, 1);
        for i=1:n-1
            if i == 1
                dCdt = (yy(i+1) - C0)/(2*h);
                dCdotdt = (yy(i+1) - 2*yy(i) + C0)/(h^2);
            else 
                dCdt = (yy(i+1) - yy(i-1))/(2*h);
                dCdotdt = (yy(i+1) - 2*yy(i) + yy(i-1))/(h^2);
            end
            zz(i) = D * dCdotdt - v * dCdt - k * yy(i);
        end
        zz(n) = (yy(n) - yy(n-1))/(h);
    end

    
end