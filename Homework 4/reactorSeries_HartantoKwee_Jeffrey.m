function conv = reactorSeries_HartantoKwee_Jeffrey(F, V, k, A0, n, showplot)
    
    X0 = zeros(1, 2*n);

    opt = odeset('RelTol', 1e-5, 'Events', @reachSS);
    [tout, Xout] = ode45(@f, [0, 1e4], X0, opt);
   
    if showplot
        figure;
        hold on;
        legends = num2str((1:n).', 'B_{%d}');
        xlabel("time");
        ylabel("B_i");
        title("hartanto kwee jeffrey");
        subtitle("params used: F="+F+" V="+V+" k="+k+" A0="+A0+" n="+n)
        for i = 1:n
            plot(tout, Xout(:,n+i), '-');
        end   
        legend(legends);
        hold off;
    end
    
    conv = Xout(end,2*n)/A0;

    function dXdt = f(t, X)
        A = X(1:n, 1).';
        B = X(n+1:2*n, 1).';
        dAdt = zeros(1, n);
        dBdt = zeros(1, n);
        dAdt(1,1) = n*A0*F/V-(k+n*F/V)*A(1);
        dBdt(1,1) = k*A(1)-n*F/V*B(1);
        for i = 2:1:n
            dAdt(1,i)=n*A(1,i-1)*F/V-(k+n*F/V)*A(1,i);
            dBdt(1,i)=n*B(1,i-1)*F/V+k*A(i)-n*F/V*B(1,i);
        end
        dXdt = [dAdt dBdt].';
    end

    function [position, isterminal, direction] = reachSS(t, X)
        dXdt = f(t, X);
        condition = max(abs(dXdt ./ X)) > 1e-4;
        position = condition - 0.5;
        isterminal = 1;
        direction = -1;
    end

end