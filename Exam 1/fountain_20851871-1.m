function Hss = fountain_20851871(Fin, h1, h2)
% REMINDER: Rename your function to fountain_<8-digit-Student_ID>
% Do not change input arguments or return value 

% To test-run your program, try h1 and h2 in the range of 0-100, 
% and Fin in the range of 0-10.

% The following parameters are fixed. Do not change them.
A = 50; % cross sectional area of cylinder (cm^2)
rho = 1; % density of water (g / cm^3)
g = 980; % acceleration due to gravity (cm / s^2)
k = 0.0001; % proportionality constant (flow rate / hydrostatic pressure) (cm^4 s / g)
H0 = 0; % initial height (cm)


opt = odeset('RelTol', 1e-6, 'Events', @reachSS); 

% Use ode45 to solve ODE below
% ----------------------------------
% Add code here

    t0 = 0;
    tf = 1e5;

    [tout, Hout] = ode45(@f, [t0, tf], H0, opt);
 
% ----------------------------------

% Plot graph H vs t below
% Label axes
% Put your name in figure title
% ----------------------------------
% Add code here
 
figure(1);
plot(tout, Hout, "b-");
xlabel("t (s): time");
ylabel("H (cm): height of water in cylinder")
title("Hartanto Kwee Jeffrey")

% ----------------------------------

% Calculate F1 and F2 and plot two curves in one graph, 
% F1 vs t and F2 vs t, below
% Label axes and add legend
% Put your name in figure title
% ----------------------------------
% Add code here
 
 F1 = [];
 F2 = [];
 for i=1:1:length(Hout)
     H = Hout(i);
     if (H < h1)
        F1(end + 1) = 0;
     else
         F1(end + 1) = k * rho * g * (H - h1);
     end
      if (H < h2)
        F2(end + 1) = 0;
     else
         F2(end + 1) = k * rho * g * (H - h2);
     end
 end

 figure(2);
 plot(tout, F1, "b-", tout, F2, "r-");
 xlabel("t (s): time");
 ylabel("Fi (cm^3/s): volumetric flow rate)");
 legend("F1", "F2");
 title("Hartanto Kwee Jeffrey"); 
 
% ----------------------------------


% Set return value below
% ----------------------------------
% Add code here

Hss = Hout(end);

% ----------------------------------


function dHdt = f(t, H)
% Compute slope dH/dt below
% ----------------------------------
% Add code here
 
     if (H <= h1)
        dHdt = Fin;
    elseif (h1 < H && H <= h2)
        dHdt = Fin - k * rho * g * (H - h1);
    else
        dHdt = Fin - k * rho * g * (H - h1) - k * rho * g * (H - h2);
    end
    dHdt = dHdt / A;
 
% ----------------------------------
end

% Nested function to detect if steady state is reached.
% Do not change.
function [position, isterminal, direction] = reachSS(t, H)
   dHdt = f(t, H);

   condition = (abs(dHdt) > 0.0001);

   position = condition - 0.5;
   isterminal = 1;
   direction = -1;

end


end