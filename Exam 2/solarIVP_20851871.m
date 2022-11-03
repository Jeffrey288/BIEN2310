function tau = solarIVP_20851871(R0, v0, showplot)
% REMINDER: Rename your function to solarIVP_<8-digit-Student-ID>
% Do not change input arguments or return value 

GM = 2.96e-4; % (AU^3 / d^2)

opt = odeset('RelTol', 1e-10, 'Events', @completeCycle); 

% Use ode45 to solve ODE. 
% Set return value tau
% ----------------------------------
leftOrigin = 0; % a flag used by event

% Add code here
% solarIVP_20851871(1, 0.0175, true)
% solarIVP_20851871(1, 0.0242, true)
% at the aphelion,
% x = R0, y = 0
% v_x = 0, v_y = v0
% store os [dxdt dydt x y]
ZZ0 = [0; v0; R0; 0];
[tt, ZZ] = ode45(@f, [0, 1e7], ZZ0, opt);
XX = ZZ(:, 3);
YY = ZZ(:, 4);
% the final tt
tau = max(tt);
 
% ----------------------------------

if (showplot)

    % Plot orbit in xy-plane in blue, sun as red circle
    % Remember to label axes

    figure(1);
    % ----------------------------------
    % Add code here
    plot(XX, YY, 'b-', 0, 0, 'ro');
    

    % ----------------------------------
    axis equal; % don't change. make it easier to see orbit shape
    grid on; 
    
    % For Part (c), create movie of planet orbiting around the sun.
    figure(2);
    % ----------------------------------
    % Add code here
    
        xmin = min(XX);
        xmax = max(XX);
        ymin = min(YY);
        ymax = max(YY);
        
        % use interp1
        for (t = linspace(1, tau, 5e2))
            figure(2);
            x = interp1(tt, XX, t);
            y = interp1(tt, YY, t);
            plot(x, y, 'bo',  0, 0, 'ro');
            axis([xmin, xmax, ymin, ymax]);
            grid on;
            xlabel('x');
            ylabel('y');
            title('Hey');
        end
 


    % ----------------------------------  
end

% Nested function to compute derivatives
function dZdt = f(t, Z)
    % ----------------------------------
    % Add code here
    
    xdot = Z(1,1);
    ydot = Z(2,1);
    x = Z(3,1);
    y = Z(4,1);
    
    dxdt = xdot;
    dydt = ydot;
    r = sqrt(x.^2+y.^2);
    dxdotdt = -GM*x./r.^3;
    dydotdt = -GM*y./r.^3;

    dZdt = [dxdotdt; dydotdt; dxdt; dydt];    

    % ----------------------------------
end

% Nested function to stop ode45 if planet completes cycle
function [position, isterminal, direction] = completeCycle(t, Z)
    % ----------------------------------
    % Add code here
    xdot = Z(1,1);
    ydot = Z(2,1);
    x = Z(3,1);
    y = Z(4,1);
        
    % stop condition
    limit = 0.05;
%     if (abs(x-R0)./R0 < 0.05 && abs(y)./R0 < 0.05)
%         hey = [abs(x-R0)./R0, abs(y-0)./R0, leftOrigin, (abs(x-R0)./R0 < limit) && (abs(y-0)./R0 < limit) && leftOrigin]
%     end
    condition = (abs(x-R0)./R0 < limit) && (abs(y-0)./R0 < limit) && leftOrigin;
    if ((abs(x-R0) > limit) || (abs(y-0) > limit)) 
        leftOrigin = 1;
    end
    position = ~condition - 0.5;
    isterminal = 1;
    direction = 0;
 
    % ----------------------------------
end
    

end