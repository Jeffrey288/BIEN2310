function v0 = rocketShooting(tground, D)
% rocket_shooting(tground, D) - Solve the 2nd-order ODE by the shooting method: 
%       d2y/dt2 = -g  * R^2 / (y^2 + R^2) - D * sign(ydot) * ydot ^ 2
%       with boundary values y(t=0) = 0 and y(t=tground) = 0
% Usage: rocket_shooting(ground, D), where 
% - tground is the time at which the rocket returns to ground (y = 0) 
% - D is the drag coefficient 
% - v0 is the initial velocity

tolerance = 0.01;
ttest = 0;

vlow = 0;
vhigh = 1000;

numTries = 0;

while (abs(tground - ttest) > tolerance && numTries < 20)

    v0 = (vlow + vhigh) / 2 % bisection search
    
    [ttest, maxy] = rocket(v0, D, false); % don't plot while searching for solution

    if (tground > ttest)
        vlow = v0;
    end
    if (tground < ttest)
        vhigh = v0;
    end
    
    numTries = numTries + 1;
end

if (abs(tground - ttest) <= tolerance) 
  [ttest, maxy] = rocket(v0, D, true); % plot the final solution
  title('Shooting method by HARTANTO KWEE, Jeffrey 20851871');
else
  error('Cannot find solution. Probably initial velocity required is greater than 1000 m/s');
end

% alternative method using fzero
%
% [v0, deviation] = fzero(@(vtest) rocket(vtest, D, false) - tground, [0, 1000]);
% % note: if fzero finds a solution, deviation will be close to zero
% if (abs(deviation) <= tolerance) 
%  [ttest, maxy] = rocket(v0, D, true); % plot the final solution
% else
%  error('Cannot find solution. Probably initial velocity required is greater than 1000 m/s');
% end


end
