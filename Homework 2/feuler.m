
function [tt, yy] = feuler(fun, y0, tf, h)
% feuler. ODE solver using the forward Euler method
% Usage: [tt yy] = feuler(fun, y0, tf, h), where
% - [tt, yy] contains data points of the solved ODE y(t).
% - fun is the function handle of dy/dt = f(t, y)
% - y0 is the initial value y(t = 0)
% - h is the step size
% - tf is the time at which to stop the simulation 
% Put the first point (0, y0) in the vectors holding the data points
tt = [0];
yy = [y0];
while (true)
    % extract the most recent data point
    ti = tt(end);
    yi = yy(end);
   
    % find next data point
    tiplus1 = ti + h;   
    slopei = fun(ti, yi);
    yiplus1 = yi + h .* slopei;    
   
    % append the new data point to the vectors 
    tt(end + 1) = tiplus1;
    yy(end + 1) = yiplus1;
    % check if we have simulated past tf already
    if (tiplus1 > tf)
        break;
    end
    
end
% uncomment below to show plot
% plot(tt, yy, 'r-+');
% xlabel('t');
% ylabel('y');
end