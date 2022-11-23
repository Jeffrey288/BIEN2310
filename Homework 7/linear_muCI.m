function ci = linear_muCI(x0, p, xdata, ydata)
% linear_muCI - Calculate the 95% confidence interval of regressed y0 at x0 in linear regression 
% Usage: ci = linear_muCI(x0, p, xdata, ydata), where 
% - x0 is the x-coordinate of the regressed point where the confidence interval is desired
% - p is the fitted polynomial from polyfit in the form [slope, intercept]
% - xdata is the x-coordinates of the data points
% - ydata is the y-coordinates of the data points
% - ci is the confidence interval in the form [ciLow, ciHigh]
%
% Also works if x0 is a column vector, in which case ci is a matrix, with 
% rows of [ciLow, ciHigh] of the corresponding regressed points.


n = length(xdata);

xbar = mean(xdata);

yhat = polyval(p, xdata);

residuals = ydata - yhat;

sigma2 = sum(residuals .^ 2) ./ (n - 2);

sem = sqrt(sigma2 .* (1 ./ n + (x0 - xbar) .^ 2 ./ (sum(xdata .^ 2) - n .* xbar .^ 2)));

tU = myTinv(0.975, n - 2);

tL = myTinv(0.025, n - 2);

y0 = polyval(p, x0);

ciU = y0 + tU .* sem;

ciL = y0 + tL .* sem;

ci = [ciU, ciL];

end