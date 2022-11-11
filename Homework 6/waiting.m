function [mu, sigma] = waiting(p, n)
% waiting - Simulate bus waiting times assuming bus arrivals follow a Poisson process 
% Usage: [mu, sigma] = waiting(p, n), where 
% - mu is the mean of the waiting times (in seconds)
% - sigma is the standard deviation of the waiting times (in seconds)
% - p is the probability of arrival in any given second
% - n is the number of trials

wt = zeros(n, 1);

for (i = 1:1:n)
    
    t = 0; % time in seconds
    while (rand(1, 1) >= p)
        t = t + 1;
    end
       
    wt(i) = t;
end

mu = mean(wt);
sigma = std(wt);

hold off;
histogram(wt, 50, 'Normalization', 'pdf');

hold on;

pd = fitdist(wt, 'exp');

tt = linspace(0, max(wt), 100);
xx = pdf(pd, tt); 
plot(tt, xx, 'r-');

xlabel('Waiting time (x) in seconds', 'Interpreter', 'latex');
ylabel('$$f_X(x)$$', 'Interpreter', 'latex');
ti = sprintf('Histogram of waiting times ($\\mu$ = %.2f; $\\sigma$ = %.2f) for n = %d, p = %.4f', mu, sigma, n, p);
title(ti, 'Interpreter', 'latex');

hold off;




end