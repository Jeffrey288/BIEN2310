function [pmean, psd] = polling(N, p, n, K)
% polling - Simulate polling of election involving two candidates 
% Usage: [pmean, psd] = polling(N, p, n, K), where 
% - pmean, psd are the mean and standard deviation of the estimate of p
% - N is the size of the population
% - p is the "real" fraction of population who prefers candidate A
% - n is the sample size (the number of people polled) 
% - K is the number of times the polling is repeated

numPreferingA = round(N .* p);

population = [ones(numPreferingA, 1); zeros(N - numPreferingA, 1)];

phats = zeros(K, 1); 

for (i = 1:1:K)

    sample = randsample(population, n);

    phat = sum(sample) ./ n;
    
    phats(i) = phat;
    
end
    
pmean = mean(phats);
psd = std(phats);

length(phats(phats <= 0.5)) ./ K

hold off;

% the histogram function is not very good at binning, and if a value is
% right at the boundary between two bins, it will always be placed in one
% side, resulting in a strange histogram.
% the following fix moves the values a little bit to avoid having
% many identical values that would sit right at the bin boundary.
% this is just for plotting purpose, as the curve fitting is 
% done with the original phats.

perturbed_phats = phats + rand(size(phats)) .* 2e-3 - 1e-3 
histogram(perturbed_phats, 'Normalization', 'pdf');

hold on;
pd = fitdist(phats, 'normal');

tt = linspace(min(phats), max(phats), 100);
xx = pdf(pd, tt); 
plot(tt, xx, 'r-');

xlabel("$ \hat{p} $, Fraction of respondents who prefer A", 'Interpreter', 'latex');
ylabel("$ f_{\hat{p}} $", 'Interpreter', 'latex');
ti = sprintf('Histogram of $ \\hat{p} $ \n $ E(\\hat{p}) = %.4f; \\sigma(\\hat{p}) = %.4f $ \n $ (N = %d, p = %.4f, n = %d, K = %d) $', pmean, psd, N, p, n, K);
title(ti, 'Interpreter', 'latex');

hold off;




end