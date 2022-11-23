function s2bar = verifyVar(mu, sigma, n, K)
% verifyVar - Simulate confidence intervals of variance estimated from samples of an underlying normal distribution 
% Usage: s2bar = verifyVar(mu, sigma, n, K), where 
% - mu is the population mean
% - sigma is the population standard deviation
% - n is the sample size
% - K is the number of experiments
% - s2bar is the expected value of s2 (its mean over K experiments)

% generate sample
sample = randn(n, K) .* sigma + mu;

xbar = mean(sample);
s = std(sample);
s2 = s .^ 2; % sample variance
var = sigma .^ 2; % population variance

hold off;
chisq = (n - 1) .* s2' ./ var;
histogram(chisq, 'Normalization', 'pdf');
hold on;

pd = fitdist(chisq, 'gamma')
tt = linspace(0, max(chisq), 100);
xx = pdf(pd, tt);
plot(tt, xx, 'r-');
xlabel('$ \chi^2 $', 'Interpreter', 'latex');
ylabel('$ f_{X^2}(\chi^2) $', 'Interpreter', 'latex');
title1 = sprintf('Sample variance, $ \\mu = %.2f, \\sigma = %.2f, n = %d $', mu, sigma, n);
title2 = sprintf('fitted to Gamma distribution $ (a = %.4f, b = %.4f) $', pd.a, pd.b);
title({title1, title2}, 'Interpreter', 'latex')

s2bar = mean(s2);

end

