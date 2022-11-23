function s2bar = verifySampleVar(mu, sigma, n, K)
    
    variances = zeros(K, 1);

    for i = 1:K
        % generate n random numbers
        nums = mu + sigma * randn(n, 1);
        variances(i) = var(nums);
    end

    s2bar = mean(variances);

    figure(1);
    mod_var = variances ./ (sigma^2) .* (n-1);
    histogram(mod_var, 'Normalization', 'pdf');
    hold on;
    pd = fitdist(mod_var, 'Gamma');
    xx = linspace(min(mod_var) - 0.1 * sigma, max(mod_var) + 0.1 * sigma, 100);
    pp = pdf(pd, xx);
    plot(xx, pp, '-m')
    title1 = sprintf('Sample variance, $ \\mu = %.2f, \\sigma = %.2f, n = %d $', mu, sigma, n);
    title2 = sprintf('fitted to Gamma distribution $ (a = %.4f, b = %.4f) $', pd.a, pd.b);
    title3 = sprintf('desired values:  $ s^2\\sim\\frac{\\sigma^2\\chi^2}{n-1} (a = %.4f, b = %.4f) $', n/2, 2);
    title({title1, title2, title3}, 'Interpreter', 'latex')
    hold off;





end