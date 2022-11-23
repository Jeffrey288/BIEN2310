function FPR = verifyFTest_HartantoKwee_Jeffrey(mu, sigma, n1, n2, K, alpha)

% verifyFTest_HartantoKwee_Jeffrey(9, 2, 10, 20, 100, 0.04)

    f_stats = zeros(K, 1);
    f_probs = zeros(1, K);
    successes = zeros(K, 1); 
    % times the null hypothesis is not wrongly rejected

    for i = 1:K
        nums1 = mu + sigma .* randn(n1, 1);
        nums2 = mu + sigma .* randn(n2, 1);
        var1 = var(nums1);
        var2 = var(nums2);
        f_stats(i) = var1/var2;

        f_probs(1, i) = fcdf(f_stats(i), n1-1, n2-1);
        successes(i) = (alpha/2 < f_probs(1, i)) && (f_probs(1, i) < 1-alpha/2);
    end
%     f_probs

%     fcdf(alpha/2, n1-1, n2-1)
%     fcdf(1-alpha/2, n1-1, n2-1)

    histogram(f_stats, 'Normalization', 'pdf');
    xx = linspace(0, max(f_stats) + 0.1, 100);
    pp = fpdf(xx, n1-1, n2-1);
    hold on;
    plot(xx, pp, 'm-');
    plot(ones(2, 1) * finv(alpha/2, n1-1, n2-1), [0; 1], '--r');
    plot(ones(2, 1) * finv(1-alpha/2, n1-1, n2-1), [0; 1], '--r');
    hold off;
    title({"F-test", ...
        sprintf("(mu, sigma, n1, n2, K, alpha)=(%.2f,%.2f,%d,%d,%d,%.4f)",mu, sigma, n1, n2, K, alpha)},...
        'Interpreter', 'latex');
    xlabel("F-statistic")
    ylabel("probability")
%     successes
    FPR = 1 - sum(successes) / K;
end