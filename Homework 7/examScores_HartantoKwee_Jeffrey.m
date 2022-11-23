function [] = examScores_HartantoKwee_Jeffrey(dataTable)
  
% examScores_HartantoKwee_Jeffrey(hw7_q4_data)

    haspaper = dataTable{:,1};
    studytimes = dataTable{:,2};
    scores = dataTable{:,3};

    cheater_scores = scores(haspaper);
    normal_scores = scores(~haspaper);

% (a), two-tailed, two-sample, equal-variance t-test

    Nc = length(cheater_scores);
    Nn = length(normal_scores);
    cheater_mean = mean(cheater_scores);
    normal_mean = mean(normal_scores);
    cheater_var = var(cheater_scores); % sample variance
    normal_var = var(normal_scores); % sample variance

    est_var = ((Nc - 1) * cheater_var + (Nn - 1) * normal_var) / (Nc + Nn - 2);
    est_n = Nc + Nn - 2;
    t_stat = (normal_mean - cheater_mean) / (sqrt(est_var * (1/Nn + 1/Nc)))
    p_value = 2 * (1 - tcdf(abs(t_stat), est_n))
    reject_hnull = p_value < 0.05;

%     🤡
%     [reject_hnull,p_value,ci,stats] = ttest2(cheater_scores, normal_scores, "Vartype","equal", 'Alpha',0.05, 'Tail','both')

% (b)

    cheater_study = studytimes(haspaper);
    normal_study = studytimes(~haspaper);
    
    [pc, Sc] = polyfit(cheater_study, cheater_scores, 1);
    [pn, Sn] = polyfit(normal_study, normal_scores, 1);
    xx = linspace(0, max(max(cheater_study), max(normal_study)) + 1, 100);
    plot(xx, pc(1).*xx+pc(2), 'b-', cheater_study, cheater_scores, 'b*');
    hold on;
    plot(xx, pn(1).*xx+pn(2), 'r-', normal_study, normal_scores, 'r*');
    hold off;
%     corrcoef(cheater_scores, cheater_study*pc(1)+pc(2)).^2
    SST_cheater = mean((cheater_scores - cheater_mean).^2);
    cheater_residuals = cheater_scores - (cheater_study*pc(1)+pc(2));
    SSE_cheater = mean(cheater_residuals.^2);
    R2_cheater = 1-SSE_cheater/SST_cheater;
    SST_normal = mean((normal_scores - normal_mean).^2);
    normal_residuals = normal_scores - (normal_study*pc(1)+pc(2));
    SSE_normal = mean(normal_residuals.^2);
    R2_normal = 1-SSE_normal/SST_normal;
    title({"Score regressed against study time", ...
        sprintf("$R^2_{cheater}=%.6f, R^2_{normal}=%.6f$", R2_cheater, R2_normal)}, ...
        'Interpreter', 'latex')
    legend("Cheaters", "", "Normal Students");
    xlabel("Study time");
    ylabel("Score")

 


end