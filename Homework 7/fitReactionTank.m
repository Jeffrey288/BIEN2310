function [] = fitReactionTank(data)
    
% fitReactionTank(hw7_q5_data)

    TT = data{:, 1};
    B1 = data{:, 2};
    meanB1 = mean(B1);
    SSTB1 = mean((B1 - meanB1).^2);
    B2 = data{:, 3};
    meanB2 = mean(B2);
    SSTB2 = mean((B2 - meanB2).^2);

    modelfunc1 = @(k1, t) reactionTank1(t, 1, 1, k1, false);
    modelfunc2 = @(k2, t) reactionTank2(t, 1, 1, k2, false);

    figure(1);
    plot(TT, B1, "bx");
    hold on;
    [k1_fitted, residuals1, jacobian1, CovB1, mse1] = ...
        nlinfit(TT, B1, modelfunc1, 0.5);
    reactionTank1(max(TT), 1, 1, k1_fitted, true);
    SSE1 = mean(residuals1 .^ 2);
    Rsq1 = 1 - SSE1 / SSTB1;

    [k2_fitted, residuals2, jacobian2, CovB2, mse2] = ...
        nlinfit(TT, B1, modelfunc2, 0.5);
    reactionTank2(max(TT), 1, 1, k2_fitted, true);
    SSE2 = mean(residuals2 .^ 2);
    Rsq2 = 1 - SSE2 / SSTB1;
    
    legend("data", "model 1", "model 2");
    title({"B1 fitted onto two reaction models", ...
        sprintf("$R^2_1 = %.6f, R^2_2 = %.6f$", Rsq1, Rsq2)}, ...
        "Interpreter","latex");
    colororder({'red', 'magenta', 'green'});
    hold off;

    % --------

    % https://stackoverflow.com/questions/13540418/convert-cell-array-of-cell-arrays-to-matrix-of-matrices
    modelfunc1 = @(k1, t) [reactionTank1(t, 1, 1, k1(1), false), reactionTank1(t, 1, 1.2, k1(2), false)];
    modelfunc2 = @(k2, t) [reactionTank2(t, 1, 1, k2(1), false), reactionTank2(t, 1, 1.2, k2(2), false)];
    [k1_fitted, resnorm, residuals1, exitflag, output, lambda, jacobian] = ...
        lsqcurvefit(modelfunc1, [0, 0], TT, [B1, B2]);
    [k2_fitted, resnorm, residuals2, exitflag, output, lambda, jacobian] = ...
        lsqcurvefit(modelfunc2, [0, 0], TT, [B1, B2]);

    figure(2);
    a1 = subplot(1,2,1);
    plot(TT, B1, 'bx');
    hold on;
    plot(TT, B2, "rx");
    reactionTank1(max(TT), 1, 1, k1_fitted(1, 1), true);
    reactionTank1(max(TT), 1, 1.2, k1_fitted(1, 2), true);
    hold off;
    legend('$k_1=1$', '$k_1=1.2$', 'Interpreter', 'latex');
    subtitle({"lsqcurvefit of Model 1", ...
        sprintf("$R^2_{1(k=1)}=%.6f, R^2_{1(k=1.2)}=%.6f$", ...
        1 - mean(residuals1(:,1) .^ 2)/SSTB1, ...
        1 - mean(residuals1(:,2) .^ 2)/SSTB2)}, ...
        'Interpreter', "latex");

    a1 = subplot(1,2,2);
    plot(TT, B1, 'bx');
    hold on;
    plot(TT, B2, "rx");
    reactionTank2(max(TT), 1, 1, k2_fitted(1, 1), true);
    reactionTank2(max(TT), 1, 1.2, k2_fitted(1, 2), true);
    hold off;
    legend('$k_1=1$', '$k_1=1.2$', 'Interpreter', 'latex');
    subtitle({"lsqcurvefit of Model 2", ...
        sprintf("$R^2_{2(k=1)}=%.6f, R^2_{2(k=1.2)}=%.6f$", ...
        1 - mean(residuals2(:,1) .^ 2)/SSTB1, ...
        1 - mean(residuals2(:,2) .^ 2)/SSTB2)}, ...
        'Interpreter', "latex");
    
    
end