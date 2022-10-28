function [] = gen_A_plots()

ns = 1:20;
As = zeros(1, 20);
for i = 1:20
    As(i) = reactorSeries_HartantoKwee_Jeffrey(1, 5, 0.01, 0.2, i, false);
end

plot(ns, As, 'bx');
title("Hartanto Kwee Jeffrey");
subtitle("F, V, k, A0 = 1, 5, 0.01, 0.2");
xlabel("n");
ylabel("conversion ratio");

end