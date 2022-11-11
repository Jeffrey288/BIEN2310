function NN = rabbitHistogram_HartantoKwee_jeffrey(m, p, d, K)

% rabbitHistogram_HartantoKwee_jeffrey(60, 1/30, 800, 100)

NN = zeros(K, 1);

for i=1:K
    NN(i) = rabbit_HartantoKwee_Jeffrey(m,p,d,false);
end

histogram(NN, 10, 'Normalization', 'pdf')
hold on;

% for i=1:4
for i=4:4
pd = [fitdist(NN, 'Normal');
    fitdist(NN, 'Lognormal');
    fitdist(NN, 'Gamma');
    fitdist(NN, 'Nakagami')];
colors=['b-', 'r-', 'g-', 'm-'];
nn = linspace(1, max(NN), 100);
probs = pdf(pd(i), nn);
plot(nn, probs, colors(i));
end

% legend('Histogram', 'Normal', 'Log-Normal', 'Gamma', 'Nakagami');
hey = sprintf("Rabbit Final Population \n" + ...
    "Fitted with \\textbf{Nakagami} distribution\n" + ...
     "$m=%d$, $p=%.2f$, $d=%d$\n $\\mu=%.2f$, $\\omega=%.2f$", m,p,d, pd(4).mu, pd(4).omega);
title(hey, 'Interpreter', 'latex')
hold off;

end