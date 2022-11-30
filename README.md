# BIEN2310 Homework


## Homework 2
Focus: basic Matlab syntax
- Foward Euler
- Backward Euler
- Cramer's Rule
- Fibonacci Series
- Stirred Tank

## Homework 3
Focus: system of first-order ODEs & population dynamics
- Lotka Volterra Model
- Bacterial Competition

## Homework 4
Focus: chemical kinetics
- Reactor Series
- Brain Drug

## Homework 5
Focus: second order ODEs & partial differential equations
- Rockets (IVP, BVP: shooting / finite difference)
- Force Directed Graph
- Dye River
- Heat Rod

## Exam 2
Improvement: maxStep could have been set to 1 or 0.1 such that the limit can be set lower for a more accurate stopping condition (but this was beyond the area I was able to edit // I could have added another odeset myself lol) \
Actually should have used speed = 0 to stop the program, and terminal flags :(

```
Some useful matlab commands

rand() % generate one random number from 0 to 1
randn() % generate one normally distributed random number

histogram(data_col, bins, 'Normalization', 'pdf')

pd = fitdist(data_col, 'exp')
pdf(pd, time) % probability density function
cdf(pd, time) % cumulative density function, equivalent to pnorm
icdf(pd, time) % quartile function, equivalent to qnorm

normpdf
normcdf
norminv

tpdf
tcdf
tinv

chi2pdf
chi2cdf
chi2inv

ttest
ttest2

modelfunc = @(params, t) rocketfunc(t, params, false); % note that t will be passed as a column
[lsbeta, residual, jacobian, CovB, mse] = nlinfit(tdata, ydata1, modelfunc, beta0);

ci = nlparci(lsbeta, residual, 'Jacobian', jacobian);

modelfunc = @(params, t) rocketfunc2(t, params, false); % note that params will be passed as a column
[lsbeta2, resnorm, residual, exitflag, output, lambda, jacobian] 
	= lsqcurvefit(modelfunc, beta0, tdata, [ydata1, ydata2]);

[pn, Sn] = polyfit(normal_study, normal_scores, 1);
yy = pn(1).*xx+pn(2)

SST_normal = mean((normal_scores - normal_mean).^2);
normal_residuals = normal_scores - (normal_study*pc(1)+pc(2));
SSE_normal = mean(normal_residuals.^2);
R2_normal = 1-SSE_normal/SST_normal;

```
