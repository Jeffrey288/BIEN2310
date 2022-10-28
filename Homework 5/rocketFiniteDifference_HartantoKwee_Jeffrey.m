function v0 = rocketFiniteDifference_Hartanto_Kwee(tground, D)

g = 9.8;
R = 6e6;
m = 10000;

n = 100;

% divide the range [0, tground] into n pieces
h = tground ./ n;
tt = linspace(0, tground, n + 1)'; % tt is [t0; t1; t2... tn]

y0 = 0;
yn = 0;

% since we know the two boundaries, we don't need to solve for y0 and yn
% initialize the yy vector, which is [y1; y2; y3... yn-1]
yyinit = ones(n - 1, 1);

% solving the algebraic equations by fsolve
yysol = fsolve(@f, yyinit);

% append the boundaries for plotting purpose
yysol = [y0; yysol; yn];

plot(tt, yysol, 'r-');
title('Finite difference method by HARTANTO KWEE, Jeffrey 20851871');

% find initial v0 by 3-point quadratic interpolation
% see "forward difference" in en.wikipedia.org/wiki/Finite_difference_coefficient
% forward difference of accuracy 2
v0 = (-1.5 .* yysol(1) + 2 .* yysol(2) - 0.5 .* yysol(3)) ./ h;

    function zz = f(yy)

        zz = ones(n-1, 1);
        
        for i = 1:n-1
        % ----------------------
        % add code below for calculating dydotdt and dydt using finite
        % differences (remember, yy is a vector of the yi's):
        if i - 1 == 0
            yiminus1 = y0;
        else
            yiminus1 = yy(i-1);
        end
        if i + 1 == n
            yiplus1 = yn;
        else
            yiplus1 = yy(i+1);
        end
        % calculate dy/dt and d^2y/dy^2 at y_i using finite difference:
        dydt = (yiplus1 - yiminus1) / 2 / h;
        dydotdt = (yiplus1 - 2 * yy(i) + yiminus1) / h^2;
        % ----------------------
        % define zz below. the system of algebraic equation we need to solve is zz = 0
        zz(i) = dydotdt + g * R^2 / (yy(i) + R)^2 + D / m * (dydt)^2 * sign(dydt);
        % ----------------------
        end

    end

end