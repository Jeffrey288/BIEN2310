function F = fiboseries_noloop_Jeffrey_HartantoKwee(n)
% fiobseries. Compute the Fibonacci series 
% Usage: F = fiboseries(n), where
% - F is a column vector containing the first n Fibonacci numbers
% - n is the number of Fibonacci numbers to compute (n must be a
%   non-negative integer)

if (n < 1 || ceil(n) ~= floor(n))
    F = [];
    return;
end

F = ones(n, 1);

low_diag = diag(ones(1, n-2));
low_rect_left = zeros(n - 2, n - 1);
low_rect_left(1:1:n-2, 1:1:n-2) = low_diag;
low_rect_right = zeros(n - 2, n - 1);
low_rect_right(1:1:n-2, 2:1:n-1) = low_diag;
low_rect = low_rect_left + low_rect_right;
A = zeros(n, n);
A(1, 1) = 1;
A(2, 2) = 1;
A(3:1:n, 1:1:n-1) = low_rect;

F = A^(n-2) * F
end