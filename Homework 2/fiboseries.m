function F = fiboseries(n)
% fiobseries. Compute the Fibonacci series 
% Usage: F = fiboseries(n), where
% - F is a column vector containing the first n Fibonacci numbers
% - n is the number of Fibonacci numbers to compute (n must be a
%   non-negative integer)
% checking to make sure n is a non-negative integer
if (n < 1 || ceil(n) ~= floor(n))
    F = [];
    return;
end
% initialize output vector to all 1's
% the first two elements are already correct
F = ones(n, 1);
% apply the recurrent relation in a loop
% note that this loop will not be executed if n < 3

% % code for part (a)
% for i = 3:1:n
%    F(i, 1) = F(i - 1, 1) + F(i - 2, 1)
% end

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

% code for part (c)
F = A^(n-2) * F
end