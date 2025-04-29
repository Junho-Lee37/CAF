function X = mylhsdesign(n, p)
%MYLHSDESIGN Generate a latin hypercube sample.
%  It is a partial implementation of LHSDESIGN in Statistics Toolbox
%  See http://www.mathworks.com/help/stats/lhsdesign.html for more detail.

X = zeros(n, p);

% Left edges of intervals
base = linspace(0, 1 - 1 / n, n);

% Permute intervals and write to each columns of X
for k = 1:p
  X(:, k) = base(randperm(n));
end

% Values should be randomly distributed within its interval
X = X + 1 / n * rand(n, p);

end

