function result = lhsampling(a, b, nsample)
%LHSAMPLING Latin Hypercube Sampling
%   Sample `nsample` points from [a, b).
%
%   Usage:
%     lhsampling(a, b, nsample)
%     lhsampling([a, b], nsample)

% input check
if nargin == 2
  nsample = b;
  b = a(2);
  a = a(1);
end

if a > b
  temp = a;
  a = b;
  b = temp;
end

h = (b - a) / nsample;
result = linspace(a, b, nsample + 1);
result = result(1:end-1) + h * rand(1, nsample);
result = result(randperm(nsample));
end

