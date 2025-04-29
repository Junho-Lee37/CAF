function [rho, pvalue] = sensitivityanalysis(lhsMatrix, interest)
%SENSITIVITYANALYSIS Sensitivity analysis
%   Return partial (ranked) correlation coefficient.
%
%     lhsMatrix: neval by nvar
%     interest: neval by ninterest
%
%     rho and pvalue: ninterest by nvar

[neval, nvar] = size(lhsMatrix);
[neval2, ninterest] = size(interest);
assert(neval == neval2, 'Dimension mismatch');

rho = zeros(ninterest, nvar);
pvalue = zeros(ninterest, nvar);
for k = 1:nvar
  [rhoTemp, pvalueTemp] = partialcorr(lhsMatrix(:, k), interest, ...
    lhsMatrix(:, [1:(k - 1) (k + 1):end]), 'type', 'Spearman');
  rho(:, k) = rhoTemp;
  pvalue(:, k) = pvalueTemp;
end
end

