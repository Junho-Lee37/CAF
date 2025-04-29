close all; clear all; 

p = baseparam_ode();

% Modify code here ===========================================================
% Number of simulations
neval = 10000;

% Informations for sampling variables for LHS
% One can set the range of sampling and its sampling scale here.
lhsInfo = { ...
    % 'variable_name', [lower_bound upper_bound], scaling
    % `variable_name` is a fieldname of a parameter variable `p`.
    % I Assume that initial value is not used for SA here and all sampling
    % parameters are coefficient of ODE system, that is, included in `p`.

    'r', [6.92e-5 0.0015], 'log'; ...
    'k1', [0.1 1], 'log'; ...
    % 'HP', [8 8], 'log'; ...
    % 'K', [1 1], 'log'; ...
    'gamma1', [4.08e-5 4.08e-3], 'log'; ...
    'gamma2', [5e-1 2e-0], 'log'; ...
    'gamma3', [1e-5 1e-3], 'log'; ...
    % 'Kr', [5e-1 1e-0], 'log'; ...
    % 'rF', [1e-3 2e-2], 'log'; ...
    'k2', [5 10], 'log'; ...
    % 'KF', [1 1], 'log'; ...
    % 'muF', [7.2963e-5 7.2963e-3], 'log'; ...
    'rT', [2e-6 2e-4], 'log'; ...
    % 'KT', [1 1], 'log'; ...
    'alpha1', [6e-2 0.6], 'log'; ...
    'alpha', [1e-3 1], 'log'; ...
    'beta', [1e-3 1], 'log'; ...
    % 'HC', [3 3], 'log'; ...
    'muR', [2.17e-3 2.17e-1], 'log'; ...
    'lambdaRF', [1.7e-2 1.7], 'log'; ...
    'lambdaR', [2e-3 2e-1], 'log'; ...
    'kp', [1.7424e-1 1.7424*2], 'log'; ...
    'km', [6.48e-1 6.48*2], 'log'; ...
    % 'HF', [0.4 0.4], 'log'; ...
    'lambdaP', [0.02 2], 'log'; ...
    'lambdaPG', [0.05 5], 'log'; ...
    'muP', [2.89e-3 2.89e-1], 'log'; ...
    % 'InjAT', [1e-5 0.3], 'log'; ...
    % 'InjAL', [1e-5 0.1], 'log'; ...
    % 'muAT', [2.1e-4 2.1e-2], 'log'; ...
    % 'muAL', [1.9e-4 1.9e-2], 'log'; ...
    % 'omega1', [0.1 10], 'log'; ...
    % 'oemga2', [0.01 1], 'log'; ...
    'gamma', [1e-3 1], 'log'; ...
    'delta', [1e-3 1], 'log'; ...
    'beta1', [2e-2 2e-0], 'log'; ...
    % 'm,', [2 200], 'log'; ...
    % 'm2', [0.2 20], 'log'; ...
  };

% varNames = {'n', 'F','Ta','Te', 'Tr', 'L', 'R', 'KR', 'P', 'AT', 'AL'};
varNames = {'n','Ta','Te', 'Tr', 'L', 'R', 'LR', 'P'};
% Set all time points of interest
% Of course, zero (or any starting time) should be included.
tspan = [0 24 24*10 24*20*1];

% Initial condition
n0 = 0.01; F0 = 0.1; Ta0 = 0.1; Te0 = 0; Tr0 = 0; L0  = 10; R0  = 10; LR0 = 0; P0  = 0.3; AT0 = 0; AL0 = 0;
% y0 = [n0; F0; Ta0; Te0; Tr0; L0; R0; LR0; P0; AT0; AL0;];
y0 = [n0; Ta0; Te0; Tr0; L0; R0; LR0; P0;];
% ============================================================================

% Variables for convinience
nvars = numel(y0);
interest = zeros(neval, nvars * (numel(tspan) - 1));
interestNames = cell(numel(varNames), (numel(tspan) - 1));
for k = 1:numel(varNames)
  for l = 2:numel(tspan)
    interestNames{k, l - 1} = sprintf('%s(%g)', varNames{k}, tspan(l));
  end
end

% Make LHS Matrix
%
% `lhsdesign` is a function included in Statistics Toolbox
% See http://www.mathworks.com/help/stats/lhsdesign.html .
%
% I wrote `lhssampling` when I did not know about `lhsdesign`,
% actually it does the similar tasks (with different syntax).
% So `lhssampling` functions was not used in this script any more.
%
% I write and send a function `mylhsdesign`, a partial implementation
% of `lhsdesign`.

nsamplingVars = size(lhsInfo, 1);
% lhsMatrix = lhsdesign(neval, nsamplingVar);
lhsMatrix = mylhsdesign(neval, nsamplingVars);
for k = 1:nsamplingVars
  % (a, b) is a range of `k`-th variable
  a = lhsInfo{k, 2}(1); % lower bound
  b = lhsInfo{k, 2}(2); % upper bound

  % Because `(my)lhsdesign` returns points lie in (0, 1),
  % we should scale each variable to lie in its own range.
  switch (lhsInfo{k, 3})
    case 'linear'
      lhsMatrix(:, k) = a + (b - a) * lhsMatrix(:, k);
    case 'log'
      % The samples distributed with log scale in the range of (a, b).
      % For example, if the range is (0.01, 100) and n = 4,
      % then
      %   one point is in (0.01, 0.1), that is, (10^-2, 10^-1)
      %   second point is in (0.1, 1),
      %   third point is in (1, 10),
      %   the last is in (10, 100),
      %   (of course, with random permutations)
      % rather than (0.01, 25.0075), (25.0075, 50.005), (50.005, 75.0025),
      % and (75.0025, 100)
      lhsMatrix(:, k) = exp(log(a) + (log(b) - log(a)) * lhsMatrix(:, k));
    otherwise
      error('Invalid scale for LHS.');
  end
end

%% Simulations

tic;
% Copy parameters to another variable (to save default values of parameters).
% If one do not use default value, ignore this.
pcopy = p; % Note here that the basic parameters were saved in the array 'p'
for k = 1:neval
  % Set parameters for the simulation of a current loop
  % If `str == 'paramname'`, then
  %     pcopy.(str)
  % is equivalent to
  %     pcopy.paramname
  % . Note that the parentheses after the comma are essential here.
  for l = 1:nsamplingVars
    pcopy.(lhsInfo{l, 1}) = lhsMatrix(k, l);
  end

  % Feb, 04, 2016, YJ comment: 
  % It is important to note that odefun takes 'pcopy' array after changing 
  % parameter values above. Therefore, one should have a spot in the input 
  % arrays for 'pcopy' and should not define parameter values in odefun 
  % subroutine
  
  [t, Y] = ode45(@(t_, y_) odefun(t_, y_, pcopy), tspan, y0);
  interest(k, :) = reshape(Y(2:end, :)', [], 1);
end
toc;

filename = sprintf('prcc-%s-neval-%d', ...
  datestr(now(), 'yyyy-mm-dd-HHMMSS'), neval);
save(filename);

%% PRCC calculation
% We need two matrices:
%   `lhsMatrix` is of size `neval` by `nsamplingVars` and
%   `interest` is of size `neval` by `ninterest`
%   where `ninterest` is the number of total interests.
%   For example, if we have four variables and we want to see at the time
%   t = 10, 20, and 100, we have 12 interests.
[rho, pvalue] = sensitivityanalysis(lhsMatrix, interest);


%graqphic 
fid = fopen([filename '.csv'], 'w');
fprintf(fid, '(rho)');
for k = 1:nsamplingVars
  fprintf(fid, ',%s', lhsInfo{k, 1});
end
fprintf(fid, '\n');
for l = 1:numel(interestNames)
  fprintf(fid, '%s', interestNames{l});
  for k = 1:nsamplingVars
    fprintf(fid, ',%f', rho(l, k));
  end
  fprintf(fid', '\n');
end

fprintf(fid, '(p-value)');
for k = 1:nsamplingVars
  fprintf(fid, ',%s', lhsInfo{k, 1});
end
fprintf(fid, '\n');
for l = 1:numel(interestNames)
  fprintf(fid, '%s', interestNames{l});
  for k = 1:nsamplingVars
    fprintf(fid, ',%g', pvalue(l, k));
  end
  fprintf(fid', '\n');
end
fclose(fid);

hf = figure();
plotcc(rho, pvalue, lhsInfo(:, 1), interestNames);
