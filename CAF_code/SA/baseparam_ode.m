function p = baseparam_ode(modelType)
%BASEPARAM

% OUTPUT = p (parameters)

p.nvars = 8;

if nargin < 1
  modelType = 'ode';
end

p.F        = 1;
p.r        = 1.2535e-4;
p.k1       = 1;
p.HP       = 8;
p.K        = 1;
p.gamma1   = 0.000408; % 0.05
p.gamma2   = 2;    % 3
p.gamma3   = 0.0001;
p.Kr       = 0.5;    % 0.5
p.rF       = 0.01;
p.k2       = 10;
p.KF       = 1;
p.muF      = 7.2963e-4;  % CHICK EMBRYO FIBROBLASTS SENESCENCE IN VITRO: PATTERN OF CELL DIVISION AND LIFE SPAN AS A FUNCTION OF CELL DENSITY
p.lambdaG  = 5;
p.muG      = 4.621;       % 9 min 
p.rT       = 0.00002;
p.KT       = 1;
p.alpha1   = 0.6;
p.alpha    = 0.25;
p.beta     = 0.25;
p.HC       = 3;
p.muR      = 0.0217;    % 32~120h Pharmacokinetics of Therapeutic Tregs
p.lambdaRF  = 0.17;
p.lambdaR  = 0.02;
p.kp       = 1.7424;    % Design, synthesis and biological evaluation of PD-1 derived peptides as inhibitors of PD-1/PD-L1 complex formation for cancer therapy
p.km       = 6.48;      % Design, synthesis and biological evaluation of PD-1 derived peptides as inhibitors of PD-1/PD-L1 complex formation for cancer therapy
p.HF       = 0.4;
p.lambdaP  = 0.2;
p.lambdaPG = 0.5;
p.muP      = 0.0289;    % 24h Akt forms an intracellular complex with heat shock protein 90 (Hsp90) and Cdc37 and is destabilized by inhibitors of Hsp90 function
p.InjAT    = 0.1;
p.InjAL    = 0.1;
p.muAT     = 0.0021;
p.muAL     = 0.0019;    % 15 day
p.omega1   = 1;
p.omega2   = 0.1;
p.gamma    = 0.25;
p.delta    = 0.25;
p.beta1    = 0.2;
p.muAP     = 0.01;
p.omega3   = 1;

p.m  = 20;
p.m2 = 2;




end

