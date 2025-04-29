function dy = odefun(t, y, p)
%ODEFUN ODE description function for simplified lung cancer model
% INPUT = p (parameter values)
% OUTPUT = y (calculated values of the right hand side of the ODE model


n  = y(1,:);  % Tumor cells
% F  = y(2,:);  % Cancer associated fibroblasts
Ta = y(2,:);  % Activated T cells
Te = y(3,:);  % Exhausted T cells
Tr = y(4,:);  % Regulatory T cells (Tregs)
L  = y(5,:);  % PDL1
R  = y(6,:);  % PD1
LR = y(7,:);  % Complex
P  = y(8,:); % PI3K
% AT = y(10,:); % Treg anti-body
% AL = y(11,:); % PDL1 anti-body
% AP = y(12,:); % PI3K anti-body

% Please do not define or set parameter values here, 
% The base parameter values were set in baseparam_ode.m already and 
% changed parameters were set in prcctest.m before calling this function.



dy = [ ...
      p.r*n*(1+p.k1*P^20/(p.HP^20+P^20))*(1-n/p.K)-p.gamma1*Ta/(p.Kr+p.gamma2*Tr)*n-p.gamma*p.gamma3*p.F*n; ...
      p.rT-p.alpha1*LR^2/(p.HC^2+LR^2)*Ta+p.beta*p.beta1*p.F*Te; ...
      p.alpha1*LR^2/(p.HC^2+LR^2)*Ta-p.beta*p.beta1*p.F*Te; ...
      p.lambdaR+p.alpha*p.lambdaRF*p.F-p.muR*Tr; ... % -p.omega1*AT; ...
      -p.kp*p.alpha*p.F*L*R+p.km*LR; ... % -p.omega2*AL; ...
      -p.kp*p.alpha*p.F*L*R+p.km*LR; ... 
      p.kp*p.alpha*p.F*L*R-p.km*LR; ...
      p.lambdaP+p.delta*p.lambdaPG*p.F-p.muP*P; ... %-p.omega3*AP*P
      % p.InjAT-p.muAT*AT; ...
      % p.InjAL-p.muAL*AL; ...
      % p.InjAP-p.muAP*AP; ...
  ];
end


