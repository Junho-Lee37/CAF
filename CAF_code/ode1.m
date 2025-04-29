%{ 
Developed by:
Junho Lee
 
Copyright (c) 2025 Junho Lee
 
All rights reserved.
 
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are 
met: 
 
1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer. 
2. Redistributions in binary form must reproduce the above copyright 
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution. 
 
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS 
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR 
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%}

function dy=ode1(t,y)

global InjAT InjAL kp rF lambdaRF alpha beta delta gamma InjAP

r        = 1.2535e-4; %0.015  r = ln(2)/T_d 0.0015 - 6.92e-5
k1       = 1;
HP       = 8;
K        = 1;
gamma1   = 0.000408; % 0.07
gamma2   = 2;    % 3
gamma3   = 0.0001; % 0.03 
Kr       = 0.5;
k2       = 10;
KF       = 1;
muF      = 7.2963e-4;  % 평균수명 57일 CHICK EMBRYO FIBROBLASTS SENESCENCE IN VITRO: PATTERN OF CELL DIVISION AND LIFE SPAN AS A FUNCTION OF CELL DENSITY
lambdaG  = 5;
muG      = 4.621;       % 9 min 
rT       = 0.00002;
KT       = 1;
alpha1    = 0.6;
HC       = 3;
muR      = 0.0217;    % 32~120h Pharmacokinetics of Therapeutic Tregs
lambdaR  = 0.02;
km       = 6.48;      % Design, synthesis and biological evaluation of PD-1 derived peptides as inhibitors of PD-1/PD-L1 complex formation for cancer therapy
HF       = 0.4;
lambdaP  = 0.2;
lambdaPG = 0.5;
muP      = 0.0289;    % 24h Akt forms an intracellular complex with heat shock protein 90 (Hsp90) and Cdc37 and is destabilized by inhibitors of Hsp90 function
muAT     = 0.0021;    % 2-3 weeks Comprehensive analysis of current approaches to inhibit regulatory T cells in cancer
muAL     = 0.0019;    % 15 day Safety and Activity of Anti–PD-L1 Antibody in Patients with Advanced Cancer
omega1   = 1;
omega2   = 0.1;
beta1    = 0.2;
muAP     = 0.0173;    % 40 h Phase I Dose Escalation Study of Taselisib (GDC-0032), an Oral PI3K Inhibitor, in Patients with Advanced Solid Tumors
omega3   = 1;

m  = 20;
m2 = 2;

n  = y(1);  % Tumor cells
F  = y(2);  % Cancer associated fibroblasts
Ta = y(3);  % Activated T cells
Te = y(4);  % Exhausted T cells
Tr = y(5);  % Regulatory T cells (Tregs)
L  = y(6);  % PDL1
R  = y(7);  % PD1
LR = y(8);  % Complex
P  = y(9);  % PI3K
AT = y(10); % Treg anti-body
AL = y(11); % PDL1 anti-body
AP = y(12); % PI3K anti-body


dy = [r*n*(1+k1*P^m/(HP^m+P^m))*(1-n/K)-gamma1*Ta/(Kr+gamma2*Tr)*n-gamma*gamma3*F*n;
      0;
      rT-alpha1*LR^m2/(HC^m2+LR^m2)*Ta+beta*beta1*F*Te;
      alpha1*LR^m2/(HC^m2+LR^m2)*Ta-beta*beta1*F*Te;
      lambdaR+alpha*lambdaRF*F-muR*Tr-omega1*AT*Tr;
      -kp*alpha*F*L*R+km*LR-omega2*AL*L;
      -kp*alpha*F*L*R+km*LR;
      kp*alpha*F*L*R-km*LR;
      lambdaP+delta*lambdaPG*F-muP*P-omega3*AP*P;
      InjAT-muAT*AT;
      InjAL-muAL*AL;
      InjAP-muAP*AP;
      ];
    
    

