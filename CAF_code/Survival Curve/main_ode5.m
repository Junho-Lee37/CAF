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
%% This script generates results of figure 10 in the article

global InjAT InjAL kp rF lambdaRF alpha beta delta gamma InjAP r gamma1

load('r_samples.mat');
load('Surviva_Probability_no_therapy_0.3.mat')

n0 = 0.01;
F0 = 1;
Ta0 = 0.1;
Te0 = 0;
Tr0 = 0;
L0  = 10;
R0  = 10;
LR0 = 0;
P0  = 0.3;
AT0 = 0;
AL0 = 0;
AP0 = 0;

alpha = 0;  % low0.2
beta  = 0.1;    % high0.1 
delta = 0;  % low0.8
gamma = 0.9;    % high0.9

lambdaRF  = 0.1;
kp       = 1.7424;
InjAT = 0;
InjAL = 0;
InjAP = 0;
rF = 0;

tfinal = 24*30*27;
numPatients = 500;




deathThreshold = 0.5;  % 사망 임계 종양 크기
deathTimes = zeros(numPatients, 1);  % 사망 시간 저장 배열


for i = 1:numPatients
    
    r = sort_samples(2*i);
    gamma1 = 0.000408;
    
    % alphap = rand;
    % betap  = rand;
    % deltap = rand;
    % gammap = rand;
    % sump   = alphap+betap+deltap+gammap;
    % 
    % alpha  = alphap/sump;
    % beta   = betap/sump;
    % delta  = deltap/sump;
    % gamma  = gammap/sump;
    % alpha = alphas(i);
    % beta  = betas(i);
    % delta = deltas(i);
    % gamma = gammas(i);

    % alphas(i) = alpha;
    % betas(i)  = beta;
    % deltas(i) = delta;
    % gammas(i) = gamma;
     
    InjAT = 0; %0.03;
    InjAL = 0; %0.01;
    InjAP = 0.0173;
    sol = ode45(@ode2,[0 tfinal],[n0, F0, Ta0, Te0, Tr0, L0, R0, LR0, P0, AT0, AL0, AP0]);
    t  = sol.x;
    C  = sol.y(1,:);
    F  = sol.y(2,:);
    Ta = sol.y(3,:);
    Te = sol.y(4,:);
    Tr = sol.y(5,:);
    L  = sol.y(6,:);
    R  = sol.y(7,:);
    LR = sol.y(8,:);
    P  = sol.y(9,:);
    AT = sol.y(10,:);
    AL = sol.y(11,:);
    AP = sol.y(12,:);

    ed = length(sol.y(1,:));

    Cancer(i) = sol.y(1,ed);

    Tact(i) = sol.y(3,ed);
    Treg(i) = sol.y(5,ed);
        
    t_day = t./24;

    % 사망시간
    deathIndex = find(C > deathThreshold, 1, 'first');
    if ~isempty(deathIndex)
        deathTimes(i) = t(deathIndex);
    else
        deathTimes(i) = 10e+10;  % 사망 시간이 없는 경우 NaN
    end

    

    
end




% 사망 시간에서 NaN 제거
deathTimes = deathTimes(~isnan(deathTimes));

% Kaplan-Meier 생존 곡선 그리기
[f, x, flo, fup] = ecdf(deathTimes, 'Function', 'survivor');
figure;
whitebg([1,1,1]);set(gcf,'Color',[1,1,1]);
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 25)
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 25)
set(0,'defaultaxeslinewidth',2.5)
stairs(x/(24*30), f, 'LineWidth', 2);
hold on;
stairs(x/(24*30), flo, 'LineStyle', '--');
stairs(x/(24*30), fup, 'LineStyle', '--');
xlabel('Time (Months)');
ylabel('Survival Probability');
title('Kaplan-Meier Survival Curve for Multiple Patients');
axis([0 27 0 1.05])
grid on;

save('Survival_Probability','deathTimes','Cancer','Tact','Treg','gamma1s','alphas','betas','deltas','gammas')


