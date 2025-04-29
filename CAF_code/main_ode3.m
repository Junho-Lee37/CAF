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

clc; clear all; close all;

global InjAT InjAL kp rF lambdaRF alpha beta delta gamma InjAP r

%% initial condition
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

%% parameter
alpha = 0.25;
beta  = 0.25;
delta = 0.25;
gamma = 0.25;

lambdaRF  = 0.1;
kp       = 1.7424;
InjAT = 0;
InjAL = 0;
rF = 0;
InjAP = 0;
r        = 1.2535e-4;

tfinal = 24*30*27;

%% simulation
for i = 1:2
    if (i==1)
        F0       = 0;
        
    elseif(i==2)
        F0       = 1;
       
        

        
    end

    sol = ode45(@ode2,[0 tfinal],[n0, F0, Ta0, Te0, Tr0, L0, R0, LR0, P0, AT0, AL0,AP0]);
    t  = sol.x;
    n  = sol.y(1,:);
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

    ed = length(sol.y(1,:));
    
    if (i==1)
        n1  = sol.y(1,:);
        F1  = sol.y(2,:);
        Ta1 = sol.y(3,:);
        Te1 = sol.y(4,:);
        Tr1 = sol.y(5,:);
        L1  = sol.y(6,:);
        R1  = sol.y(7,:);
        LR1 = sol.y(8,:);
        P1  = sol.y(9,:);
        AT1 = sol.y(10,:);
        AL1 = sol.y(11,:);
        t_day1 = t./24;
    elseif (i==2)
        n2  = sol.y(1,:);
        F2  = sol.y(2,:);
        Ta2 = sol.y(3,:);
        Te2 = sol.y(4,:);
        Tr2 = sol.y(5,:);
        L2  = sol.y(6,:);
        R2  = sol.y(7,:);
        LR2 = sol.y(8,:);
        P2  = sol.y(9,:);
        AT2 = sol.y(10,:);
        AL2 = sol.y(11,:);
        t_day2 = t./24;
  
    end
        

%% plot time courses of cancer cells
    figure(6)
    whitebg([1,1,1]);set(gcf,'Color',[1,1,1]);
    set(0,'DefaultAxesFontName', 'Times New Roman')
    set(0,'DefaultAxesFontSize', 25)
    set(0,'DefaultTextFontname', 'Times New Roman')
    set(0,'DefaultTextFontSize', 25)
    set(0,'defaultaxeslinewidth',2.5)
    plot(t_day/(30),n,'linewidth',2.5); hold on;
    xlabel('Months')
    ylabel('Cancer cells')
    axis([0 27 0 0.5])
    legend({'CAF-','CAF+'},'location','east')
    
end

%% plot results

factorLR = 1/LR2(end);
LR1s     = LR1(end)*factorLR;
LR2s     = LR2(end)*factorLR;

Complex_bar  = [LR1s; LR2s;];

figure(1)
whitebg([1,1,1]);set(gcf,'Color',[1,1,1]);
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 20)
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 20)
set(0,'defaultaxeslinewidth',2.5)
bar(Complex_bar)
ax = gca;
ax.XTick = [1 2 3 4];
ax.XTickLabel = {'CAF-','CAF+','TZB','ATV+TZB'};
axis([0.5 2.5 0 1.1])
ylabel('Complex')
%legend({'CSC','IL-6'},'location','east')

factorTa = 1/Ta1(end);
Ta1s     = Ta1(end)*factorTa;
Ta2s     = Ta2(end)*factorTa;

factorTe = 1/Te2(end);
Te1s     = Te1(end)*factorTe;
Te2s     = Te2(end)*factorTe;

T_bar  = [Ta1s Ta2s; Te1s Te2s];

figure(2)
whitebg([1,1,1]);set(gcf,'Color',[1,1,1]);
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 20)
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 20)
set(0,'defaultaxeslinewidth',2.5)
bar(T_bar)
ax = gca;
ax.XTick = [1 2 3 4];
ax.XTickLabel = {'Activated T cell','Exhausted T cell'};
axis([0.5 2.5 0 1.1])
ylabel('Normalized level')
legend({'CAF-','CAF+'},'location','east')

factorP = 1/P2(end);
P1s     = P1(end)*factorP;
P2s     = P2(end)*factorP;

factorTr = 1/Tr2(end);
Tr1s     = Tr1(end)*factorTr;
Tr2s     = Tr2(end)*factorTr;

TrP_bar  = [Tr1s Tr2s; P1s P2s];

figure(3)
whitebg([1,1,1]);set(gcf,'Color',[1,1,1]);
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 20)
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 20)
set(0,'defaultaxeslinewidth',2.5)
bar(TrP_bar)
ax = gca;
ax.XTick = [1 2 3 4];
ax.XTickLabel = {'Regulatory T cell','PI3K'};
axis([0.5 2.5 0 1.1])
ylabel('Normalized level')
legend({'CAF-','CAF+'},'location','east')
