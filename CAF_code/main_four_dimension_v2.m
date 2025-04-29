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

% --- Initialization ---
clc; clear all; close all;  % Clear command window, variables, and close figures

% --- Define global parameters used across ODE functions ---
global InjAT InjAL kp rF lambdaRF alpha beta delta gamma InjAP

% --- Set initial conditions for model variables ---
n0 = 0.01;  % initial tumor cell density
F0 = 1;     % initial fibroblast population
Ta0 = 0.1;  % activated T cell population
Te0 = 0;    % exhausted T cell population
Tr0 = 0;    % Treg cell population
L0  = 10;   % PD-L1 level
R0  = 10;   % PD-1 level
LR0 = 0;    % PD-1/PD-L1 complex
P0  = 0.3;  % PI3K activity
AT0 = 0; AL0 = 0; AP0 = 0;  % therapeutic agents (Treg inhibitor, PD-L1 inhibitor, PI3K inhibitor)

% --- Define CAF phenotype proportions and treatment status ---
alpha = 0.25;  % anti-immune CAF
beta  = 0.25;  % pro-immune CAF
delta = 0.25;  % pro-cancer CAF
gamma = 0.25;  % anti-cancer CAF

InjAT = 0; InjAL = 0; InjAP = 0;  % no treatment injected initially


% --- Set additional model parameters ---
lambdaRF = 0.1;   % fibroblast decay rate
kp = 1.7424;      % PI3K activation constant

tfinal = 24*30*12; % final simulation time (in hours)

n = 11; % number of divisions per parameter axis

% --- Start parameter sweep over alpha, beta, gamma values ---
for i = 1:n
    for j = 1:n
        for k = 1:n
            % Set CAF phenotype parameters
            alpha = (1/(n-1))*(i-1);
            beta  = (1/(n-1))*(j-1);
            gamma = (1/(n-1))*(k-1);
            delta = 1-alpha-beta-gamma; % ensure sum to 1

            if i+j+k < 13.5
                % --- Solve ODE system without treatment --
                InjAT = 0;
                InjAL = 0;
                InjAP = 0;
                sol = ode45(@ode1,[0 tfinal],[n0, F0, Ta0, Te0, Tr0, L0, R0, LR0, P0, AT0, AL0, AP0]);
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

                Cancer(i,j,k) = sol.y(1,ed);
        
                Tact(i,j,k) = sol.y(3,ed);
                Treg(i,j,k) = sol.y(5,ed);
                
                % --- Solve ODE system with treatment ---
                % InjAT = 0.03;
                % InjAL = 0.01;
                % InjAP = 0.0173;
                sol = ode45(@ode3,[0 tfinal],[n0, F0, Ta0, Te0, Tr0, L0, R0, LR0, P0, AT0, AL0, AP0]);
                t2  = sol.x;
                C2  = sol.y(1,:);
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
                AP2 = sol.y(12,:);

                ed2 = length(sol.y(1,:));

                Cancer2(i,j,k) = sol.y(1,ed2);
                Tact2(i,j,k) = sol.y(3,ed2);
                Treg2(i,j,k) = sol.y(5,ed2);
                index(i,j,k) = 0;
            else
                Cancer(i,j,k) = 0;
                Cancer2(i,j,k) = 0;
                Tact(i,j,k) = 0;
                Treg(i,j,k) = 0;
                Tact2(i,j,k) = 0;
                Treg2(i,j,k) = 0;
                index(i,j,k) = 1;


            end

            % --- Calculate cancer burden reduction due to treatment ---
            Reduction(i,j,k) = (Cancer(i,j,k)-Cancer2(i,j,k))/Cancer(i,j,k);
        end
    end
end


% --- Calculate normalized change in final cancer cell count ---
for i = 1:n
    for j = 1:n
        for k = 1:n
            Change(i,j,k) = (Cancer2(i,j,k)-0.025)/0.025; % normalize by baseline
        end
    end
end

% --- Volume rendering of cancer reduction ---
figure;
ax = axes;
whitebg([1,1,1]);set(gcf,'Color',[1,1,1]);
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 25)
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 25)
set(0,'defaultaxeslinewidth',2.5)
h = vol3d('CData',Reduction(:,[2:11],:)*100,'texture','3D');
view(3);  % Set to 3D view
axis tight;  % Fit axes to data
xlim([0 11]);  % X axis
ylim([0 11]);  % Y axis
zlim([0 11]);  % Z axis
daspect([1 1 1]);  % Set data aspect ratio to 1:1:1
colormap(ax, jet);  % Apply jet color map
colorbar;  % Add color bar to visualize values
% CLim set
ax.CLim = [0 100];
xlabel('\beta');
ylabel('\alpha');
zlabel('\gamma');

% --- Volume rendering of change in tumor burden ---
figure;
ax = axes;
whitebg([1,1,1]);set(gcf,'Color',[1,1,1]);
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 25)
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 25)
set(0,'defaultaxeslinewidth',2.5)
h = vol3d('CData',Change(:,[2:11],:)*100,'texture','3D');
view(3);  % Set to 3D view
axis tight;  % Fit axes to data
xlim([0 11]);  % X axis
ylim([0 11]);  % Y axis
zlim([0 11]);  % Z axis
daspect([1 1 1]);  % Set data aspect ratio to 1:1:1
colormap(ax, jet);  % Apply jet color map
colorbar;  % Add color bar to visualize values
% CLim set
ax.CLim = [-100 100];
xlabel('\beta');
ylabel('\alpha');
zlabel('\gamma');
% --- Calculate ratio of activated T cells to Tregs ---
for i = 1:n
    for j = 1:n
        for k = 1:n
            Tratio2(i,j,k) = Tact2(i,j,k)/Treg2(i,j,k)          
        end
    end
end

% --- Volume rendering of T cell ratio ---
figure;
ax = axes;
whitebg([1,1,1]);set(gcf,'Color',[1,1,1]);
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 25)
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 25)
set(0,'defaultaxeslinewidth',2.5)
h = vol3d('CData',Tratio2([2:10],[2:10],:),'texture','3D');
view(3);  % Set to 3D view
axis tight;  % Fit axes to data
daspect([1 1 1]);  % Set data aspect ratio to 1:1:1
colormap(ax, jet);  % Apply jet color map
colorbar;  % Add color bar to visualize values
alphamap('rampup');  % Apply transparency gradient
alphamap('decrease', 0.1);  % Reduce overall transparency
ax.CLim = [0 0.1];  % Set color axis limits for T cell ratio
xlabel('\beta');
ylabel('\alpha');
zlabel('\gamma');
   