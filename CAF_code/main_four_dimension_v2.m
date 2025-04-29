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

global InjAT InjAL kp rF lambdaRF alpha beta delta gamma InjAP

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

alpha = 0.25;
beta  = 0.25;
delta = 0.25;
gamma = 0.25;

InjAT = 0;
InjAL = 0;
InjAP = 0;


rF       = 0.01;
lambdaRF  = 0.1;
kp       = 1.7424;

tfinal = 24*30*12;

n = 11;

for i = 1:n
    for j = 1:n
        for k = 1:n
            alpha = (1/(n-1))*(i-1);
            beta  = (1/(n-1))*(j-1);
            gamma = (1/(n-1))*(k-1);
            delta = 1-alpha-beta-gamma;

            if delta >= 0
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

            Reduction(i,j,k) = (Cancer(i,j,k)-Cancer2(i,j,k))/Cancer(i,j,k);
        end
    end
end

for i = 1:n
    for j = 2:n
        for k = 1:n
            if index(i,j,k) == 1
                Change(i,j,k) = (Cancer2(i,j,k)-0.02)/0;
            else
            Change(i,j,k) = (Cancer2(i,j,k)-0.02)/0.02;
            end
        end
    end
end

for i = 1:n
    for j = 1:n
        for k = 1:n
            Change(i,j,k) = (Cancer2(i,j,k)-0.025)/0.025;
        end
    end
end

% 볼륨 렌더링
figure;
ax = axes;
whitebg([1,1,1]);set(gcf,'Color',[1,1,1]);
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 25)
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 25)
set(0,'defaultaxeslinewidth',2.5)
h = vol3d('CData',Reduction(:,[2:11],:)*100,'texture','3D');
view(3);  % 3D 뷰 설정
axis tight;  % 축 조정
xlim([0 10]);  % X축 범위
ylim([0 9]);  % Y축 범위
zlim([0 10]);  % Z축 범위
daspect([1 1 1]);  % 데이터 종횡비 설정
colormap(ax, jet);  % 색상 맵 설정
colorbar;  % 색상 막대 추가
% CLim 속성 설정
ax.CLim = [0 100];
xlabel('\beta');
ylabel('\alpha');
zlabel('\gamma');

figure;
ax = axes;
whitebg([1,1,1]);set(gcf,'Color',[1,1,1]);
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 25)
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 25)
set(0,'defaultaxeslinewidth',2.5)
h = vol3d('CData',Change(:,[2:11],:)*100,'texture','3D');
view(3);  % 3D 뷰 설정
axis tight;  % 축 조정
xlim([0 10]);  % X축 범위
ylim([0 9]);  % Y축 범위
zlim([0 10]);  % Z축 범위
daspect([1 1 1]);  % 데이터 종횡비 설정
colormap(ax, jet);  % 색상 맵 설정
colorbar;  % 색상 막대 추가
% CLim 속성 설정
ax.CLim = [-100 100];
xlabel('\beta');
ylabel('\alpha');
zlabel('\gamma');

for i = 1:n
    for j = 1:n
        for k = 1:n
            Tratio2(i,j,k) = Tact2(i,j,k)/Treg2(i,j,k)          
        end
    end
end


figure;
ax = axes;
whitebg([1,1,1]);set(gcf,'Color',[1,1,1]);
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 25)
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 25)
set(0,'defaultaxeslinewidth',2.5)
h = vol3d('CData',Tratio2([2:10],[2:10],:),'texture','3D');
view(3);  % 3D 뷰 설정
axis tight;  % 축 조정
daspect([1 1 1]);  % 데이터 종횡비 설정
colormap(ax, jet);  % 색상 맵 설정
colorbar;  % 색상 막대 추가
% 투명도 맵 설정
alphamap('rampup');
alphamap('decrease', 0.1);  % 전체 투명도 감소
% CLim 속성 설정
ax.CLim = [0 0.1];
xlabel('\beta');
ylabel('\alpha');
zlabel('\gamma');

   