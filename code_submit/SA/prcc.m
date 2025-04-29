clear all; close all;

load rho1.mat
load pvalue1.mat
load lhsinfo.mat
load interestNames1.mat

rho_new(1:5,:) = rho(1:5,:);
rho_new(6:10,:) = rho(11:15,:);
pvalue_new(1:5,:) = pvalue(1:5,:);
pvalue_new(6:10,:) = pvalue(6:10,:);
interestNames_new(:,1) = interestNames(:,1);
interestNames_new(:,2) = interestNames(:,3);

load rho2.mat
load pvalue2.mat
load interestNames2.mat

rho_new(11:15,:) = rho(6:10,:);
pvalue_new(11:15,:) = pvalue(6:10,:);
interestNames_new(:,3) = interestNames(:,3);

hf = figure();
plotcc(rho_new, pvalue_new, lhsInfo(:, 1), interestNames_new);