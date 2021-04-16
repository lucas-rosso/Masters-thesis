% ---------------------------------------------- %
%% Two asset Portfolio Allocation Model %%
% Calibration code including Type Dependence
% Author: Lucas Rosso %
% Date: 15-04-2021 %
% Extremely Preliminar %
% ---------------------------------------------- %

clear all; close all; clc;

% parameters to calibrate: kappa (adjustment cost)
kappa_0 = 0.2323; %guess
varthe_0 = 0.02; 
sigma_0 = 0.25;

params = [kappa_0, varthe_0,sigma_0];
LB = [0,0,0]; % bounds
UB = [inf,inf,inf];

options = optimset('PlotFcns',@optimplotfval,'Display','iter');
optimal_params = fminsearchbnd(@SMM, params, LB, UB, options)

save('calibration','optimal_params')

% run MAIN
