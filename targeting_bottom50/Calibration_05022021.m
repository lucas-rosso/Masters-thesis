% ---------------------------------------------- %
%% Two asset Portfolio Allocation Model %%
% Based on Liquid-Illiquid LCP code (from Benjamin Moll) %
% Calibration %
% Author: Lucas Rosso %
% Date: 05-02-2021 %
% Extremely Preliminar %
% ---------------------------------------------- %

clear all; close all; clc;

% parameters to calibrate: kappa (adjustment cost) and rho (discount rate)
kappa_0 = 0.21; 
% rho_0   = 0.04; 

% params = [kappa_0, rho_0];

% loss = SMM(params)

% LB = [0, 0];
% UB = [inf,inf];

options = optimset('PlotFcns',@optimplotfval,'Display','iter');
optimal_params = fminsearchbnd(@SMM, kappa_0, 0, inf, options)

save('calibration','optimal_params')

% run Main_code_05022021
