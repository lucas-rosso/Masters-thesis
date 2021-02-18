% ---------------------------------------------- %
%% Two asset Portfolio Allocation Model %%
% Based on Liquid-Illiquid LCP code (from Benjamin Moll) %
% Calibration %
% Author: Lucas Rosso %
% Date: 17-02-2021 %
% Extremely Preliminar %
% ---------------------------------------------- %

clear all; close all; clc;

% parameters to calibrate: kappa (adjustment cost) and rho (discount rate)
kappa_0 = 0.23; 
bmin_0   = -1; 

params = [kappa_0, bmin_0];

% loss = SMM(params)

LB = [0, -inf];
UB = [inf,inf];

options = optimset('PlotFcns',@optimplotfval,'Display','iter');
optimal_params = fminsearchbnd(@SMM, params, LB, UB, options)

save('calibration','optimal_params')

run Main_code_05022021
