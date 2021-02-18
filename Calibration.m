% ---------------------------------------------- %
%% Two asset Portfolio Allocation Model %%
% Based on Liquid-Illiquid LCP code (from Benjamin Moll) %
% Calibration %
% Author: Lucas Rosso %
% Date: 18-02-2021 %
% Extremely Preliminar %
% ---------------------------------------------- %

clear all; close all; clc;

% parameters to calibrate: kappa (adjustment cost)
kappa_0 = 0.23; %guess

options = optimset('PlotFcns',@optimplotfval,'Display','iter');
optimal_params = fminsearchbnd(@SMM, kappa_0, 0, inf, options)

save('calibration','optimal_params')

% run MAIN
