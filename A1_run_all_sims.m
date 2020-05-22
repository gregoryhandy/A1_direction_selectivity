%%
% This code runs the code associated figures 4 and 5 of Aponte et al., 2020
%
% Calls the function A1_freq_sweep_fn, which requires at minimum two inputs
%   1) wave speed (unites of octaves/ms)
%   2) parameter set (see below for avaliable options)
%
% Runs both an upward and downward sweep
%
% The numerical tolerance of ode45, param.dx, and the timestep in
% param.tspan has been adjusted to optimize speed and size fo files, but
% results in numerical artifacts in DSI curves
%   These values can be adjusted directly in the A1_freq_sweep_fn.m and
%   A1_params_Kato_200211.m files
%   The upward and downward sweeps can also be run in parallel if desired
%
% Run create_fig4 and create_fig5 to produce the corresponding figures
%
% Codes completes on the order of ~2 minutes on my laptop
%   Each segment of code is completely self-contained and can be run
%   individually
%
% This code was writen by Gregory Handy (2020)
% Please email gregoryhandy@pitt.edu with any questions
%%
clear; clc; close all;

%% Run the model with default parameters
fprintf('Running simulations for default parameters \n')
[param, delta_r_up, delta_r_down, DSI_values] = ...
    A1_freq_sweep_fn(20*10^(-3),'Default'); %#ok<*ASGLU>
filename = sprintf('./A1_Sim_Data/freq_sweep_%.f_default.mat',20);
save(filename)

%% Run the model with a partial SOM block
% Two additional inputs: SOM_block_factor and SOM_inhibitory_curr
fprintf('Running simulations for partial SOM block parameters \n')
[param, delta_r_up, delta_r_down, DSI_values] = ...
    A1_freq_sweep_fn(20*10^(-3),'SOM Block', 0.8, -1.5);

filename = sprintf('./A1_Sim_Data/freq_sweep_%.f_SOM_Block.mat',20);
save(filename,'param', 'delta_r_up', 'delta_r_down', 'DSI_values')

%% Run the model with a partial PV block
% Two additional inputs: PV_block_factor and PV_inhibitory_curr
fprintf('Running simulations for partial PV block parameters \n')
[ param, delta_r_up, delta_r_down, DSI_values] = ...
    A1_freq_sweep_fn(20*10^(-3),'PV Block', 0.8, -1.5);
filename = sprintf('./A1_Sim_Data/freq_sweep_%.f_PV_Block.mat',20);
save(filename,'param', 'delta_r_up', 'delta_r_down', 'DSI_values')

%% Run the model with different SOM projections
fprintf('Running simulations for shorter SOM projections \n')
[ param, delta_r_up, delta_r_down, DSI_values] = ...
    A1_freq_sweep_fn(20*10^(-3),'SOM Projections',1);
filename = sprintf('./A1_Sim_Data/freq_sweep_%.f_SOM_Projections_%.1f.mat',20,1);
save(filename,'param', 'delta_r_up', 'delta_r_down', 'DSI_values')

[ param, delta_r_up, delta_r_down, DSI_values] = ...
    A1_freq_sweep_fn(20*10^(-3),'SOM Projections',0.5);
filename = sprintf('./A1_Sim_Data/freq_sweep_%.f_SOM_Projections_%.1f.mat',20,0.5);
save(filename,'param', 'delta_r_up', 'delta_r_down', 'DSI_values')

%% Run the model in the non-ISN state
fprintf('Running simulations for the non-ISN state \n')
[ param, delta_r_up, delta_r_down, DSI_values] = ...
    A1_freq_sweep_fn(20*10^(-3),'nonISN', 0.5);
filename = sprintf('./A1_Sim_Data/freq_sweep_%.f_nonISN_%.1f.mat',20,0.5);
save(filename,'param', 'delta_r_up', 'delta_r_down', 'DSI_values')

[ param, delta_r_up, delta_r_down, DSI_values] = ...
    A1_freq_sweep_fn(20*10^(-3),'nonISN', 0.2);
filename = sprintf('./A1_Sim_Data/freq_sweep_%.f_nonISN_%.1f.mat',20,0.2);
save(filename,'param', 'delta_r_up', 'delta_r_down', 'DSI_values')

%% Run the model with default parameters at different wave speeds

fprintf('Running simulations for default parameters at different wave speeds\n')

[ param, delta_r_up, delta_r_down, DSI_values] = ...
    A1_freq_sweep_fn(2.5*10^(-3),'Default');
filename = sprintf('./A1_Sim_Data/freq_sweep_%.1f_default.mat',2.5);
save(filename,'param', 'delta_r_up', 'delta_r_down', 'DSI_values')

[ param, delta_r_up, delta_r_down, DSI_values] = ...
    A1_freq_sweep_fn(5*10^(-3),'Default');
filename = sprintf('./A1_Sim_Data/freq_sweep_%.f_default.mat',5);
save(filename,'param', 'delta_r_up', 'delta_r_down', 'DSI_values')

[ param, delta_r_up, delta_r_down, DSI_values] = ...
    A1_freq_sweep_fn(10*10^(-3),'Default');
filename = sprintf('./A1_Sim_Data/freq_sweep_%.f_default.mat',10);
save(filename,'param', 'delta_r_up', 'delta_r_down', 'DSI_values')

[ param, delta_r_up, delta_r_down, DSI_values] = ...
    A1_freq_sweep_fn(40*10^(-3),'Default');
filename = sprintf('./A1_Sim_Data/freq_sweep_%.f_default.mat',40);
save(filename,'param', 'delta_r_up', 'delta_r_down', 'DSI_values')

[ param, delta_r_up, delta_r_down, DSI_values] = ...
    A1_freq_sweep_fn(80*10^(-3),'Default');
filename = sprintf('./A1_Sim_Data/freq_sweep_%.f_default.mat',80);
save(filename,'param', 'delta_r_up', 'delta_r_down', 'DSI_values')

[ param, delta_r_up, delta_r_down, DSI_values] = ...
    A1_freq_sweep_fn(160*10^(-3),'Default');
filename = sprintf('./A1_Sim_Data/freq_sweep_%.f_default.mat',160);
save(filename,'param', 'delta_r_up', 'delta_r_down', 'DSI_values')

fprintf('Simulations completed \n')