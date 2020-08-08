%%
% This code runs the code associated figures 4 and 5 of Aponte et al., 2020
%
% Calls the function A1_freq_sweep_fn, which requires at minimum three inputs
%   1) function handle to parameter function (e.g., @A1_params_200722)
%   2) wave speed (unites of octaves/ms)
%   3) parameter set (see below for avaliable options)
%
% Runs both an upward and downward sweep
%
% The numerical tolerance of ode45, param.dx, and the timestep in
% param.tspan have been adjusted to optimize speed and size of files
%   These values can be adjusted directly in the A1_freq_sweep_fn.m and
%   A1_params_200722.m files
%   The upward and downward sweeps can also be run in parallel if desired
%
% Run create_fig4 and create_fig5 to produce the corresponding figures
%
% Codes completes on the order of ~1 minute on my laptop
%   Each segment of code is completely self-contained and can be run
%   individually
%
% This code was writen by Gregory Handy (2020)
% Please email ghandy@uchicago.edu with any questions
%%
clear; clc; close all;
tic;
%% Set up the direction and load the appropriate dataset

% For default figures
param_fn = @A1_params_200722;
dir_name = './A1_Sim_Data_Default/';
wave_speed = 20*10^(-3);
loop_wave_speeds = 1;

% For Supplemental figures
% param_fn = @A1_params_200722;
% dir_name = './A1_Sim_Data_SuppFig5/';
% wave_speed = 10*10^(-3);
% loop_wave_speeds = 0;

% For adjusted feedforward figures
% param_fn = @A1_params_200722_FFWD;
% dir_name = './A1_Sim_Data_FFWD/';
% wave_speed = 20*10^(-3);
% loop_wave_speeds = 0;

mkdir(dir_name)

%% Run the model with default parameters
fprintf('Running simulations for default parameters \n')
[param, delta_r_up, delta_r_down, DSI_values] = ...
    A1_freq_sweep_fn(param_fn, wave_speed,'Default'); %#ok<*ASGLU>
temp_filename = sprintf('freq_sweep_%.1f.mat',wave_speed*10^3);
filename = strcat(dir_name,temp_filename);
save(filename)

%% Run the model with a partial SOM block
% Two additional inputs: SOM_block_factor and SOM_inhibitory_curr
fprintf('Running simulations for partial SOM block parameters \n')
[param, delta_r_up, delta_r_down, DSI_values] = ...
    A1_freq_sweep_fn(param_fn, wave_speed,'SOM Block', 0.8, -1.5);

temp_filename = sprintf('freq_sweep_%.1f_SOM_Block.mat',wave_speed*10^3);
filename = strcat(dir_name,temp_filename);
save(filename,'param', 'delta_r_up', 'delta_r_down', 'DSI_values')

%% Run the model with a partial PV block
% Two additional inputs: PV_block_factor and PV_inhibitory_curr
fprintf('Running simulations for partial PV block parameters \n')
[ param, delta_r_up, delta_r_down, DSI_values] = ...
    A1_freq_sweep_fn(param_fn, wave_speed,'PV Block', 0.8, -1.5);
temp_filename = sprintf('freq_sweep_%.1f_PV_Block.mat',wave_speed*10^3);
filename = strcat(dir_name,temp_filename);
save(filename,'param', 'delta_r_up', 'delta_r_down', 'DSI_values')
% 
%% Run the model with different SOM projections
fprintf('Running simulations for shorter SOM projections \n')
[ param, delta_r_up, delta_r_down, DSI_values] = ...
    A1_freq_sweep_fn(param_fn, wave_speed,'SOM Projections',1);

temp_filename = sprintf('freq_sweep_%.1f_SOM_Projections_%.1f.mat',wave_speed*10^3,1);
filename = strcat(dir_name,temp_filename);
save(filename,'param', 'delta_r_up', 'delta_r_down', 'DSI_values')

[ param, delta_r_up, delta_r_down, DSI_values] = ...
    A1_freq_sweep_fn(param_fn, wave_speed,'SOM Projections',0.5);
temp_filename = sprintf('freq_sweep_%.1f_SOM_Projections_%.1f.mat',wave_speed*10^3,0.5);
filename = strcat(dir_name,temp_filename);
save(filename,'param', 'delta_r_up', 'delta_r_down', 'DSI_values')

%% Run the model in the non-ISN state
fprintf('Running simulations for the non-ISN state \n')
[ param, delta_r_up, delta_r_down, DSI_values] = ...
    A1_freq_sweep_fn(param_fn, wave_speed,'nonISN', 0.5);
temp_filename = sprintf('freq_sweep_%.1f_nonISN_%.1f.mat',wave_speed*10^3,0.5);
filename = strcat(dir_name,temp_filename);
save(filename,'param', 'delta_r_up', 'delta_r_down', 'DSI_values')

[ param, delta_r_up, delta_r_down, DSI_values] = ...
    A1_freq_sweep_fn(param_fn, wave_speed,'nonISN', 0.2);
temp_filename = sprintf('freq_sweep_%.1f_nonISN_%.1f.mat',wave_speed*10^3,0.2);
filename = strcat(dir_name,temp_filename);
save(filename,'param', 'delta_r_up', 'delta_r_down', 'DSI_values')

%% Run the model with default parameters at different wave speeds

if loop_wave_speeds == 1
    
    fprintf('Running simulations for default parameters at different wave speeds\n')
    
    [param, delta_r_up, delta_r_down, DSI_values] = ...
        A1_freq_sweep_fn(param_fn, 2.5*10^(-3),'Default');
    temp_filename = sprintf('freq_sweep_%.1f.mat',2.5);
    filename = strcat(dir_name,temp_filename);
    save(filename,'param', 'delta_r_up', 'delta_r_down', 'DSI_values')
    
    [param, delta_r_up, delta_r_down, DSI_values] = ...
        A1_freq_sweep_fn(param_fn, 5*10^(-3),'Default'); %#ok<*ASGLU>
    temp_filename = sprintf('freq_sweep_%.1f.mat',5);
    filename = strcat(dir_name,temp_filename);
    save(filename,'param', 'delta_r_up', 'delta_r_down', 'DSI_values')
    
    [param, delta_r_up, delta_r_down, DSI_values] = ...
        A1_freq_sweep_fn(param_fn, 10*10^(-3),'Default'); %#ok<*ASGLU>
    temp_filename = sprintf('freq_sweep_%.1f.mat',10);
    filename = strcat(dir_name,temp_filename);
    save(filename,'param', 'delta_r_up', 'delta_r_down', 'DSI_values')
    
    [param, delta_r_up, delta_r_down, DSI_values] = ...
        A1_freq_sweep_fn(param_fn, 40*10^(-3),'Default'); %#ok<*ASGLU>
    temp_filename = sprintf('freq_sweep_%.1f.mat',40);
    filename = strcat(dir_name,temp_filename);
    save(filename,'param', 'delta_r_up', 'delta_r_down', 'DSI_values')
    
    [param, delta_r_up, delta_r_down, DSI_values] = ...
        A1_freq_sweep_fn(param_fn, 80*10^(-3),'Default'); %#ok<*ASGLU>
    temp_filename = sprintf('freq_sweep_%.1f.mat',80);
    filename = strcat(dir_name,temp_filename);
    save(filename,'param', 'delta_r_up', 'delta_r_down', 'DSI_values')
    
    [param, delta_r_up, delta_r_down, DSI_values] = ...
        A1_freq_sweep_fn(param_fn, 160*10^(-3),'Default'); %#ok<*ASGLU>
    temp_filename = sprintf('freq_sweep_%.1f.mat',160);
    filename = strcat(dir_name,temp_filename);
    save(filename,'param', 'delta_r_up', 'delta_r_down', 'DSI_values')
end

fprintf('Simulations completed \n')
toc;