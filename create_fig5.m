%%
% Recreate Figure 5 using the data produced from A1_run_all_sims.m 
%
% Note: some numerical artifacts in DSI curves might appear due to
% adjustments to the code to optimize speed and sizes of data
%
% This code was writen by Gregory Handy (2020)
% Please email gregoryhandy@pitt.edu with any questions
%%
clear; close all; clc;

restoredefaultpath;
folder = fileparts(which('create_fig5.m')); 
addpath(genpath(folder));
rmpath(folder)

fig_handle = figure(66);
total_rows = 4;

%% Load the SOM data

load('./A1_Sim_Data/freq_sweep_20_default.mat');
DSI_default = DSI_values;

load('./A1_Sim_Data/freq_sweep_20_SOM_Block.mat');
DSI_param_set_2 = DSI_values;
delta_r_up_param_set_2 = delta_r_up;
delta_r_down_param_set_2 = delta_r_down;

num_data_sets = 2;

row_num = 1; 
create_fig5_row(num_data_sets,DSI_default, [], DSI_param_set_2,...
    delta_r_up_param_set_2,delta_r_down_param_set_2, param, total_rows, ...
    row_num, fig_handle)

%% Load the PV data


load('./A1_Sim_Data/freq_sweep_20_PV_Block.mat')
DSI_param_set_2 = DSI_values;
delta_r_up_param_set_2 = delta_r_up;
delta_r_down_param_set_2 = delta_r_down;

num_data_sets = 2;

row_num = 2;
create_fig5_row(num_data_sets,DSI_default, [], DSI_param_set_2,...
    delta_r_up_param_set_2,delta_r_down_param_set_2, param, total_rows, ...
    row_num, fig_handle)

%% Load the SOM projections data

load('./A1_Sim_Data/freq_sweep_20_SOM_Projections_1.0.mat')
DSI_param_set_1 = DSI_values;

load('./A1_Sim_Data/freq_sweep_20_SOM_Projections_0.5.mat')
DSI_param_set_2 = DSI_values;
delta_r_up_param_set_2 = delta_r_up;
delta_r_down_param_set_2 = delta_r_down;

num_data_sets = 3;

row_num = 3;
create_fig5_row(num_data_sets, DSI_default, DSI_param_set_1, DSI_param_set_2,...
    delta_r_up_param_set_2,delta_r_down_param_set_2, param, total_rows, ...
    row_num, fig_handle)

%% Load the nonISN data

load('./A1_Sim_Data/freq_sweep_20_nonISN_0.5.mat')
DSI_param_set_1 = DSI_values;

load('./A1_Sim_Data/freq_sweep_20_nonISN_0.2.mat')
DSI_param_set_2 = DSI_values;
delta_r_up_param_set_2 = delta_r_up;
delta_r_down_param_set_2 = delta_r_down;

num_data_sets = 3;

row_num = 4;
create_fig5_row(num_data_sets, DSI_default, DSI_param_set_1, DSI_param_set_2,...
    delta_r_up_param_set_2,delta_r_down_param_set_2, param, total_rows, ...
    row_num, fig_handle)

%% Make row-specific adjustments here

row_num = 1;
subplot(total_rows,3,1+3*(row_num-1))
legend('Default','SOM Partial Block')
xticks([]);
xlabel('')

subplot(total_rows,3,2+3*(row_num-1))
xticks([]);
xlabel('')
legend('Upward Sweep','Downward Sweep')
title('BF = 5 (kHz)')

subplot(total_rows,3,3+3*(row_num-1))
xticks([]);
xlabel('')
legend('Upward','Downward')

%%
row_num = 2;
subplot(total_rows,3,1+3*(row_num-1))
legend('Default','PV Partial Block')
xticks([]);
xlabel('')

subplot(total_rows,3,2+3*(row_num-1))
xticks([]);
xlabel('')

subplot(total_rows,3,3+3*(row_num-1))
xticks([]);
xlabel('')

%%
row_num = 3;
subplot(total_rows,3,1+3*(row_num-1))
legend('Default','\lambda = 1','\lambda = 0.5')
xticks([]);
xlabel('')

subplot(total_rows,3,2+3*(row_num-1))
xticks([]);
xlabel('')

subplot(total_rows,3,3+3*(row_num-1))
xticks([]);
xlabel('')

%%

row_num = 4;
subplot(total_rows,3,1+3*(row_num-1))
legend('Default','\beta = 0.5','\beta = 0.2')

% To export to Adobe illustrator, use this:
% print -painters -depsc Figure_6.eps
