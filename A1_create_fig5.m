%%
% Recreate Figure 5 using the data produced from A1_run_all_sims.m 
%
% Note: some numerical artifacts in DSI curves might appear due to
% adjustments to the code to optimize speed and sizes of data
%
% This code was writen by Gregory Handy (2020)
% Please email ghandy@uchicago.edu with any questions
%%
clear; close all; clc;

restoredefaultpath;
folder = fileparts(which('A1_create_fig5.m')); 
addpath(genpath(folder));
rmpath(folder)

fig_handle = figure(66);
total_rows = 4;

%% Choose the name of the data directory and wave speed

% Default
data_dir = './A1_Sim_Data_Default';
wave_speed = 20*10^-3;

% Supplemental
% data_dir = './A1_Sim_Data_SuppFig5';
% wave_speed = 10*10^-3;

% Feedforward
% data_dir = './A1_Sim_Data_FFWD';
% wave_speed = 20*10^-3;


%% Load the SOM data

temp_file_name = sprintf('/freq_sweep_%.1f.mat',wave_speed*10^3);
load(strcat(data_dir,temp_file_name));
DSI_default = DSI_values;

temp_file_name = sprintf('/freq_sweep_%.1f_SOM_Block.mat',wave_speed*10^3);
load(strcat(data_dir,temp_file_name));
DSI_param_set_2 = DSI_values;
delta_r_up_param_set_2 = delta_r_up;
delta_r_down_param_set_2 = delta_r_down;

num_data_sets = 2;

row_num = 1; 
create_fig5_row(num_data_sets,DSI_default, [], DSI_param_set_2,...
    delta_r_up_param_set_2,delta_r_down_param_set_2, param, total_rows, ...
    row_num, fig_handle)

%% Load the PV data

temp_file_name = sprintf('/freq_sweep_%.1f_PV_Block.mat',wave_speed*10^3);
load(strcat(data_dir,temp_file_name));
DSI_param_set_2 = DSI_values;
delta_r_up_param_set_2 = delta_r_up;
delta_r_down_param_set_2 = delta_r_down;

num_data_sets = 2;

row_num = 2;
create_fig5_row(num_data_sets,DSI_default, [], DSI_param_set_2,...
    delta_r_up_param_set_2,delta_r_down_param_set_2, param, total_rows, ...
    row_num, fig_handle)

%% Load the SOM projections data

temp_file_name = sprintf('/freq_sweep_%.1f_SOM_Projections_1.0.mat',wave_speed*10^3);
load(strcat(data_dir,temp_file_name));
DSI_param_set_1 = DSI_values;

temp_file_name = sprintf('/freq_sweep_%.1f_SOM_Projections_0.5.mat',wave_speed*10^3);
load(strcat(data_dir,temp_file_name));
DSI_param_set_2 = DSI_values;
delta_r_up_param_set_2 = delta_r_up;
delta_r_down_param_set_2 = delta_r_down;

num_data_sets = 3;

row_num = 3;
create_fig5_row(num_data_sets, DSI_default, DSI_param_set_1, DSI_param_set_2,...
    delta_r_up_param_set_2,delta_r_down_param_set_2, param, total_rows, ...
    row_num, fig_handle)

%% Load the nonISN data

temp_file_name = sprintf('/freq_sweep_%.1f_nonISN_0.5.mat',wave_speed*10^3);
load(strcat(data_dir,temp_file_name));
DSI_param_set_1 = DSI_values;

temp_file_name = sprintf('/freq_sweep_%.1f_nonISN_0.2.mat',wave_speed*10^3);
load(strcat(data_dir,temp_file_name));
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
legend box off
xticks([]);
xlabel('')

subplot(total_rows,3,2+3*(row_num-1))
xticks([]);
xlabel('')
legend box off

subplot(total_rows,3,3+3*(row_num-1))
xticks([]);
xlabel('')
legend('Upward','Downward')
legend box off

%%
row_num = 2;
subplot(total_rows,3,1+3*(row_num-1))
legend('Default','PV Partial Block')
legend box off
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
legend('\lambda = 2 (Default)','\lambda = 1','\lambda = 0.5 (narrow)')
legend box off
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
legend('\beta =1.0 (Default, ISN)','\beta = 0.5','\beta = 0.2 (non-ISN)')
legend box off
subplot(total_rows,3,2+3*(row_num-1))

% To export to Adobe illustrator, use this:
% print -painters -depsc Fig5.eps