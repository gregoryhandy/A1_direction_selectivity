
clear; clc; close all;

[ param, delta_r_up, delta_r_down, DSI_values] = A1_freq_sweep_fn(2.5*10^(-3),'Default');
save('./A1_Sim_Data/freq_sweep_2.5_Default.mat')

clear;
[ param, delta_r_up, delta_r_down, DSI_values] = A1_freq_sweep_fn(5*10^(-3),'Default');
save('./A1_Sim_Data/freq_sweep_5_Default.mat')

clear;
[ param, delta_r_up, delta_r_down, DSI_values] = A1_freq_sweep_fn(10*10^(-3),'Default');
save('./A1_Sim_Data/freq_sweep_10_Default.mat')

clear;
[ param, delta_r_up, delta_r_down, DSI_values] = A1_freq_sweep_fn(20*10^(-3),'Default');
save('./A1_Sim_Data/freq_sweep_20_Default.mat')

clear;
[ param, delta_r_up, delta_r_down, DSI_values] = A1_freq_sweep_fn(40*10^(-3),'Default');
save('./A1_Sim_Data/freq_sweep_40_Default.mat')

clear;
[ param, delta_r_up, delta_r_down, DSI_values] = A1_freq_sweep_fn(80*10^(-3),'Default');
save('./A1_Sim_Data/freq_sweep_80_Default.mat')

clear;
[ param, delta_r_up, delta_r_down, DSI_values] = A1_freq_sweep_fn(160*10^(-3),'Default');
save('./A1_Sim_Data/freq_sweep_160_Default.mat')

